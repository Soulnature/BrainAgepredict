import torch
import torch.nn as nn
import torchvision.models as models
from torch.optim import Adam
import pandas as pd
import utiles
from tqdm import tqdm
import os
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from torch.autograd import Variable
from torch.utils.data import DataLoader
from SACN import Age_regression,Encoder,Age_mae
from itertools import chain
from datetime import datetime
from  utiles import image_domian_gender,my_KLDivLoss,num2vect
### first load the encoder and fine tune the Age regression and age mae ###
# 加载预训练模型
torch.multiprocessing.set_sharing_strategy('file_system')

class FineTune:
    def __init__(self, encoder_path, age_path,train_dir, num_epochs,bin_num,best_model_path,retrain):
        self.encoder_path = encoder_path
        self.train_dir = train_dir
        self.num_epochs = num_epochs
        self.bin_num=bin_num
        self.age_path=age_path
        self.min_mae=10000
        self.r2=-1000
        self.retrain=retrain
        self.best_model_path = best_model_path
        # Load the encoder and fine tune the Age regression and age mae
        self.Encoder_m = Encoder()
        self.encoder = torch.nn.DataParallel(self.Encoder_m, device_ids=[0])
        self.encoder.load_state_dict(torch.load(self.encoder_path))
        self.Age_regression = Age_regression(output_dim=90)
        self.Age_regression = torch.nn.DataParallel(self.Age_regression, device_ids=[0])
        self.Age_regression.load_state_dict(torch.load(self.age_path))
        self.optimizer = Adam(chain(self.Age_regression.parameters(), self.encoder.parameters()), lr=0.01)
        self.prepare_data()
        if retrain:
            self.fine_tune()

    def prepare_data(self):
        # Load the data
        self.data = pd.read_csv(self.train_dir)
        self.train_data = self.data[0:2400]
        self.val_data = self.data[2400:2700]
        self.test_data=self.data[2700:self.data.shape[0]]
        self.train_dataset = image_domian_gender(self.train_data['Loc'].values,  self.train_data['Gender'].values,  self.train_data['Scanner'].values,
                                                  self.train_data['Age'].values)
        self.val_data = image_domian_gender(self.val_data['Loc'].values,  self.val_data['Gender'].values,  self.test_data['Scanner'].values,
                                                  self.val_data['Age'].values)
        self.test_data = image_domian_gender(self.test_data['Loc'].values,  self.test_data['Gender'].values,  self.test_data['Scanner'].values,
                                                  self.test_data['Age'].values)
        self.train_dataloader = DataLoader(self.train_dataset, batch_size=12, num_workers=4)
        self.val_dataloader = DataLoader(self.val_data, batch_size=1, num_workers=1)
        self.test_dataloader = DataLoader(self.test_data, batch_size=1, num_workers=1)
    def load_new_data(self,new_loc):
        new_data = pd.read_csv(new_loc)
        new_data_dataset = image_domian_gender(new_data['Loc'].values,  new_data['Gender'].values,  new_data['Scanner'].values,
                                                  new_data['Age'].values)
        return(DataLoader(new_data_dataset, batch_size=1, num_workers=4))
        
    def save_model(self):
            # Save the best model
            date_str = datetime.now().strftime('%Y%m%d')
            self.model_path = f"{self.best_model_path}_{date_str}"
            torch.save({
                'encoder': self.encoder.state_dict(),
                'Age_regression': self.Age_regression.state_dict(),
            }, self.model_path)
    
    def train(self, use_gpu):
        # setup models
        self.Age_regression.train()
        self.encoder.train()
        bin_step = 1
        sigma = 1
        maeloss = nn.L1Loss(reduction='mean')
        loss_all=0
        mae=0
        for batch_idx, data in enumerate(self.train_dataloader):
            # setup hyperparameters
            # prepare the data #
            input_data, label_all = data
            gender = label_all[0]
            domain = label_all[1]
            age_or = label_all[2].float()
            ##encoder the age label and domain label
            age, bin1 = num2vect(age_or, self.bin_num, bin_step, sigma)
            age = torch.tensor(age, dtype=torch.float32)
            ## domian encoder #
            if use_gpu:
                input_data, age, age_or, gender, domain = Variable(input_data.cuda()), Variable(age.cuda()), Variable(
                    age_or.cuda()), Variable(gender.cuda()), Variable(domain.cuda())
            self.optimizer.zero_grad()
            # compute the output of source domain and target domain
            original_feature = self.encoder(input_data)
            # compute the class loss of age
            age_preds,age_one = self.Age_regression(original_feature)
            age_preds = age_preds.reshape(age_one.size()[0],90)
            class_loss = my_KLDivLoss(age_preds, age) 
            # mae lass
            mae_los = maeloss(age_one[:, 0], age_or)
            ##
            loss =mae_los+class_loss
            loss_all += loss
            mae+=mae_los
            loss.backward()
            self.optimizer.step()
            if (batch_idx + 1) % 10 == 0: 
                 print(
                    '[{}/{} ({:.0f}%)]\t allLoss: {:.6f}\t ageClass Loss: {:.6f} \t MAE Loss: {:.6f}'.format(
                        batch_idx * len(input_data), len(self.train_dataloader.dataset),
                        100. * batch_idx / len(self.train_dataloader), loss.item(), class_loss.item(),
                        mae_los.item()
                        ))
        ##save the best model #
        return  float(mae) / len(self.train_dataloader), float(loss_all)/len(self.train_dataloader)
    def test(self, use_gpu):
        self.load_model()
        test_res,domain,_=self.val_model(True,self.test_dataloader)
        return test_res

    def load_model(self):
    # Load the best model
        checkpoint = torch.load(self.model_path)
        self.encoder.load_state_dict(checkpoint['encoder'])
        self.Age_regression.load_state_dict(checkpoint['Age_regression'])
        
    def val_model(self, use_gpu,datamoudle):
        MAE = 0.0
        self.Age_regression.eval()
        self.encoder.eval()
        correct_domain = 0.0
        bin_step = 1
        sigma = 1
        pred_all_age = []
        original_age = []
        maeloss = nn.L1Loss(reduction='mean')
        with torch.no_grad():
            for batch_idx, data in enumerate(tqdm(datamoudle)):
                input_data, label_all = data
                gender = label_all[0]
                domain = label_all[1]
                age_or = label_all[2]
                age, bin1 = num2vect(age_or, self.bin_num, bin_step, sigma)
                age = torch.tensor(age, dtype=torch.float32)
                gender = torch.tensor(gender, dtype=torch.int64)
                domain = torch.tensor(domain, dtype=torch.int64)
                if use_gpu:
                    input_data, age, gender, domain = Variable(input_data.cuda()), Variable(age.cuda()), Variable(
                        gender.cuda()), Variable(domain.cuda())
                else:
                    input_data, age, gender, domain = Variable(input_data), Variable(age), Variable(gender), Variable(
                        domain)
                original_feature = self.encoder(input_data)
                age_preds,age_one = self.Age_regression(original_feature)
                age_preds = age_preds.reshape(age_one.size()[0],90)
                output_data2 = age_one[0].detach().cpu().numpy().reshape([1, -1])
                mae_data_re = np.abs(output_data2 - age_or.numpy())
                out_data = age_preds[0].detach().cpu().numpy().reshape([1, -1])
                prob = np.exp(out_data)
                pred = prob @ bin1
                original_age.append(age_or.numpy())
                mae_data = np.abs(pred - age_or.numpy())
                pred_age_nu = np.asarray(((output_data2 + pred) / 2))
                MAE += (mae_data+mae_data_re[0])/2
                pred_all_age.append(pred_age_nu)
        pred_age_all = np.asarray(pred_all_age)[:, 0]
        pred_age_all = pred_age_all.astype(np.float)
        original_age = np.asarray(original_age)[:, 0]
        original_age = original_age.astype(np.float)
        original_age = original_age.astype(float)
        pred_age_all = pred_age_all.astype(float)  
        mae=np.mean(np.abs(pred_age_all[:,0]-original_age))
        r_value = pearsonr(original_age, pred_age_all)
        r_square = r2_score(original_age, pred_age_all)
        if(r_square>self.r2):
            self.r2=r_square
            self.min_mae=np.mean(float(MAE) / len(datamoudle.dataset))
        print(r_value)
        print(r_square)
        print(mae)
        print('\n Brain Age  Accuracy: {} ({:.4f})'
              'Domain Accuracy: {}/{} ({:.4f}%)\n'.
              format(
                  MAE, float(MAE) / len(datamoudle.dataset),
                  correct_domain, len(datamoudle.dataset),
                  100. * float(correct_domain) / len(datamoudle.dataset)
              ))
        return mae, float(correct_domain) / len(datamoudle.dataset),pred_age_all

    def predict_data(self,path_model,other,inden_loc):
        checkpoint = torch.load(path_model)
        self.encoder.load_state_dict(checkpoint['encoder'])
        self.Age_regression.load_state_dict(checkpoint['Age_regression'])      
        ###predicted the brain age of the data ##
        if other:
            indenp_data=self.load_new_data(inden_loc)
            test_res,domain,pred_age=self.val_model(True,indenp_data)
        else:
            test_res,domain,pred_age=self.val_model(True,self.test_dataloader)
        return pred_age
        ###
    def fine_tune(self):
        best_val_metric = float('inf')  # or -float('inf') depending on whether lower is better or higher is better
        for epoch in range(self.num_epochs):
            print(epoch)
            self.train(True)  # Assuming use_gpu=True and bin_range=range(0, 90, 1)
            val_metric,_,_ = self.val_model(True,self.val_dataloader)  # Assuming use_gpu=True and bin_range=range(0, 90, 1)
            if val_metric < best_val_metric:  # or > for higher is better
                best_val_metric = val_metric
                self.save_model()
                print('save model was finished')
        #self.test(True)  # Assuming use_gpu=True and bin_range=range(0, 90, 1)
        print(self.min_mae)
        print(self.r2)

encoder_path='/home1/zhaoxz/reCOns/save_Models/2023-11-02SACN1SACN_encoder'
os.environ["CUDA_VISIBLE_DEVICES"] = '5'
age_path='/home1/zhaoxz/reCOns/save_Models/2023-11-02SACN1SACNage_predictor'
UKB='/home1/zhaoxz/reCOns/UKB.csv'
FT=FineTune(encoder_path,age_path,'/home1/zhaoxz/reCOns/UKB_match.csv',150,[5,95],'/home1/zhaoxz/reCOns/UKB_fine_tune/',True)
FT=FineTune(encoder_path,age_path,UKB,150,[5,95],'/home1/zhaoxz/reCOns/UKB_fine_tune/',False)
## AD path ##

pred_ahe=FT.predict_data('/home1/zhaoxz/reCOns/UKB_fine_tune/_20231129',True,UKB)

#np.savetxt('PNC_pred.txt',pred_ahe[:,0])
#print(pred_ahe)