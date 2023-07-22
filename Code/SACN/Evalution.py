## evaluate the res ## 300 epoch of constant #
import torch
import pandas as pd
from   utiles import image_domian_gender
from utiles import my_KLDivLoss
from model import Encoder, Domain_classifier, Age_regression,Age_mae
#from Resnet import resnet50, Domain_classifier,Age_regression,Age_mae,Gende_net
from utiles import num2vect
import numpy as np
from tqdm import tqdm
import os
from scipy.stats import pearsonr
from torch.autograd import Variable
from torch.utils.data import DataLoader
# compute the metrics for predicted age
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model  import LinearRegression
import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system')
class Metrics_age(object):
    def __init__(self,age,predict_age,correction=True):
        self.age=age
        self.predict_age=predict_age
        self.correction=correction
    def checkdata(self):
        age=np.asarray(self.age)
        predict_age=np.asarray(self.predict_age)
        if len(age.shape)>1:
            age=age[0:]
        if len(predict_age.shape)>1:
            predict_age=predict_age[0,:]
        return age,predict_age
    def correction_age(self):
        age_real,predict_age_real=self.checkdata()
        age_realt=age_real.reshape(len(age_real),1)
        model_re=LinearRegression()
        model_re.fit(age_realt,predict_age_real)
        alph=model_re.coef_[0]
        beta=model_re.intercept_
        cor_age=np.zeros(len(age_real))
        for i in range(len(age_real)):
            cor_age[i]=predict_age_real[i]+(age_real[i]-(alph*age_real[i]+beta))
        return cor_age,age_real
    def metrics_com(self):
        if self.correction==True:
            cor_age,age_real=self.correction_age()
        else:
            age_real,cor_age=self.checkdata()
        MAE=mean_absolute_error(cor_age,age_real)
        R_square=r2_score(cor_age,age_real)
        r=pearsonr(cor_age,age_real)
        return MAE,R_square,r,cor_age
class evaluate_data:
    def __init__(self,fea_model,reg_model,age_cla_model,dataloader,use_gpu,save_encoder_st):
        self.fea=fea_model
        self.reg=reg_model
        self.age_c=age_cla_model
        self.dataloader=dataloader
        self.use_gpu=use_gpu
        self.save_encoder_st=save_encoder_st
    def __call__(self):
        # fucntion: load model and  set the model on CUDA
        #encoder = Encoder()
        encoder = Encoder()
        encoder = torch.nn.DataParallel(encoder, device_ids=[0])
        age_pred = torch.nn.DataParallel(Age_regression(output_dim=90), device_ids=[0])
        age_mae = torch.nn.DataParallel(Age_mae(), device_ids=[0])
        age_mae.load_state_dict(torch.load(self.reg))
        encoder.load_state_dict(torch.load(self.fea))
        age_pred.load_state_dict(torch.load(self.age_c))
        ##set the model run on cuda
        encoder.cuda()
        age_pred.cuda()
        age_mae.cuda()
        #set the model model as the evaluation mode #
        encoder.eval()
        age_pred.eval()
        age_mae.eval()
        allage_value=[]
        bin_step = 1
        sigma = 1
        bin_range = [5, 95]
        use_gpu = True
        MAE = 0.0
        real_age=[]
        save_encoder=[]
        ### test the model #
       # model = gcam.inject(encoder, output_dir="./attention_maps/", backend="gcam", layer='conv_5.0', save_maps=True)
        with torch.no_grad():
            for batch_idx, data in enumerate(tqdm(self.dataloader)):
                # setup hyperparameters
                # prepare the data
                input_data, label_all = data
                gender = label_all[0]
                domain = label_all[1]
                # print(domain)
                age_or = label_all[2]
                real_age.append(age_or.numpy())
                # print(age_or)
                ##encoder the age label and domain label
                age, bin1 = num2vect(age_or.numpy().tolist(), bin_range, bin_step, sigma)
                age = torch.tensor(age, dtype=torch.float32)
                gender = torch.tensor(gender, dtype=torch.int64)
                domain = torch.tensor(domain, dtype=torch.int64)
                ## domian encoder #
                if use_gpu:
                    input_data, age, gender, domain = Variable(input_data.cuda()), Variable(age.cuda()), Variable(
                        gender.cuda()), Variable(domain.cuda())
                else:
                    input_data, age, gender, domain = Variable(input_data), Variable(age), Variable(gender), Variable(
                        domain)
                # setup optimizer
                if self.save_encoder_st:
                    encoder_data = encoder(input_data)
                    encoder_dim = encoder_data.detach().cpu().numpy().reshape([-1, 1])
                    save_encoder.append(encoder_dim)
                    output1 = age_pred(encoder(input_data))
                    output2 = age_mae(encoder(input_data))
                    # define the MAE of the data
                    out_data = output1[0].detach().cpu().numpy().reshape([1, -1])
                    output_data2 = output2[0].detach().cpu().numpy().reshape([1, -1])
                    prob = np.exp(out_data)
                    pred = prob @ bin1
                    # compute the MAE #
                    #print(pred)
                    mae_data = np.abs(pred - age_or.numpy())
                    mae_data_re = np.abs(output_data2 - age_or.numpy())
                    MAE += ((mae_data + mae_data_re) / 2)
                   # print(pred)
                    #age_two = np.asarray(output_data2)
                    age_two = np.asarray(((output_data2 + pred) / 2))
                    allage_value.append(age_two[:, 0])
               # allage_value

        # real_age=real_age.detach().cpu().numpy().reshape([1, -1])
        real_age=np.asarray(real_age)
        real_age=real_age.reshape((real_age.shape[0],1))
        #real_age = real_age.detach().cpu().numpy().reshape([1, -1])
        allage_value=np.asarray(allage_value)
     #   print(allage_value)
        conbind_df=np.concatenate((real_age,allage_value), axis=1)
        conbind_df = conbind_df[~np.isnan(conbind_df).any(axis=1), :]
        me_trix=Metrics_age(conbind_df[:,0],conbind_df[:,1],False)
        MAE,R_square, r,cor_age=me_trix.metrics_com()
       # np.save("PPMI_input_fea.npy",np.asarray(save_encoder))
      #  np.savetxt('ADNI_pred.txt',np.asarray(cor_age))
        #np.savetxt('original_pred.txt', np.asarray(cor_age))
        print(MAE)
        print(R_square)
        print(r)
def spilt_dataset(csv_object):
    train_cat=csv_object[0:1]
    val_cat=csv_object[0:1]
    scanner_index=np.unique(csv_object['Scanner'].values)
    for i in scanner_index:
        oneset=csv_object[(csv_object['Scanner']==scanner_index[i])]
        train_data,test_data=train_test_split(oneset,test_size=0.2,random_state=0)
        train_cat=pd.concat([train_cat,train_data],axis=0)
        val_cat=pd.concat([val_cat,test_data])
    val_cat=val_cat.drop([0])
    train_cat=train_cat.drop([0])
    return train_cat,val_cat
if __name__ == '__main__':
    # set gpu number for evaluted #
    ##
    os.environ["CUDA_VISIBLE_DEVICES"] = "3"
   # fea_model = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_fea__val_dann_1.24.pth'
   # age_clas = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_sfcn_all_val_dann_end_agecls_com_1.24.pth'
    #age_mae = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_sfcn_reg_val_dann_ukb_cons_end_mae_1.24.pth'
    #fea_model = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_fea__val_dann_cons_300.pth'
 #  age_clas = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_sfcn_all_val_dann_cons_300.pth'
   # age_mae = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_sfcn_reg_val_dann_cons_300.pth'

    fea_model = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_fea__val_dann_2.39.pth'
    age_clas = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_sfcn_all_val_dann_end_agecls_com_2.39.pth'
    age_mae = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_sfcn_reg_val_dann_ukb_cons_end_mae_2.39.pth'
     ## Resnet path ###
    #fea_model = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_fea_resnet50_domain.pth'
    #age_clas = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_resnet50_class_domain.pth'
    #age_mae = '/data/home/zhaoxz/Documents/domain_adaption/Models/best_model_reg_resnet50_domain.pth'
    #root_dir='/data/home/zhaoxz/UKB_match.csv'
   # root_dir = '/data/home/zhaoxz/ABIDE/real_ABIDE.csv'
   # root_dir='/data/home/zhaoxz/ABCD_MNI.csv'
   # root_dir = '/data/home/zhaoxz/ABIDE_data.csv'
    root_dir='/data/home/zhaoxz/Documents/domain_adaption/PNC.csv'
    data = pd.read_csv(root_dir)
   # data=data[6000:6500]
   # data=data.loc[data["Group"] ==2]
    # file_info,test_set=spilt_dataset(data)
    # train_d, test_d = train_test_split(file_info, test_size=0.1,random_state=0)
    #data['Age']= data['Age'].values
    train_dataset = image_domian_gender(data['Loc'].values, data['Gende'].values, data['Scanner'].values,
                                        data['Age'].values)
    # val_dataset = image_domian_gender(test_d['Loc'].values, test_d['Gende'].values, test_d['Scanner'].values,
    #                                    test_d['Age'].values)
    # ##
    # test_dataset = image_domian_gender(test_set['Loc'].values, test_set['Gende'].values, test_set['Scanner'].values,
    #                                    test_set['Age'].values)
    test_dataloader=DataLoader(train_dataset,batch_size=1, num_workers=4)
    evalue_model_data=evaluate_data(fea_model,age_mae,age_clas,test_dataloader,True,True)
    evalue_model_data()
###########