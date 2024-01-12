import os
import copy
import torch
import pandas as pd
import numpy as np
import argparse
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from SACN import Encoder, Domain_classifier, Age_regression, Gende_net
from Resnet import resnet50,Age_regression_re,Gende_net_re,Domain_classifier_re,resnet101
from utiles import image_domian_gender,EarlyStopping
from config import params

from train_c import train, test_model
import datetime
MODEL_PATH='./save_Models/'
def split_dataset(csv_object):
    train_cat = csv_object[0:1]
    val_cat = csv_object[0:1]
    scanner_index = np.unique(csv_object['Scanner'].values)

    for i in scanner_index:
        one_set = csv_object[(csv_object['Scanner'] == i)]
        train_data, test_data = train_test_split(one_set, test_size=0.3, random_state=0)
        train_cat = pd.concat([train_cat, train_data], axis=0)
        val_cat = pd.concat([val_cat, test_data])
    return train_cat, val_cat

def initialize_model(base,lr):
    
    if base=='SACN':
        feature_extractor = Encoder().cuda()
        gender_predictor = Gende_net().cuda()
        age_predictor = Age_regression(output_dim=90).cuda()
        domain_predictor = Domain_classifier(params.data_num).cuda()
    elif base=='Resnet':
        feature_extractor = resnet50().cuda()
        gender_predictor = Gende_net_re().cuda()
        age_predictor = Age_regression_re(output_dim=90).cuda()
        domain_predictor = Domain_classifier_re(params.data_num).cuda()

    models = {
        'feature_extractor': nn.DataParallel(feature_extractor),
        'age_predictor': nn.DataParallel(age_predictor),
        'gender_predictor': nn.DataParallel(gender_predictor),
        'domain_predictor': nn.DataParallel(domain_predictor),
    }

    optimizer = optim.Adam([
            {'params': feature_extractor.parameters()},
            {'params': age_predictor.parameters()},
            {'params': domain_predictor.parameters()},
            {'params': gender_predictor.parameters()}
        ], lr=lr)
 

    return models, optimizer

def main(root_dir, network_type, con_st,DOMAIN_L,model_name,base,epoches,batch_size,lr):
    min_loss = 10000
    best_model_fea = None
    best_model_just = None
    #get the data time #
    print(model_name)
    models, optimizer = initialize_model(base,lr)

    data = pd.read_csv(root_dir)
   # data=data[0:3000]
    train_d, test_set = split_dataset(data)
    test_d, test_e = train_test_split(test_set, test_size=0.66, random_state=0)

    train_dataset = image_domian_gender(train_d['Loc'].values, train_d['Gender'].values, train_d['Scanner'].values, train_d['Age'].values)
    val_dataset = image_domian_gender(test_d['Loc'].values, test_d['Gender'].values, test_d['Scanner'].values, test_d['Age'].values)
    test_dataset = image_domian_gender(test_e['Loc'].values, test_e['Gender'].values, test_e['Scanner'].values, test_e['Age'].values)

    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, num_workers=4)
    val_dataloader = DataLoader(val_dataset, batch_size=1, num_workers=4)
    test_dataloader = DataLoader(test_dataset, batch_size=1, num_workers=4)

    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.1)

    for epoch in range(epoches):
        print(f'Epoch: {epoch}')
        scheduler.step()
        train_mae, train_loss = train(network_type, models, train_dataloader, optimizer, epoch, True, params.binnum, con_st,base,DOMAIN_L,epoches)
        loss, cls_loss = test_model(models, params.binnum,val_dataloader, use_gpu=True)
       # early_stopping(loss, models)
        print(f'MAE_loss: {train_mae}\t train_loss: {train_loss}')

        if loss < min_loss:
            min_loss = loss
            best_model_fea = copy.deepcopy(models['feature_extractor'])
            best_model_just = copy.deepcopy(models['age_predictor'])
            best_model_domain = copy.deepcopy(models['domain_predictor'])

    print("save model")
    save_model(best_model_fea, best_model_just,model_name,network_type,base)
    model={'feature_extractor':best_model_fea,
           "age_predictor":best_model_just,
           "domain_predictor":best_model_domain}
    loss, cls_loss = test_model(model, params.binnum, test_dataloader, use_gpu=True)
    print(loss)
    return loss

def save_model(best_model_fea, best_model_just,model_name,network_type,base):
        
        torch.save(best_model_fea.state_dict(), os.path.join(MODEL_PATH,model_name+str(base)+str(network_type)+'_encoder'))
        torch.save(best_model_just.state_dict(), os.path.join(MODEL_PATH,model_name+str(base)+str(network_type)+'age_predictor'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--root_dir', type=str, required=True, help='Root directory')
    parser.add_argument('--network_type', type=str, required=True, help='Network type')
    parser.add_argument('--model_name', type=str, required=True, help='Model name')
    parser.add_argument('--base', type=str, required=True, help='Base')
    parser.add_argument('--epochs', type=int, required=True, help='Number of epochs')
    parser.add_argument('--batch_size', type=int, required=True, help='Batch size')
    parser.add_argument('--lr', type=float, required=True, help='Learning rate')
    args = parser.parse_args()

    ISOTIMEFORMAT = '%Y-%m-%d'
    theTime = datetime.datetime.now().strftime(ISOTIMEFORMAT)
    os.environ["CUDA_VISIBLE_DEVICES"] = '4'
    
    tem_res=main(root_dir=args.root_dir, network_type=args.network_type,con_st=1,DOMAIN_L=1,model_name=theTime+args.model_name,base=args.base,epoches=args.epochs,batch_size=args.batch_size,lr=args.lr)

