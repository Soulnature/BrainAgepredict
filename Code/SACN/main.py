# -*- coding: utf-8 -*-
"""
Created on Wed May 12 12:31:27 2022

@author: Xingzhong
"""
import torch
import copy
from torch.autograd import Variable
import torch.nn as nn
from utiles import my_KLDivLoss,plot_embeding_res
from model import Encoder, Domain_classifier,Age_regression,Age_mae,Gende_net
import torch.optim as optim
from torch.utils.data import DataLoader
import pandas as pd
from utiles import image_domian_gender
from sklearn.model_selection import train_test_split
from config import params
from train_c import train, test_model
import os
import numpy as np
from torch.utils.tensorboard import SummaryWriter
writer = SummaryWriter(params.loss_path)
def spilt_dataset(csv_object):
    train_cat=csv_object[0:1]
    val_cat=csv_object[0:1]
    scanner_index=np.unique(csv_object['Scanner'].values)
    print(scanner_index)
    #scanner_index=[0,1,2,3,5]
    for i in scanner_index:
        oneset=csv_object[(csv_object['Scanner']==i)]
        print(oneset)
        train_data,test_data=train_test_split(oneset,test_size=0.2,random_state=0)
        train_cat=pd.concat([train_cat,train_data],axis=0)
        val_cat=pd.concat([val_cat,test_data])
    return train_cat,val_cat

def main(root_dir,network_type,con_st):
    min_loss = 10000
    best_model_fea =None
    best_model_just=None
    #torch.cuda.set_device(1)
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"
    # init the  model #
    feature_extractor = Encoder()
    Gende_pre=Gende_net()
    Age_pre = Age_regression(output_dim=90)
    Age_reg=Age_mae()
    domain_pre = Domain_classifier(params.data_num)
    # set the model to the CUDA
    feature_extractor.cuda()
    Age_pre.cuda()
    Gende_pre.cuda()
    domain_pre.cuda()
    Age_reg.cuda()
    feature_extractor = nn.DataParallel(feature_extractor)
    Age_pre = nn.DataParallel(Age_pre)
    domain_pre = nn.DataParallel(domain_pre)
    Age_reg=nn.DataParallel(Age_reg)
    Gende_pre=nn.DataParallel(Gende_pre)
    domain_criterion = nn.CrossEntropyLoss()
    if network_type=='dann':
        optimizer = optim.Adam([{'params': feature_extractor.parameters()},
                                {'params': Age_pre.parameters()},
                                {'params': domain_pre.parameters()},
                                {'params': Age_reg.parameters()},
                                {'params': Gende_pre.parameters()}], lr=params.lr)
    elif network_type=='sfcn':
         optimizer = optim.Adam([{'params': feature_extractor.parameters()},
                                {'params': Age_pre.parameters()}], lr=params.lr)
    elif network_type=='Resnet':
         optimizer = optim.Adam([{'params': feature_extractor.parameters()},
                                {'params': Age_pre.parameters()}], lr=params.lr)
    # loader the data . define the csv dir path #
    data = pd.read_csv(root_dir)
    #data=data[:1904]
  #  print(data['Loc'].values)
    train_d,test_set=spilt_dataset(data)
    ## split the dataset ##
    test_d, test_e = train_test_split(test_set, test_size=0.4,random_state=0)
    ###############
    train_dataset = image_domian_gender(train_d['Loc'].values, train_d['Gender'].values, train_d['Scanner'].values,
                                        train_d['Age'].values)
    val_dataset = image_domian_gender(test_d['Loc'].values, test_d['Gender'].values, test_d['Scanner'].values,
                                       test_d['Age'].values)
    ##
    test_dataset = image_domian_gender(test_e['Loc'].values, test_e['Gender'].values, test_e['Scanner'].values,
                                       test_e['Age'].values)
    # package the data set #
    train_dataloader = DataLoader(train_dataset, batch_size=12, num_workers=4)
    val_dataloader = DataLoader(val_dataset, batch_size=4, num_workers=4)
    ##
    test_dataloader=DataLoader(test_dataset,batch_size=4, num_workers=4)
    # start the train process
    scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.1)
    for epoch in range(params.epochs):
        print('Epoch: {}'.format(epoch))
        scheduler.step()
        train_mae,train_loss=train(network_type, feature_extractor, Age_pre, domain_pre,Age_reg,Gende_pre,domain_criterion, train_dataloader,
              optimizer, epoch, True, params.data_num, params.binnum,con_st)
        loss,cls_loss,encoder_data,label_res = test_model(network_type,feature_extractor, Age_pre, domain_pre,Age_reg,params.binnum, val_dataloader, use_gpu=True)
        print('MAE_loss:  ' + str(train_mae) + '\t train_loss:  ' + str(train_loss))
       ## write the val loss to the tensorboard #
        writer.add_scalar('./train_mae_end'+str(con_st), train_mae, epoch)
        writer.add_scalar('./test_mae_end'+str(con_st),loss,epoch)
        if loss < min_loss:
            min_loss = loss
            best_model_fea = copy.deepcopy(feature_extractor)
            best_model_just=copy.deepcopy(Age_pre)
            best_model_regre=copy.deepcopy(Age_reg)
   # writer.close()
    print("save model")
    print(min_loss)
    if network_type=='dann':
        torch.save(best_model_fea.state_dict(), './Models/best_model_fea__val_SFCN_o_5.22'+str(con_st)+'_UKB.pth')
        torch.save(best_model_just.state_dict(), './Models/best_model_sfcn_all_val_SFCN_o_end_agecls_com_5.22'+str(con_st)+'_UKB_.pth')
        torch.save(best_model_regre.state_dict(), './Models/best_model_sfcn_reg_val_SFCN_o_ukb_cons_end_mae_5.22'+str(con_st)+'_UKB_.pth')
    elif network_type=='sfcn':
        torch.save(best_model_fea.state_dict(), './Models/best_model_fea_all_test_val_sfcn.pth')
        torch.save(best_model_just.state_dict(), './Models/best_model_sfcn_all_val_sfcn.pth')
    elif network_type=='Resnet':
        torch.save(best_model_fea.state_dict(), './Models/best_model_fea_all_test_val_resNet.pth')
        torch.save(best_model_regre.state_dict(), './Models/best_model_sfcn_all_val_resNet.pth')
    ##print the test the result ##
    loss,cls_loss,encoder_data,label_res = test_model(network_type,best_model_fea, best_model_just, domain_pre,best_model_regre,params.binnum, test_dataloader, use_gpu=True)
    # plt.plot(h.history['mean_absolute_error'])
    # plt.plot(h.history['val_mean_absolute_error'])
    # plt.title('ResNet MAE')
    # plt.ylabel('MAE')
    print(loss)
#
if __name__ == '__main__':
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"
    #main(root_dir="/data/home/zhaoxz/Documents/domain_adaption/full_age1.csv", network_type='dann',con_st=con_st)
    con_st=10
    main(root_dir="/data/home/zhaoxz/UKB_match_1.csv", network_type='sacn', con_st=con_st)
    #main(root_dir="/data/home/zhaoxz/Documents/domain_adaption/full_age_pnc_com.csv",network_type='dann',con_st=con_st)
#main(root_dir="/data/home/zhaoxz/domain_adaption/ukb_train_samle.csv", network_type='dann')