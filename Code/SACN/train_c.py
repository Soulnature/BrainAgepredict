# -*- coding: utf-8 -*-
"""
Created on Thu May  6 18:26:06 2021

@author: admin
"""
import torch
from torch.autograd import Variable
import numpy as np
import torch.nn as nn
from utiles import my_KLDivLoss
from model import Encoder, Domain_classifier, Age_regression
from utiles import num2vect
from config import params
import utiles
from tqdm import tqdm
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from torch.utils.tensorboard import SummaryWriter
## init the tensorboard path
import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system')
def train(training_mode, feature_extractor, age_classifier, domain_classifier, age_regression, Gende_net,domain_criterion,
          dataloader_all, optimizer, epoch, use_gpu, domain_num, bin_range,con_st):
    """
    Execute target domain adaptation
    :param training_mode:
    :param feature_extractor:
    :param class_classifier:
    :param domain_classifier:
    :param class_criterion:
    :param domain_criterion:
    :param source_dataloader:
    :param target_dataloader:
    :param optimizer:
    :return:
    """
    # setup models
    feature_extractor.train()
    age_classifier.train()
    domain_classifier.train()
    age_regression.train()
    Gende_net.train()
    bin_step = 1
    sigma = 1
    maeloss = nn.L1Loss(reduction='mean')
    #  torch.cuda.set_device(1)
    # steps
    start_steps = epoch * len(dataloader_all.dataset)
    total_steps = params.epochs * len(dataloader_all.dataset)
    loss_all=0
    mae=0
    for batch_idx, data in enumerate(dataloader_all):
            # setup hyperparameters
          #  print('sadasdsad')
            p = float(batch_idx + start_steps) / total_steps
            constant = 2. / (1. + np.exp(-con_st* p)) - 1
            # prepare the data #
            input_data, label_all = data
            gender = label_all[0]
            domain = label_all[1]
            age_or = label_all[2].float()
            # print(domain)
            ##encoder the age label and domain label
            age, bin1 = num2vect(age_or, bin_range, bin_step, sigma)
            age = torch.tensor(age, dtype=torch.float32)
            ## domian encoder #
            if use_gpu:
                input_data, age, age_or, gender, domain = Variable(input_data.cuda()), Variable(age.cuda()), Variable(
                    age_or.cuda()), Variable(gender.cuda()), Variable(domain.cuda())
            if training_mode == 'sacn':
            # setup optimizer
                optimizer = utiles.optimizer_scheduler(optimizer, p)
                optimizer.zero_grad()
                # compute the output of source domain and target domain
                original_feature = feature_extractor(input_data)
                #print(original_feature.shape)
                # compute the class loss of age
                # real age
                age_one = age_regression(original_feature)
                age_preds = age_classifier(original_feature)
                # mae lass
                mae_los = maeloss(age_one[:, 0], age_or)
                # reshape the predicted label
                age_preds = age_preds.reshape(age_one.size()[0],90)
               ###
                out_data = age_preds[0].detach().cpu().numpy().reshape([1, -1])
                prob = np.exp(out_data)
                pred = prob @ bin1
                class_loss = my_KLDivLoss(age_preds, age)
                domian_pre = domain_classifier(original_feature, constant*constant)
                domain_loss = domain_criterion(domian_pre, domain)
                #gender loss
                Gender_pre=Gende_net(original_feature,constant*(1-constant))
                Gender_loss=domain_criterion(Gender_pre,gender)
                loss = class_loss +mae_los+constant*(domain_loss+Gender_loss)
                loss_all += loss
                mae+=mae_los
                loss.backward()
                optimizer.step()
                if (batch_idx + 1) % 10 == 0:
                    print(
                        '[{}/{} ({:.0f}%)]\t allLoss: {:.6f}\t ageClass Loss: {:.6f}\t Domain Loss: {:.6f} \t MAE Loss: {:.6f} \t Gender Loss: {:.6f}'.format(
                            batch_idx * len(input_data), len(dataloader_all.dataset),
                            100. * batch_idx / len(dataloader_all), loss.item(), class_loss.item(),
                            domain_loss.item(), mae_los.item(),Gender_loss
                        ))
        ##save the best model #
            elif training_mode == 'resnet':
                # only using the SFCN struceture predicted age
                # setup hyperparameters
                optimizer = utiles.optimizer_scheduler(optimizer, p)
                optimizer.zero_grad()
                original_feature = feature_extractor(input_data)
                # compute the class loss of age
                age_preds = age_regression(original_feature)
                mae_los = maeloss(age_preds[:, 0], age_or)
                # reshape the predicted label
                # add the MAE loss ##
                loss = mae_los
                loss_all += loss
                # print(loss)
                loss.backward()
                optimizer.step()
                # print loss
                if (batch_idx + 1) % 10 == 0:
                    print('[{}/{} ({:.0f}%)]\n allLoss: {:.6f}'.format(
                        batch_idx * len(input_data), len(dataloader_all.dataset),
                        100. * batch_idx / len(dataloader_all), loss.item()
                    ))
    return  float(mae) / len(dataloader_all), float(loss_all)/len(dataloader_all)

## test process ##
def test_model(training_mode, feature_extractor, class_classifier, domain_classifier, age_regression, bin_range,
               test_dataloader, use_gpu):
    """
    Test the performance of the model
    :param feature_extractor: network used to extract feature from target samples
    :param class_classifier: network used to predict age
    :bin_range: the fucntion to divide the age to smome class
    :param domain_classifier: network used to predict domain
    :param domain_criterion: the loss function of domain
    :return: loss
    """
    MAE = 0.0
    correct_domain = 0.0
    loss_sum = 0.0
    bin_step = 1
    sigma = 1
    save_encoder=[]
    class_label=[]
    # setup the network
    if training_mode == 'sacn':
        feature_extractor.eval()
        class_classifier.eval()
        domain_classifier.eval()
        age_regression.eval()
    elif training_mode == 'resnet':
        feature_extractor.eval()
        age_regression.eval()
    pred_all_age=[]
    original_age=[]
    with torch.no_grad():
        for batch_idx, data in enumerate(test_dataloader):
            # setup hyperparameters
            p = float(batch_idx) / len(test_dataloader)
            constant = 2. / (1. + np.exp(-10 * p)) - 1.
            # prepare the data
            input_data, label_all = data
            gender = label_all[0]
            domain = label_all[1]
            class_label.append(domain)
            # print(domain)
            age_or = label_all[2]
            ##encoder the age label and domain label
            age, bin1 = num2vect(age_or, bin_range, bin_step, sigma)
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
            encoder_data=feature_extractor(input_data)
            encoder_dim = encoder_data.detach().cpu().numpy().reshape([-1, 1])
            save_encoder.append(encoder_dim[:,0])
            if training_mode=="sacn":
                output1 = class_classifier(encoder_data)
                output2 = age_regression(encoder_data)
                output_data2=output2[0].detach().cpu().numpy().reshape([1, -1])
                domian_pre = domain_classifier(encoder_data, constant)
                ##
                domian_pre = domian_pre.data.max(1, keepdim=True)[1]
                correct_domain += domian_pre.eq(domain.data.view_as(domian_pre)).cpu().sum()
                mae_data_re = np.abs(output_data2 - age_or.numpy())
            elif training_mode=="resnet":
                output2 = age_regression(encoder_data)
                output_data2 = output2[0].detach().cpu().numpy().reshape([1, -1])
                mae_data_re = np.abs(output_data2 - age_or.numpy())
                # Compute the MAE #
            #class_loss = my_KLDivLoss(output1, age)
           # define the MAE of the data
            out_data = output1[0].detach().cpu().numpy().reshape([1, -1])
            prob = np.exp(out_data)
            pred = prob @ bin1
            original_age.append(age_or.numpy())
            mae_data = np.abs(pred - age_or.numpy())
            if training_mode=="sacn":
                pred_age_nu=np.asarray(((output_data2+pred)/2))
                #pred_age_nu = np.asarray(output_data2)
                #MAE+=mae_data_re
                MAE += (mae_data + mae_data_re) / 2
                pred_all_age.append(pred_age_nu)
            elif training_mode == 'resnet':
                pred_age_nu = np.asarray(output_data2)
                MAE += mae_data_re
                pred_all_age.append(pred_age_nu[:, 0])

    ## return th mae and r index
   # print(pred_all_age)
    pred_age_all=np.asarray(pred_all_age)[:,0]
    pred_age_all=pred_age_all.astype(np.float)
    original_age=np.asarray(original_age)[:,0]
    original_age = original_age.astype(np.float)
    #print(original_age)
    r_value=pearsonr(original_age,pred_age_all)
    r_square=r2_score(original_age,np.asarray(pred_all_age)[:,0])
    if training_mode == 'sacn':
        print(r_value)
        print(r_square)
        print('\n Source Accuracy: {} ({:.4f})'
              'Domain Accuracy: {}/{} ({:.4f}%)\n'.
            format(
            MAE, float(MAE) / len(test_dataloader.dataset),
            correct_domain, len(test_dataloader.dataset),
                 100. * float(correct_domain) / len(test_dataloader.dataset)
        ))
    else:
        print(float(MAE) / len(test_dataloader.dataset))
        print(r_value)
        print(r_square)
    if training_mode=='resnet':
        return float(MAE) / len(test_dataloader.dataset)
    else:
        return float(MAE) / len(test_dataloader.dataset),  float(correct_domain) / len(test_dataloader.dataset),np.asarray(save_encoder),class_label
