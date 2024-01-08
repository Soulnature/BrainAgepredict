# -*- coding: utf-8 -*-
"""
Created on Thu May  6 18:26:06 2021

@author: admin
"""
import torch
from torch.autograd import Variable
import numpy as np
import torch.nn as nn
from utiles import my_KLDivLoss, num2vect, optimizer_scheduler
from config import params
from tqdm import tqdm
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
import torch.multiprocessing as mp

# Set the multiprocessing start method to 'spawn' to avoid sharing issues
mp.set_start_method('spawn', force=True)


def calculate_loss(training_mode, base_net, class_loss, mae_los, original_feature, gender, domain, constant, Gende_model, domain_model, clc_corss, DOMAIN_L):
    """
    Calculate the loss based on the training mode.

    Args:
    - training_mode (str): The current training mode.
    - base_net (str): The base network used.
    - class_loss (Tensor): The class loss.
    - mae_los (Tensor): The MAE loss.
    - original_feature (Tensor): The feature extracted by the feature extractor.
    - gender (Tensor): The gender labels.
    - domain (Tensor): The domain labels.
    - constant (float): A constant value used for computations in some modes.
    - Gende_model (nn.Module): The gender prediction model.
    - domain_model (nn.Module): The domain prediction model.
    - clc_corss (function): The cross-entropy loss function.
    - DOMAIN_L (float): The domain loss weight.

    Returns:
    - loss (Tensor): The computed loss.
    """
    if training_mode == 'base':
        if base_net == 'Resnet':
            loss = mae_los
        else:
            loss = mae_los + class_loss
    elif training_mode == 'Gender_AD':
        Gender_pre = Gende_model(original_feature, constant)
        Gender_loss = clc_corss(Gender_pre, gender)
        loss = class_loss + mae_los + (1 - DOMAIN_L) * Gender_loss
    elif training_mode == 'Domain_AD':
        domain_pre = domain_model(original_feature, constant)
        domain_loss = clc_corss(domain_pre, domain)
        loss = class_loss + mae_los + DOMAIN_L * domain_loss
    elif training_mode == 'SACN':
        Gender_pre = Gende_model(original_feature, constant)
        Gender_loss = clc_corss(Gender_pre, gender)
        domain_pre = domain_model(original_feature, constant)
        domain_loss = clc_corss(domain_pre, domain)
        loss = class_loss + mae_los + DOMAIN_L * domain_loss + (1 - DOMAIN_L) * Gender_loss
    else:
        raise ValueError(f"Unknown training mode: {training_mode}")

    return loss


def train(training_mode, models, dataloader_all, optimizer, epoch, use_gpu, bin_range, con_st, base_net, DOMAIN_L, epoche_all):
    """
    Execute training with domain adaptation.

    Args:
        training_mode (str): Training mode (base, Gender_AD, Domain_AD, SACN).
        models (dict): Dictionary containing the models (feature_extractor, age_predictor, etc.).
        dataloader_all (DataLoader): DataLoader for training data.
        optimizer (Optimizer): Optimizer for training.
        epoch (int): Current epoch.
        use_gpu (bool): Flag to use GPU if available.
        bin_range (array): Range of bins for age classification.
        con_st (float): Constant start value (unused in this function).
        base_net (str): The base network used.
        DOMAIN_L (float): Domain loss weight.
        epoche_all (int): Total number of epochs for training.

    Returns:
        float: Mean absolute error over the dataset.
        float: Total loss over the dataset.
    """

    feature_extractor = models['feature_extractor']
    age_model = models['age_predictor']
    Gende_model = models['gender_predictor']
    domain_model = models['domain_predictor']

    # Set models to training or evaluation mode based on the training_mode
    model_modes = {'Gender_AD': ('train', 'eval'), 'Domain_AD': ('eval', 'train'), 'SACN': ('train', 'train'), 'base': ('eval', 'eval')}
    Gende_model_mode, domain_model_mode = model_modes.get(training_mode, ('eval', 'eval'))
    
    getattr(Gende_model, Gende_model_mode)()
    getattr(domain_model, domain_model_mode)()
    feature_extractor.train()
    age_model.train()

    bin_step = 1
    sigma = 1
    maeloss = nn.L1Loss(reduction='mean')
    clc_corss = nn.CrossEntropyLoss()

    # Initialize loss accumulators
    loss_all = 0.0
    mae = 0.0

    # Calculate total number of steps
    start_steps = epoch * len(dataloader_all.dataset)
    total_steps = epoche_all * len(dataloader_all.dataset)

    # Initialize tqdm progress bar
    pbar = tqdm(enumerate(dataloader_all), total=len(dataloader_all), desc='Loss: 0.0000')

    for batch_idx, data in pbar:
        # Setup hyperparameters
        p = float(batch_idx + start_steps) / total_steps
        constant = 2. / (1. + np.exp(-1 * p))

        # Prepare the data
        input_data, label_all = data
        gender = label_all[0]
        domain = label_all[1]
        age_or = label_all[2].float()

        # Encode the age label and domain label
        age, _ = num2vect(age_or, bin_range, bin_step, sigma)
        age = torch.tensor(age, dtype=torch.float32).to('cuda' if use_gpu else 'cpu')

        # Reset optimizer
        optimizer = optimizer_scheduler(optimizer, p)
        optimizer.zero_grad()

        # Compute the output of source domain and target domain
        original_feature = feature_extractor(input_data.to('cuda' if use_gpu else 'cpu'))
        age_preds, age_one = age_model(original_feature)
        age_preds = age_preds.view(age_one.size()[0], -1)

        # Compute the class loss of age and MAE loss
        class_loss = my_KLDivLoss(age_preds, age)
        mae_los = maeloss(age_one[:, 0], age_or.to('cuda' if use_gpu else 'cpu'))

        # Calculate loss based on training mode
        loss = calculate_loss(training_mode, base_net, class_loss, mae_los, original_feature, gender.to('cuda' if use_gpu else 'cpu'), domain.to('cuda' if use_gpu else 'cpu'), constant, Gende_model, domain_model, clc_corss, DOMAIN_L)

        # Update overall loss and MAE
        loss_all += loss.item()
        mae += mae_los.item()

        # Backward and optimize
        loss.backward()
        optimizer.step()

        # Update tqdm progress bar description with the current loss
        pbar.set_description(f'Loss: {loss.item():.4f} | MAE: {mae_los.item():.4f}')

    # Print final average losses
    final_loss = loss_all / len(dataloader_all)
    final_mae = mae / len(dataloader_all)
    print(f'Final Loss: {final_loss:.6f}')
    print(f'Final MAE: {final_mae:.6f}')

    # Return mean absolute error and loss
    return final_mae, final_loss


## test process ##
def test_model(models, bin_range, test_dataloader, use_gpu):
    """
    Test the performance of the model
    :param models: dictionary of models including feature_extractor, age_predictor, and domain_predictor
    :param bin_range: function to divide the age into classes
    :param test_dataloader: dataloader for the test data
    :param use_gpu: boolean indicating whether to use GPU
    :return: mean absolute error (MAE) and domain accuracy
    """
    MAE = 0.0
    correct_domain = 0.0
    bin_step = 1
    sigma = 1
    save_encoder = []
    class_label = []
    maeloss = nn.L1Loss(reduction='mean')
    # Set models to evaluation mode
    for key in models:
        models[key].eval()
    
    pred_all_age = []
    original_age = []
    
    with torch.no_grad():
        for batch_idx, data in tqdm(enumerate(test_dataloader), total=len(test_dataloader)):
            # Setup hyperparameters
            p = float(batch_idx) / len(test_dataloader)
            constant = 2. / (1. + np.exp(-1 * p)) - 1.
            
            # Prepare the data
            input_data, label_all = data
            gender = label_all[0]
            domain = label_all[1]
            class_label.append(domain)
            age_or = label_all[2]
            
            # Encode the age label and domain label
            age, bin1 = num2vect(age_or, bin_range, bin_step, sigma)
            age = torch.tensor(age, dtype=torch.float32)
            gender = torch.tensor(gender, dtype=torch.int64)
            domain = torch.tensor(domain, dtype=torch.int64)
            
            # Encode domain
            if use_gpu:
                input_data, age, gender, domain,age_or = Variable(input_data.cuda()), Variable(age.cuda()), Variable(
                    gender.cuda()), Variable(domain.cuda()),Variable(age_or.cuda())
            else:
                input_data, age, gender, domain,age_or = Variable(input_data), Variable(age), Variable(gender), Variable(
                    domain)
            
            # Extract features
            encoder_data = models['feature_extractor'](input_data)
            encoder_dim = encoder_data.detach().cpu().numpy().reshape([-1, 1])
            save_encoder.append(encoder_dim[:, 0])
            # Predict age and domain
            output1, output2 = models['age_predictor'](encoder_data)
            output_data2 = output2[0].detach().cpu().numpy().reshape([1, -1])
           # age_preds = output1.reshape(output2.size()[0],90)
            class_loss = my_KLDivLoss(output1, age) 
            ##
            mae_los = maeloss(output2[:, 0], age_or)
            domian_pre = models['domain_predictor'](encoder_data, constant)
            domian_pre = domian_pre.data.max(1, keepdim=True)[1]
            correct_domain += domian_pre.eq(domain.data.view_as(domian_pre)).cpu().sum()
            mae_data_re = np.abs(output_data2 - age_or.cpu().numpy())
            out_data = output1[0].detach().cpu().numpy().reshape([1, -1])
            prob = np.exp(out_data)
            pred = prob @ bin1
            original_age.append(age_or.cpu().numpy())
            #class_loss = my_KLDivLoss(age_preds, age) 
            mae_data = np.abs(pred - age_or.cpu().numpy())
            pred_age_nu = np.asarray(((output_data2 + pred) / 2))
            #pred_age_nu = np.asarray(output_data2)
            MAE += np.mean(mae_data,mae_data_re[0])
            pred_all_age.append(pred_age_nu)

    pred_age_all = np.asarray(pred_all_age)[:, 0]
    pred_age_all = pred_age_all.astype(np.float)
    original_age = np.asarray(original_age)[:, 0]
    original_age = original_age.astype(np.float)
    
    original_age = original_age.astype(float)
    pred_age_all = pred_age_all.astype(float)   
    # Compute correlation and accuracy
   # print(pred_age_all)
    r_value = pearsonr(original_age, pred_age_all)
    r_square = r2_score(original_age, pred_age_all)
    print(r_value)
    print(r_square)
    print('\n Brain Age  Accuracy: {} ({:.4f})'
          'Domain Accuracy: {}/{} ({:.4f}%)\n'.
          format(
              MAE, float(MAE) / len(test_dataloader.dataset),
              correct_domain, len(test_dataloader.dataset),
              100. * float(correct_domain) / len(test_dataloader.dataset)
          ))
    
    return float(MAE) / len(test_dataloader.dataset), float(correct_domain) / len(test_dataloader.dataset) 
