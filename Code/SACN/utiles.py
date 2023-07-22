# -*- coding: utf-8 -*-
"""
Created on Thu May  6 17:00:17 2021

@author: admin
"""
import torch
import  os
import numpy as np
from torch.utils.data import Dataset
from nilearn.image import get_data
import torch.nn as nn
from scipy.stats import norm
import logging
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

class image_dataset(Dataset):  # Inherit from Dataset class
    def __init__(self, input_data, target, transform=None):
        self.data = input_data
        self.target = torch.from_numpy(target).float()
        self.transform = transform

    def __getitem__(self, index):
        # read the data in this function
        x = get_data(self.data[index])
        # normalzie the data
        x = (x - np.mean(x)) / np.std(x)
        # scale the data ##
        x = (x - np.min(x)) / (np.max(x) - np.min(x))
        x = torch.from_numpy(x).float()
        y = self.target[index]
        if self.transform:
            x = self.transform(x)
        return x, y

    def __len__(self):
        return len(self.data)


class image_dataset_gender(Dataset):
    def __init__(self, input_data, gender, target, transform=None):
        self.data = input_data
        self.gender_label = torch.from_numpy(gender).float()
        self.target = torch.from_numpy(target).float()

    def __getitem__(self, index):
        x = get_data(self.data[index])
        # normalzie the data
        x = (x - np.mean(x)) / np.std(x)
        x = torch.from_numpy(x).float()
        gender = self.gender_label[index]
        y = self.target[index]
        if self.transform:
            x = self.transform(x)
        return x, [gender, y]

    def __len__(self):
        return len(self.data)


class image_domian_gender(Dataset):

    def __init__(self, input_data, gender, domain, target, transform=None):
        self.data = input_data
        self.gender_label = torch.from_numpy(gender).long()
        self.domain = torch.from_numpy(domain).long()
        self.target = torch.from_numpy(target).long()
        self.transform = transform

    def __getitem__(self, index):
        x = get_data(self.data[index])
        # normalzie the data
       # x = (x - np.mean(x)) / np.std(x)
        # scale the data ##
        x = (x - np.min(x)) / (np.max(x) - np.min(x))
        x = x.reshape((1,)+x.shape)
        x = torch.from_numpy(x).float()
        gender = self.gender_label[index]
        domain = self.domain[index]
        y = self.target[index]
        if self.transform:
            x = self.transform(x)
        return x, [gender, domain, y]

    def __len__(self):
        return len(self.data)


def my_KLDivLoss(x, y):
    """Returns K-L Divergence loss
    Different from the default PyTorch nn.KLDivLoss in that
    a) the result is averaged by the 0th dimension (Batch size)
    b) the y distribution is added with a small value (1e-16) to prevent log(0) problem
    """
    loss_func = nn.KLDivLoss(reduction='sum')
    y += 1e-16
    n = y.shape[0]
    loss = loss_func(x, y) / n
    # print(loss)
    return loss


def num2vect(x, bin_range, bin_step, sigma):
    """
    v,bin_centers = number2vector(x,bin_range,bin_step,sigma)
    bin_range: (start, end), size-2 tuple
    bin_step: should be a divisor of |end-start|
    sigma:
    = 0 for 'hard label', v is index
    > 0 for 'soft label', v is vector
    < 0 for error messages.
    """
    bin_start = bin_range[0]
    bin_end = bin_range[1]
    bin_length = bin_end - bin_start
    if not bin_length % bin_step == 0:
        print("bin's range should be divisible by bin_step!")
        return -1
    bin_number = int(bin_length / bin_step)
    bin_centers = bin_start + float(bin_step) / 2 + bin_step * np.arange(bin_number)

    if sigma == 0:
        x = np.array(x)
        i = np.floor((x - bin_start) / bin_step)
        i = i.astype(int)
        return i, bin_centers
    elif sigma > 0:
        if np.isscalar(x):
            v = np.zeros((bin_number,))
            for i in range(bin_number):
                x1 = bin_centers[i] - float(bin_step) / 2
                x2 = bin_centers[i] + float(bin_step) / 2
                cdfs = norm.cdf([x1, x2], loc=x, scale=sigma)
                v[i] = cdfs[1] - cdfs[0]
            return v, bin_centers
        else:
            v = np.zeros((len(x), bin_number))
            for j in range(len(x)):
                for i in range(bin_number):
                    x1 = bin_centers[i] - float(bin_step) / 2
                    x2 = bin_centers[i] + float(bin_step) / 2
                    cdfs = norm.cdf([x1, x2], loc=x[j], scale=sigma)
                    v[j, i] = cdfs[1] - cdfs[0]
            return v, bin_centers

def optimizer_scheduler(optimizer, p):
    """
    Adjust the learning rate of optimizer
    :param optimizer: optimizer for updating parameters
    :param p: a variable for adjusting learning rate
    :return: optimizer
    """
    for param_group in optimizer.param_groups:
        param_group['lr'] = 0.01 / (1. + 10 * p) ** 0.75

    return optimizer
## log info define #

# -*-coding:utf-8-*-
# 设置格式
def get_log(file_name):
    logger = logging.getLogger('train')
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    fh = logging.FileHandler(file_name, mode='a')
    fh.setLevel(logging.INFO)


    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger
## plot the emneding result #
class plot_embeding_res:
    def __init__(self,data,label,title):
        self.data=data
        self.label=label
        self.title=title
    @staticmethod
    def plot_embedding(label,resdata):

        # “data为n*2矩阵，label为n*1向量，对应着data的标签,title未使用”
        fig = plt.figure()
        ax = plt.subplot(111)
        type1_x = []
        type1_y = []
        type2_x = []
        type2_y = []
        type3_x = []
        type3_y = []
        type4_x = []
        type4_y = []
        type5_x = []
        type5_y = []
        for i in range(resdata.shape[0]):
            if label[i] == 0:
                type1_x.append(resdata[i][0])
                type1_y.append(resdata[i][1])
            if label[i] == 1:
                type2_x.append(resdata[i][0])
                type2_y.append(resdata[i][1])
            if label[i] == 2:
                type3_x.append(resdata[i][0])
                type3_y.append(resdata[i][1])
            if label[i] == 3:
                type4_x.append(resdata[i][0])
                type4_y.append(resdata[i][1])
            if label[i] == 4:
                type5_x.append(resdata[i][0])
                type5_y.append(resdata[i][1])
        color = plt.cm.Set3(0)
        color = np.array(color).reshape(1, 4)
        color1 = plt.cm.Set3(1)
        color1 = np.array(color1).reshape(1, 4)
        color2 = plt.cm.Set3(2)
        color2 = np.array(color2).reshape(1, 4)
        color3 = plt.cm.Set3(3)
        color3 = np.array(color3).reshape(1, 4)

        type1 = plt.scatter(type1_x, type1_y, s=10, c='r')
        type2 = plt.scatter(type2_x, type2_y, s=10, c='g')
        type3 = plt.scatter(type3_x, type3_y, s=10, c='b')
        type4 = plt.scatter(type4_x, type4_y, s=10, c='k')
        type5 = plt.scatter(type5_x, type5_y, s=10, c='c')
        # type8 = plt.scatter(type8_x, type8_y, s=10, c=color)
        # type9 = plt.scatter(type9_x, type9_y, s=10, c=color1)
        # type10 = plt.scatter(type10_x, type10_y, s=10, c=color2)
        # type11 = plt.scatter(type11_x, type11_y, s=10, c='r')
        plt.legend((type1, type2, type3, type4, type5),
                   ('IXI', 'AIBL', 'south', 'BIC', 'OASIS',),
                   loc=(0.97, 0.5))
        #plt.xticks(np.linspace(int(x_min[0]), math.ceil(x_max[0]), 5))
        #plt.yticks(np.linspace(int(x_min[1]), math.ceil(x_max[1]), 5))
        plt.xticks()
        plt.yticks()
       # plt.title(self.title)
        ax.spines['right'].set_visible(False)  # 去除右边框
        ax.spines['top'].set_visible(False)  # 去除上边框
        return fig

    def plot_2D(self):
        #
        # n_samples, n_features = data.shape
        print('Computing t-SNE embedding')
        tsne = TSNE(n_components=2, init='pca')  # TSNE process
        # t0 = time()
        result = tsne.fit_transform(self.data)  # 降维后的数据
        # print(result.shape)
        # 画图
        fig =self.plot_embedding(self.label,result)
        fig.subplots_adjust(right=0.8)  # 图例过大，保存figure时无法保存完全，故对此参数
        fig.savefig(self.title, dpi=300)

