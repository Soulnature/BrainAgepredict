# -*- coding: utf-8 -*-
import logging
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Function
from typing import Optional, Any, Tuple

logger = logging.getLogger(__name__)



class ChannelSELayer3D(nn.Module):
    """
    3D extension of Squeeze-and-Excitation (SE) block described in:
        *Hu et al., Squeeze-and-Excitation Networks, arXiv:1709.01507*
        *Zhu et al., AnatomyNet, arXiv:arXiv:1808.05238*
    """

    def __init__(self, num_channels, reduction_ratio=4):
        """
        :param num_channels: No of input channels
        :param reduction_ratio: By how much should the num_channels should be reduced
        """
        super(ChannelSELayer3D, self).__init__()
        self.avg_pool = nn.AdaptiveAvgPool3d(1)
        num_channels_reduced = num_channels // reduction_ratio
        self.reduction_ratio = reduction_ratio
        self.fc1 = nn.Linear(num_channels, num_channels_reduced, bias=True)
        self.fc2 = nn.Linear(num_channels_reduced, num_channels, bias=True)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, input_tensor):
        """
        :param input_tensor: X, shape = (batch_size, num_channels, D, H, W)
        :return: output tensor
        """
        batch_size, num_channels, D, H, W = input_tensor.size()
        # Average along each channel
        squeeze_tensor = self.avg_pool(input_tensor)

        # channel excitation
        fc_out_1 = self.relu(self.fc1(squeeze_tensor.view(batch_size, num_channels)))
        fc_out_2 = self.sigmoid(self.fc2(fc_out_1))

        output_tensor = torch.mul(input_tensor, fc_out_2.view(batch_size, num_channels, 1, 1, 1))

        return output_tensor




class GradReverse(torch.autograd.Function):
    """
    Extension of grad reverse layer
    """

    @staticmethod
    def forward(ctx, x, constant):
        ctx.constant = constant
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        grad_output = grad_output.neg() * ctx.constant
        return grad_output, None

    def grad_reverse(x, constant):
        return GradReverse.apply(x, constant)




## SE BLOCK #
class Encoder(nn.Module):
    # 32, 64, 128, 256, 256, 128
    def __init__(self, channel_number=[32, 64, 128, 256,256,128], dropout=True):
        super(Encoder, self).__init__()
        n_layer = len(channel_number)
        self.feature_extractor = nn.Sequential()  # Define the feature extractor using sample convolution structure
        for i in range(n_layer):
            if i == 0:
                in_channel = 1
            else:
                in_channel = channel_number[i - 1]
            out_channel = channel_number[i]
            if i < n_layer - 1:
                self.feature_extractor.add_module('conv_%d' % i,
                                                  self.conv_layer(in_channel,
                                                                  out_channel,
                                                                  maxpool=True,
                                                                  kernel_size=3,
                                                                  padding=1))
            else:
                self.feature_extractor.add_module('conv_%d' % i,
                                                  self.conv_layer(in_channel,
                                                                  out_channel,
                                                                  maxpool=False,
                                                                  kernel_size=1,
                                                                  padding=0))
    @staticmethod
    def conv_layer(in_channel, out_channel, maxpool=True, kernel_size=3, padding=0, maxpool_stride=2):
        if maxpool is True:
            layer = nn.Sequential(
                nn.Conv3d(in_channel, out_channel, padding=padding, kernel_size=kernel_size),
                nn.BatchNorm3d(out_channel),
                nn.MaxPool3d(2, stride=maxpool_stride),
                nn.ReLU()

            )
        else:
            layer = nn.Sequential(
                nn.Conv3d(in_channel, out_channel, padding=padding, kernel_size=kernel_size),
                nn.BatchNorm3d(out_channel),
                nn.ReLU(),
            )
        return layer

    def forward(self, x):
        x_f = self.feature_extractor(x)
        return x_f

### 32, 64, 128, 256, 256, 128]
class Age_regression(nn.Module):
    def __init__(self, channel_number=[32, 64, 128, 256,256,128], output_dim=90, dropout=True):
        super(Age_regression, self).__init__()
        n_layer = len(channel_number)
        ##define the age classification model (1) regresison model (2) classification model #
        self.age_regression = nn.Sequential()
        avg_shape = [3, 4, 3]
        self.age_regression.add_module('average_pool', nn.AvgPool3d(avg_shape))
        if dropout is True:
            self.age_regression.add_module('dropout', nn.Dropout(0.5))
        i = n_layer
        in_channel = channel_number[-1]
        self.age_regression.add_module('conv_%d' % i,
                                   nn.Conv3d(in_channel,output_dim, padding=0, kernel_size=1))
        self.agemae=Age_mae()

    def forward(self, x):
        x_out = self.age_regression(x)
        y=self.agemae(x)
        x = F.log_softmax(x_out, dim=1)
        return x,y
class Age_mae(nn.Module):
    def __init__(self):
        super(Age_mae,self).__init__()
        self.fc1 = nn.Linear(128*3*4*3, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, 1)
    def forward(self, x):
        x = x.view(-1, 128*3*4*3)
        x=F.relu(self.fc1(F.dropout(x)))
        x=F.relu(self.fc2(x))
        out=F.relu(self.fc3(x))
        return out
class Gende_net(nn.Module):
    def __init__(self):
        super(Gende_net,self).__init__()
        self.fc1 = nn.Linear(128*3*4*3, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, 2)
    def forward(self, x,constant):
        x = x.view(-1, 128*3*4*3)
        x = GradReverse.grad_reverse(x, constant)
        x=F.relu(self.fc1(F.dropout(x)))
        x=F.relu(self.fc2(x))
        out=F.log_softmax(self.fc3(x))
        return out
class Domain_classifier(nn.Module):

    def __init__(self, domain_num):
        super(Domain_classifier, self).__init__()  ## for the domian adaption and updata the gradient
        self.fc1 = nn.Linear(128*3*4*3, 128)
        self.fc3 = nn.Linear(128, domain_num)
    def forward(self, x, constant):
        #flattern the cnn blocks
        x = x.view(-1, 128*3*4*3)
        x = GradReverse.grad_reverse(x, constant)
        logits = F.relu(self.fc1(x))
        logits =self.fc3(logits)
        return logits
