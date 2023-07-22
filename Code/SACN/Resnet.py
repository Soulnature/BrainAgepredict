# -*- coding: utf-8 -*-
import logging
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.model_zoo as model_zoo
logger = logging.getLogger(__name__)


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


__all__ = ['ResNet', 'resnet18', 'resnet34', 'resnet50', 'resnet101',
           'resnet152']



def conv3x3(in_planes, out_planes, stride=1):
    """3x3 convolution with padding"""
    return nn.Conv3d(in_planes, out_planes, kernel_size=3, stride=stride,
                     padding=1, bias=False)


class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(BasicBlock, self).__init__()
        self.conv1 = conv3x3(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm3d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv3x3(planes, planes)
        self.bn2 = nn.BatchNorm3d(planes)
        self.downsample = downsample
        self.stride = stride

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            residual = self.downsample(x)

        out += residual
        out = self.relu(out)

        return out


class Bottleneck(nn.Module):
    expansion = 4

    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(Bottleneck, self).__init__()
        self.conv1 = nn.Conv3d(inplanes, planes, kernel_size=1, bias=False)
        self.bn1 = nn.BatchNorm3d(planes)
        self.conv2 = nn.Conv3d(planes, planes, kernel_size=3, stride=stride,
                               padding=1, bias=False)
        self.bn2 = nn.BatchNorm3d(planes)
        self.conv3 = nn.Conv3d(planes, planes * self.expansion, kernel_size=1, bias=False)
        self.bn3 = nn.BatchNorm3d(planes * self.expansion)
        self.relu = nn.ReLU(inplace=True)
        self.downsample = downsample
        self.stride = stride

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)
        out = self.relu(out)

        out = self.conv3(out)
        out = self.bn3(out)

        if self.downsample is not None:
            residual = self.downsample(x)

        out += residual
        out = self.relu(out)

        return out


class ResNet(nn.Module):
    def __init__(self, block, layers, num_classes=1000,
                 channel_size=[64,64,128,256,512],
                 dropout=False):
        c = channel_size
        self.inplanes = c[0]
        super(ResNet, self).__init__()
        net = nn.Sequential()
        net.add_module('conv1', nn.Conv3d(1, c[0], kernel_size=7, stride=2, padding=0,
                               bias=False))
        net.add_module('bn1', nn.BatchNorm3d(c[0]))
        net.add_module('relu', nn.ReLU(inplace=True))
        net.add_module('maxpool', nn.MaxPool3d(kernel_size=3, stride=2, padding=1))
        net.add_module('layer1', self._make_layer(block, c[1], layers[0]))
        net.add_module('layer2', self._make_layer(block, c[2], layers[1], stride=2))
        net.add_module('layer3', self._make_layer(block, c[3], layers[2], stride=2))
        net.add_module('layer4', self._make_layer(block, c[4], layers[3], stride=2))
        net.add_module('avgpool', nn.AvgPool3d([3,4,3], stride=1))
        if dropout is True:
            net.add_module('dropout', nn.Dropout(0.5))
        self.feature_extractor = net
        for m in self.modules():
            if isinstance(m, nn.Conv3d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
            elif isinstance(m, nn.BatchNorm3d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)

    def _make_layer(self, block, planes, blocks, stride=1):
        downsample = None
        if stride != 1 or self.inplanes != planes * block.expansion:
            downsample = nn.Sequential(
                nn.Conv3d(self.inplanes, planes * block.expansion,
                          kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm3d(planes * block.expansion),
            )

        layers = []
        layers.append(block(self.inplanes, planes, stride, downsample))
        self.inplanes = planes * block.expansion
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes))

        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.feature_extractor(x)
        return x


def resnet18(**kwargs):
    """Constructs a ResNet-18 model.
    Args:
    """
    model = ResNet(BasicBlock, [2, 2, 2, 2], **kwargs)
    return model


def resnet34(**kwargs):
    """Constructs a ResNet-34 model.
    Args:
    """
    model = ResNet(BasicBlock, [3, 4, 6, 3], **kwargs)
    return model


def resnet50(**kwargs):
    """Constructs a ResNet-50 model.
    Args:
    """
    model = ResNet(Bottleneck, [3, 4, 6, 3], **kwargs)
    return model


def resnet101(**kwargs):
    """Constructs a ResNet-101 model.
    Args:
    """
    model = ResNet(Bottleneck, [3, 4, 23, 3], **kwargs)
    return model


def resnet152(**kwargs):
    """Constructs a ResNet-152 model.
    Args:
    """
    model = ResNet(Bottleneck, [3, 8, 36, 3], **kwargs)
    return model

class Age_regression(nn.Module):
    def __init__(self, channel_number=[32, 64, 128, 256, 256, 128], output_dim=90, dropout=True):
        super(Age_regression, self).__init__()
        n_layer = len(channel_number)
        ##define the age classification model (1) regresison model (2) classification model #
        self.age_regression = nn.Sequential()
        avg_shape = [2, 2, 2]
        self.age_regression.add_module('average_pool', nn.AvgPool3d(avg_shape))
        if dropout is True:
            self.age_regression.add_module('dropout', nn.Dropout(0.5))
        i = n_layer
        in_channel = 2048
        self.age_regression.add_module('conv_%d' % i,
                                   nn.Conv3d(in_channel,output_dim, padding=0, kernel_size=1))

    def forward(self, x):
        # out = list()
        x_out = self.age_regression(x)
        x = F.log_softmax(x_out, dim=1)
        # out.append(x)
        return x
class Age_mae(nn.Module):
    def __init__(self):
        super(Age_mae,self).__init__()
        self.fn=nn.Flatten()
        self.fc1 = nn.Linear(2048*8, 512)
        self.fc2 = nn.Linear(512, 1)
       # self.fc3 = nn.Linear(64, 1)
    def forward(self, x):
        x = self.fn(x)
        x=F.relu(self.fc1(F.dropout(x)))
        x=F.relu(self.fc2(x))
     #   out=self.fc3(x)
        return x
class Gende_net(nn.Module):
    def __init__(self):
        super(Gende_net,self).__init__()
        self.fc1 = nn.Linear(2048*8, 128)
        self.fc2 = nn.Linear(128, 64)
        self.fc3 = nn.Linear(64, 2)
    def forward(self, x,constant):
        x = x.view(-1, 2048*8)
        x = GradReverse.grad_reverse(x, constant)
        x=F.relu(self.fc1(F.dropout(x)))
        x=F.relu(self.fc2(x))
        out=F.log_softmax(self.fc3(x))
        return out
class Domain_classifier(nn.Module):

    def __init__(self, domain_num):
        super(Domain_classifier, self).__init__()  ## for the domian adaption and updata the gradient
        self.fc1 = nn.Linear(2048*8, 128)
        self.fc3 = nn.Linear(128, domain_num)
    def forward(self, x, constant):
        #flattern the cnn blocks
        x = x.view(-1,2048*8)
        x = GradReverse.grad_reverse(x, constant)
        logits = F.relu(self.fc1(x))
        logits =self.fc3(logits)
        return logits
