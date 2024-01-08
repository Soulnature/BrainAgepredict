# -*- coding: utf-8 -*-

# utility params
import os
fig_mode = None
embed_plot_epoch=20

# model params
use_gpu = True
dataset_mean = (0.5, 0.5, 0.5)
dataset_std = (0.5, 0.5, 0.5)
loss_path='/data/home/zhaoxz/domain_adaption/Pic/'
gamma = 10
theta = 1
binnum=[5,95]
data_num=5
BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
