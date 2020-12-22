import pandas as pd
import os
import scipy.ndimage as nd
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score,mean_absolute_error
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
from DataLoader import dataGenerator
from scipy.stats import pearsonr
from tensorflow.keras.models import load_model
import random
import glob
import nibabel as nib
from collections import defaultdict
from tensorflow.keras.optimizers import Adam, SGD,Adagrad
from tensorflow.compat.v1 import reset_default_graph
from sklearn.model_selection import KFold
from Util import plotData,getPredictions,loadMR,loadHeader,calculateMeanImg
from ResNet import generateAgePredictionResNet
from tensorflow.keras.callbacks import ModelCheckpoint,EarlyStopping
from tensorflow.python.keras import backend as K
import tensorflow.compat.v1 as tf
import matplotlib.pyplot as plt
tf.disable_v2_behavior()
tf.enable_eager_execution()
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
dataShape=(121,145,121)
csvdir="/home2/HWGroup/zhengrc/Xingzhong/PPMI/PPMI_test.csv"
#validation data#
ukb=pd.read_csv(csvdir)
#ukb=ukb[:2000]
imageType="RawT1"
ukb = ukb.rename(columns={'Gende':'Gender'})
default_parameters = [0.001,1e-6,'RawImg','IncludeGender','IncludeScanner',0.00005,0.3,40,10]
lr, decayRate, meanImg, gender, scanner,regAmount, dropRate, maxAngle,maxShift = default_parameters
if gender == 'RandomInput':
    #gender_train = np.random.rand(train.Gender.shape[0])
    gender_val = np.random.rand(ukb.Gender.shape[0])
else:
    #gender_train = train.Gender.values
    gender_val = ukb.Gender.values
if scanner == 'RandomInput':
    #scanner_train = np.random.rand(train.Scanner.shape[0])
    scanner_val = np.random.rand(ukb.Scanner.shape[0])
else:
    #scanner_train = train.Scanner.values
    scanner_val = ukb.Scanner.values
if meanImg == 'SubtractMean':
    meanImg = icelandicMeanImg
else:
    meanImg = None
#batchExample = dataGenerator([ukb.Loc.values,ukb.Scanner.values,ukb.Gender.values],ukb.Age.values, batch_size = 4, meanImg=None,dim=dataShape,shuffle=False,augment=False,maxAngle=40,maxShift=10)
model = load_model('../Models/BrainAgeSFCN_8_90')
val_prediction = model.predict(dataGenerator([ukb.Loc.values,scanner_val,gender_val],ukb.Age.values, batch_size = 1, meanImg=meanImg,dim=dataShape,shuffle=False,augment=False),
                        verbose=1,
                        max_queue_size=32,
                        workers=4,
                        use_multiprocessing=False,
                        )
predictions = val_prediction[:,0]
yVal = ukb.Age.values
print('Validation R^2: ',r2_score(yVal,predictions))
print('Validation MAE: ',mean_absolute_error(yVal,predictions))
print(pearsonr(yVal,predictions))
y_range = np.arange(np.min(yVal),np.max(yVal))
plt.show()
plt.scatter(yVal,predictions,label='T1 Prediction')
plt.plot(y_range,y_range,c='black',ls='dashed',label='45 deg line')
plt.xlabel('Age')
plt.ylabel('Predicted Age')
plt.legend()
plt.savefig("val_SENet_PPMI_raw.png")
