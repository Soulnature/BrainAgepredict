import pandas as pd
import numpy as np
import os
import scipy.ndimage as nd
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score,mean_absolute_error
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
from tensorflow.keras.models import load_model
from tensorflow.keras.models import Model
import random
import glob
import nibabel as nib
from collections import defaultdict
from tensorflow.keras.optimizers import Adam, SGD,Adagrad
from tensorflow.compat.v1 import reset_default_graph
from sklearn.model_selection import KFold
from DataLoader import dataGenerator,getIcelandicData,getIXIData,getUKBData
from Util import plotData,getPredictions,loadMR,loadHeader,calculateMeanImg
from ResNet import generateAgePredictionResNet
from tensorflow.keras.callbacks import ModelCheckpoint,EarlyStopping
from tensorflow.python.keras import backend as K
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
tf.enable_eager_execution()
from tensorflow.keras.models import load_model
dataShape = (119, 143, 119)
#csvdir="/home2/HWGroup/zhengrc/Xingzhong/BrainAgePredictionResNet-master/Code/oldUkbDatainfoGrey.csv"
csvdir="/home2/HWGroup/zhengrc/Xingzhong/BrainAgePredictionResNet-master/Code/IXIdata_info.csv"
ukb=pd.read_csv(csvdir)
imageType="GrayMatter"
#ukb=ukb[:1000]
ukb = ukb.rename(columns={'Gende':'Gender'})
# default_parameters = [0.001,1e-6,'GrayMatter','IncludeGender','IncludeScanner',0.00005,0.2,40,10]
# lr, decayRate, meanImg, gender, scanner,regAmount, dropRate, maxAngle,maxShift = default_parameters
#model= generateAgePredictionResNet(dataShape,regAmount=regAmount,dropRate=dropRate)
batchExample = dataGenerator([ukb.Loc.values,ukb.Scanner.values,ukb.Gender.values],ukb.Age.values, batch_size = 1, meanImg=None,dim=dataShape,shuffle=False,augment=False,maxAngle=40,maxShift=10)
model = load_model('../Models/BrainAgeResNet(GrayMatter-Ice-TransferLearningOnIXI)')
# model = load_model('../Models/BrainAgeResNetUKB1(GrayMatter-Ice)')
# model.summary()
#model.load_weights('../Models/BrainAgeResNetUKB1(GrayMatter-Ice)')
#model = load_model('../Models/BrainAgeResNetUKB1(GrayMatter-Ice)')
model.summary()
ukb_prediction_TL = model.predict(batchExample,
                        verbose=1,
                        max_queue_size=32,
                        workers=4,
                        use_multiprocessing=True,
                        )
np.save("IXI_predict.npy",ukb_prediction_TL)
# intermediate_layer_model = Model(inputs=model.input,outputs=model.get_layer('flatten1').output)
#intermediate_layer_model2 = Model(inputs=model.input,outputs=model.get_layer('FullyConnectedLayer').output)
# inter = intermediate_layer_model.predict(batchExample, verbose=1,max_queue_size=32, workers=4,)
# #inter2 = intermediate_layer_model2.predict(batchExample, verbose=1,max_queue_size=32, workers=4,)
# print(inter)
#np.save("lastLayeroutNewmodelFlateen.npy",inter)
#np.save("lastLayeroutNewmodelFunc.npy",inter2)
plt.clf()
plt.cla()
predictions= ukb_prediction_TL[:,0]
y = ukb.Age.values
print('R^2: ',r2_score(y,predictions))
print('MAE: ',mean_absolute_error(y,predictions))
plt.scatter(y,predictions,label='Predictions')
y_range = np.arange(20,np.max(y))
plt.plot(y_range,y_range,c='black',ls='dashed',label='45 deg line')
plt.xlabel('Age')
plt.ylabel('Predicted Age')
plt.title('Prediction with transfer learning (IXIdata)')
plt.legend()
plt.savefig("IXIdata.png")
plt.show()
