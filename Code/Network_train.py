import pandas as pd
import os
import scipy.ndimage as nd
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score,mean_absolute_error
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
from scipy.stats import pearsonr
from tensorflow.keras.models import load_model
import random
import glob
import nibabel as nib
from DataLoader import dataGenerator
from collections import defaultdict
from tensorflow.keras.optimizers import Adam, SGD,Adagrad
from tensorflow.compat.v1 import reset_default_graph
from sklearn.model_selection import KFold
from Util import plotData,getPredictions,loadMR,loadHeader,calculateMeanImg
from SFCN import generateSFCN
from tensorflow.keras.callbacks import ModelCheckpoint,EarlyStopping
from tensorflow.python.keras import backend as K
import tensorflow.compat.v1 as tf
import matplotlib.pyplot as plt
tf.disable_v2_behavior()
tf.enable_eager_execution()
class train(object):
    def __init__(self, nEpochs,path,dataShape=(121,145,121),batchSize=300,decayRate=1e-6,meanImg='RawImg',gender='IncludeGender',scanner='IncludeScanner',dropRate=0.00005,Netname='SFCN',test_size=0.2,lr=0.001,maxAngle=40):
        self.nEpochs=nEpochs
        self.path=path
        self.dataShape=dataShape
        self.batchSize=batchSize
        self.decayRate=decayRate
        self.meanImg=meanImg
        self.gender=gender
        self.scanner=scanner
        self.dropRate=dropRate
        self.Netname=Netname
        self.test_size=test_size
        self.lr=lr
        self.maxAngle=maxAngle
    def __call__(self):
        data_table=pd.read_csv(self.path)
        train,val = train_test_split(data_table, test_size = self.test_size,random_state=1)
        steps_per_epoch= train.shape[0]//self.batchSize
        validation_steps = val.shape[0]//self.batchSize
        default_parameters = [0.001,1e-6,'RawImg','IncludeGender','IncludeScanner',0.00005,0.3,40,10]
        #lr, decayRate, meanImg, gender, scanner,regAmount, dropRate, maxAngle,maxShift = default_parameters
        if self.gender == 'RandomInput':
            gender_train = np.random.rand(self.train.Gender.shape[0])
            gender_val = np.random.rand(self.val.Gender.shape[0])
        else:
            gender_train = self.train.Gender.values
            gender_val = self.val.Gender.values
        if scanner == 'RandomInput':
            scanner_train = np.random.rand(self.train.Scanner.shape[0])
            scanner_val = np.random.rand(self.val.Scanner.shape[0])
        else:
            scanner_train = self.train.Scanner.values
            scanner_val = self.val.Scanner.values
        if meanImg == 'SubtractMean':
            meanImg = icelandicMeanImg
        else:
            meanImg = None
        if self.Netname=='SFCN':
            model = generateSFCN(self.dataShape,dropRate=self.dropRate)
        elif self.Netname=='ResNet':
            model = generateSFCN(self.dataShape,dropRate=self.dropRate)
        elif self.Netname=='SEResNet':
            model = generateSFCN(self.dataShape,dropRate=self.dropRate)
        #model.compile(loss='kullback_leibler_divergence',optimizer=adam, metrics=['mae','mse'])
        adam = Adam(lr=self.lr, decay=self.decayRate)
        model.compile(loss='mean_absolute_error',optimizer=adam, metrics=['mae','mse'])
        mc = ModelCheckpoint('../Models/'+self.Netname,verbose=1,mode='min',save_best_only=True)
        early = EarlyStopping(patience=200, verbose=1)
        h = model.fit(dataGenerator([train.Loc.values,scanner_train,gender_train],train.Age.values, batch_size = self.batchSize, meanImg=meanImg,dim=self.dataShape,shuffle=True,augment=True,maxAngle=self.maxAngle,maxShift=10),
                        validation_data=dataGenerator([val.Loc.values,scanner_val,gender_val],val.Age.values, batch_size = self.batchSize, meanImg=meanImg,dim=dataShape,shuffle=False,augment=False),
                        validation_steps=validation_steps,
                        steps_per_epoch=self.steps_per_epoch,
                        epochs=self.nEpochs,
                        verbose=1,
                        max_queue_size=32,
                        workers=4,
                        use_multiprocessing=False,
                        callbacks=[mc,early]
                           )
        plt.plot(h.history['loss'])
        plt.plot(h.history['val_loss'])
        plt.title('ResNet Loss')
        plt.ylabel('Loss')
        plt.xlabel('Epoch')
        plt.legend(['Train', 'Valid'], loc='upper left')
        plt.savefig("loss1_raSFCn_scanner.png")
        plt.clf()
        plt.cla()
        plt.plot(h.history['mean_absolute_error'])
        plt.plot(h.history['val_mean_absolute_error'])
        plt.title('ResNet MAE')
        plt.ylabel('MAE')
        plt.xlabel('Epoch')
        plt.legend(['Train', 'Valid'], loc='upper left')
        plt.savefig("MAE_raw"+self.Netname+".png") #save the image
        meanImg = None
        model = load_model('../Models/'+self.Netname)
        val_prediction = model.predict(dataGenerator([val.Loc.values,scanner_val,gender_val],val.Age.values, batch_size = 1, meanImg=meanImg,dim=self.dataShape,shuffle=False,augment=False),
                                verbose=1,
                                max_queue_size=32,
                                workers=4,
                                use_multiprocessing=False,
                                )
        predictions = val_prediction[:,0]
        yVal = val.Age.values
        # plt.hist(predictions)
        # plt.savefig("val_IXI_PrehisySFCN.png")
        plt.clf()
        plt.cla()
        #plt.hist(yVal)
        #plt.savefig("val_IXI_SFCn.png")
        print('Validation R^2: ',r2_score(yVal,predictions))
        print('Validation MAE: ',mean_absolute_error(yVal,predictions))
        print(pearsonr(yVal,predictions))
        plt.clf()
        plt.cla()
        y_range = np.arange(np.min(yVal),np.max(yVal))
        plt.scatter(yVal,predictions,label='T1 Prediction')
        plt.plot(y_range,y_range,c='black',ls='dashed',label='45 deg line')
        plt.xlabel('Age')
        plt.ylabel('Predicted Age')
        plt.legend()
        plt.savefig("val_SFCN_ulb_test_raw_8-90.png")

