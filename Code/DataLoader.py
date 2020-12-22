import threading
from transformations import rotation_matrix
from random import gauss
from scipy.ndimage.interpolation import map_coordinates
import operator
from scipy.ndimage.interpolation import shift,rotate
import pandas as pd
import os
import scipy.ndimage as nd
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score,mean_absolute_error
import numpy as np
import matplotlib.pyplot as plt
import math
import sys
from tensorflow.keras.models import load_model
import random
import glob
import nibabel as nib
from collections import defaultdict
from tensorflow.keras.utils import Sequence

def loadMR(path):
    #img = nib.load(path).get_data()
    img = nib.load(path).get_fdata()
    img=(img - np.mean(img)) / np.std(img)
   # img=img[1:120,1:144,1:120]
    return img
class dataGenerator(Sequence):
    'Generates data for Keras'
    def __init__(self, features, labels, batch_size=32,meanImg=None, dim=(121, 145, 121),maxAngle=40,maxShift=10, shuffle=True,augment=False,includeScannerGender=True):
        'Initialization'
        self.batch_size = batch_size
        self.features = features
        self.labels = labels
        self.dim = dim
        self.meanImg = meanImg
        self.augment = augment
        self.maxAngle = maxAngle
        self.maxShift = maxShift
        self.shuffle = shuffle
        self.IncludeScannerGender = includeScannerGender
        self.on_epoch_end()
        
    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(self.features[0].shape[0] / self.batch_size))
    
    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.features[0]))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)
            
    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        #print(index)
        index = index%self.__len__()
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        #features_temp = [self.features[k] for k in indexes]

        # Generate data
        X, y = self.__data_generation(indexes)

        return X, y

    def __data_generation(self, indexes):
        'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
        # Initialization
        #X = np.empty((self.batch_size, self.dim[0],self.dim[1],self.dim[2]),dtype=np.uint8)
        X = np.empty((self.batch_size, self.dim[0],self.dim[1],self.dim[2], 1))
        age = np.empty((self.batch_size))
        if self.IncludeScannerGender: 
            sex = np.empty((self.batch_size))
            scanner = np.empty((self.batch_size))
        # Generate data
        for i, index in enumerate(indexes):
            X[i,:,:,:,:] = processing(self.features[0][index],self.dim,self.meanImg,augment=self.augment)
            age[i] = self.labels[index]
            if self.IncludeScannerGender: 
                scanner[i] = self.features[1][index]
                sex[i] = self.features[2][index]
                
        
        if self.IncludeScannerGender: 
            return [X,scanner,sex], [age]
        else:
            return [X], [age]


def resize3d(image,new_shape,order=3):
    real_resize_factor = tuple(map(operator.truediv, new_shape, image.shape))
    return nd.zoom(image, real_resize_factor, order=order)

def make_rand_vector(dims):
    vec = [gauss(0, 1) for i in range(dims)]
    mag = sum(x**2 for x in vec) ** .5
    return [x/mag for x in vec]

def coordinateTransformWrapper(X_T1,maxDeg=40,maxShift=7.5):
    randomAngle = np.radians(maxDeg*2*(random.random()-0.5))
    unitVec = tuple(make_rand_vector(3))
    shiftVec = [maxShift*2*(random.random()-0.5),
                maxShift*2*(random.random()-0.5),
                maxShift*2*(random.random()-0.5)]
    X_T1 = coordinateTransform(X_T1,randomAngle,unitVec,shiftVec)
    return X_T1

def coordinateTransform(vol,randomAngle,unitVec,shiftVec,order=1,mode='constant'):
    #from transformations import rotation_matrix
    ax = (list(vol.shape))
    ax = [ ax[i] for i in [1,0,2]]
    coords=np.meshgrid(np.arange(ax[0]),np.arange(ax[1]),np.arange(ax[2]))

    # stack the meshgrid to position vectors, center them around 0 by substracting dim/2
    xyz=np.vstack([coords[0].reshape(-1)-float(ax[0])/2,     # x coordinate, centered
               coords[1].reshape(-1)-float(ax[1])/2,     # y coordinate, centered
               coords[2].reshape(-1)-float(ax[2])/2,     # z coordinate, centered
               np.ones((ax[0],ax[1],ax[2])).reshape(-1)])    # 1 for homogeneous coordinates
    
    # create transformation matrix
    mat=rotation_matrix(randomAngle,unitVec)

    # apply transformation
    transformed_xyz=np.dot(mat, xyz)

    # extract coordinates, don't use transformed_xyz[3,:] that's the homogeneous coordinate, always 1
    x=transformed_xyz[0,:]+float(ax[0])/2+shiftVec[0]
    y=transformed_xyz[1,:]+float(ax[1])/2+shiftVec[1]
    z=transformed_xyz[2,:]+float(ax[2])/2+shiftVec[2]
    x=x.reshape((ax[1],ax[0],ax[2]))
    y=y.reshape((ax[1],ax[0],ax[2]))
    z=z.reshape((ax[1],ax[0],ax[2]))
    new_xyz=[y,x,z]
    new_vol=map_coordinates(vol,new_xyz, order=order,mode=mode)
    return new_vol


def processing(features,inputShape,meanImg,maxAngle=40,maxShift=10,resizeImg=False,augment=False,training=False):
    X_T1 = loadMR(features)
    if meanImg is not None:
        X_T1 = X_T1-meanImg
    if augment:
        X_T1 = coordinateTransformWrapper(X_T1,maxDeg=maxAngle,maxShift=maxShift)
    if resizeImg:
        inputShape = (121, 145, 121)
        X_T1 = resize3d(X_T1,inputShape)
    return X_T1.reshape(inputShape+(1,))
