from tensorflow.keras.models import Sequential
from tensorflow.keras import backend as K
from tensorflow.keras.callbacks import TensorBoard, EarlyStopping
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, add, concatenate,GlobalAveragePooling3D,Conv3D,LeakyReLU,ELU, MaxPooling3D,AveragePooling3D,ReLU
from tensorflow.keras.regularizers import l2 as L2
from tensorflow.keras.layers import BatchNormalization as BatchNorm
from tensorflow.keras.models import Model
from tensorflow.keras.wrappers.scikit_learn import KerasRegressor
from tensorflow.keras.optimizers import Adam, SGD,Adagrad
#import tensorflow.compat.v1 as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
#@staticmethod
def littlebuildblock(inputs,feature,kernalsize,maxpool):
    if maxpool==True:
        hidden=Conv3D(feature,kernalsize,padding='same',kernel_initializer='he_uniform')(inputs)
        hidden=BatchNorm(renorm=True)(hidden)
        hidden=MaxPooling3D(pool_size=(2,2,2),strides=(2,2,2), padding='same')(hidden)
        output=ReLU()(hidden)
    else:
        hidden=Conv3D(feature,kernalsize,padding='same',kernel_initializer='he_uniform')(inputs)
        hidden=BatchNorm(renorm=True)(hidden)
        output=ReLU()(hidden)
    return output
def generateSFCN(inputShape,dropRate=0.3,includeScannerGender=True,output_dim=1):

    channel_number=[32, 64, 128, 256, 256, 64]
    n_layer = len(channel_number)
    T1input = Input(inputShape+(1,), name='T1_Img')
    X=littlebuildblock(T1input,32,(3,3,3),True)
    for i in range(1,n_layer):
            if i == 0:
                in_channel = 1
            else:
                in_channel = channel_number[i-1]
            out_channel = channel_number[i]
            if i < n_layer-1:
                X=littlebuildblock(X,out_channel,(3, 3, 3),True)
            else:
                X=littlebuildblock(X,out_channel,(1, 1, 1),False)
    ##classify layers
    avg_shape = (2, 2, 2)
    X=AveragePooling3D(avg_shape, data_format='channels_first')(X)
    output=Dropout(dropRate)(X)
    outputchannel = output_dim
    #print(output.shape)
    output=Conv3D(outputchannel,(1,1,1))(output)
    output = Flatten(name='flatten1')(output)
    print(output.shape)
    if includeScannerGender:
        scanner  = Input((1,), name='Scanner')
        gender  = Input((1,), name='Gender')
        output = Dense(128,name='FullyConnectedLayer')(output)
        output = concatenate([scanner,gender,output])
        prediction = Dense(output_dim,name='Regression')(output)
    else:
        prediction = Dense(output_dim,name='Regression')(output)
    if includeScannerGender:
        model = Model(inputs=[T1input,scanner,gender],outputs=prediction)
    else:
        model = Model(inputs=[T1input],outputs=prediction)
    return model
