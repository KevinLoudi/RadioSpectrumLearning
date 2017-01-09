# -*- coding: utf-8 -*-
"""
Created on Wed Jan 04 11:06:14 2017
Propose: DNN freamwork with theano packages
@author: kevin
Environment: Python 3.6
"""
#here would show how to build a simple ConvNet architecture 
#with some convolutional and pooling layers.
#ConvNets works as edge detectors

#first, namespace
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import urllib
import pickle
import os
import gzip
import numpy as np
import theano
import lasagne
from lasagne import layers
from lasagne.updates import nesterov_momentum
from nolearn.lasagne import NeuralNet
from nolearn.lasagne import visualize
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix

#downloading the MNIST pickled dataset and then unpacking 
#it into the three different datasets: train, validation and test.
def load_dataset():
    url = 'http://deeplearning.net/data/mnist/mnist.pkl.gz'
    filename = 'mnist.pkl.gz'
    if not os.path.exists(filename):
        print("Downloading MNIST dataset...")
        urllib.request.urlretrieve(url, filename)
    with gzip.open(filename, 'rb') as f:
        u=pickle._Unpickler(f)
        u.encoding='latin1'
        data = u.load()
    X_train, y_train = data[0]
    X_val, y_val = data[1]
    X_test, y_test = data[2]
    X_train = X_train.reshape((-1, 1, 28, 28))
    X_val = X_val.reshape((-1, 1, 28, 28))
    X_test = X_test.reshape((-1, 1, 28, 28))
    y_train = y_train.astype(np.uint8)
    y_val = y_val.astype(np.uint8)
    y_test = y_test.astype(np.uint8)
    return X_train, y_train, X_val, y_val, X_test, y_test
    
X_train, y_train, X_val, y_val, X_test, y_test = load_dataset()
#output an image
plt.imshow(X_train[0][0], cmap=cm.binary)

# define our ConvNet architecture
net1 = NeuralNet(
    layers=[('input', layers.InputLayer),
            ('conv2d1', layers.Conv2DLayer),
            ('maxpool1', layers.MaxPool2DLayer),
            ('conv2d2', layers.Conv2DLayer),
            ('maxpool2', layers.MaxPool2DLayer),
            ('dropout1', layers.DropoutLayer),
            ('dense', layers.DenseLayer),
            ('dropout2', layers.DropoutLayer),
            ('output', layers.DenseLayer),
            ],
    # input layer
    input_shape=(None, 1, 28, 28),
    # layer conv2d1
    conv2d1_num_filters=32,
    conv2d1_filter_size=(5, 5),
    conv2d1_nonlinearity=lasagne.nonlinearities.rectify,
    conv2d1_W=lasagne.init.GlorotUniform(),  
    # layer maxpool1
    maxpool1_pool_size=(2, 2),    
    # layer conv2d2
    conv2d2_num_filters=32,
    conv2d2_filter_size=(5, 5),
    conv2d2_nonlinearity=lasagne.nonlinearities.rectify,
    # layer maxpool2
    maxpool2_pool_size=(2, 2),
    # dropout1
    dropout1_p=0.5,    
    # dense
    dense_num_units=256,
    dense_nonlinearity=lasagne.nonlinearities.rectify,    
    # dropout2
    dropout2_p=0.5,    
    # output
    output_nonlinearity=lasagne.nonlinearities.softmax,
    output_num_units=10,
    # optimization method params
    update=nesterov_momentum,
    update_learning_rate=0.01,
    update_momentum=0.9,
    max_epochs=10,
    verbose=1,
    )
# Train the network
# a fully connected layer or 
nn = net1.fit(X_train, y_train)

# Prediction and Confusion Matrix
# use the model to predict the entire testing dataset
preds = net1.predict(X_test)

# plot a confusion matrix to check the performance of 
# the neural network classification
cm = confusion_matrix(y_test, preds)
plt.matshow(cm)
plt.title('Confusion matrix')
plt.colorbar()
plt.ylabel('True label')
plt.xlabel('Predicted label')
plt.show()

# visualize the 32 filters from the first convolutional layer
visualize.plot_conv_weights(net1.layers_['conv2d1'])

dense_layer = layers.get_output(net1.layers_['dense'], deterministic=True)
output_layer = layers.get_output(net1.layers_['output'], deterministic=True)
input_var = net1.layers_['input'].input_var
f_output = theano.function([input_var], output_layer)
f_dense = theano.function([input_var], dense_layer)

pred = f_output(instance)
N = pred.shape[1]
plt.bar(range(N), pred.ravel())

#from theano.tensor.nnet import conv  
#import theano.tensor as T  
#import numpy, theano  
#  
#  
#rng = numpy.random.RandomState(23455)  
#  
## symbol variable  
#input = T.tensor4(name = 'input')  
#  
## initial weights  
#w_shape = (2,3,9,9) #2 convolutional filters, 3 channels, filter shape: 9*9  
#w_bound = numpy.sqrt(3*9*9)  
#W = theano.shared(numpy.asarray(rng.uniform(low = -1.0/w_bound, high = 1.0/w_bound,size = w_shape),  
#                                dtype = input.dtype),name = 'W')  
#  
#b_shape = (2,)  
#b = theano.shared(numpy.asarray(rng.uniform(low = -.5, high = .5, size = b_shape),  
#                                dtype = input.dtype),name = 'b')  
#                                  
#conv_out = conv.conv2d(input,W)  
#  
##T.TensorVariable.dimshuffle() can reshape or broadcast (add dimension)  
##dimshuffle(self,*pattern)  
## >>>b1 = b.dimshuffle('x',0,'x','x')  
## >>>b1.shape.eval()  
## array([1,2,1,1])  
#output = T.nnet.sigmoid(conv_out + b.dimshuffle('x',0,'x','x'))  
#f = theano.function([input],output)  
#  
## demo  
#import pylab  
#from PIL import Image  
##minibatch_img = T.tensor4(name = 'minibatch_img')  
#  
##-------------img1---------------  
#img1 = Image.open(open("level_day1.jpg",mode='rb'))  
#width1,height1 = img1.size  
#img1 = numpy.asarray(img1, dtype = 'float32')/256. # (height, width, 3)  
#  
## put image in 4D tensor of shape (1,3,height,width)  
#img1_rgb = img1.swapaxes(0,2).swapaxes(1,2).reshape(1,3,height1,width1) #(3,height,width)  
#  
#  
##-------------img2---------------  
#img2 = Image.open(open("level_day1.jpg",mode='rb'))  
#width2,height2 = img2.size  
#img2 = numpy.asarray(img2,dtype = 'float32')/256.  
#img2_rgb = img2.swapaxes(0,2).swapaxes(1,2).reshape(1,3,height2,width2) #(3,height,width)  
#  
#  
#  
##minibatch_img = T.join(0,img1_rgb,img2_rgb)  
#minibatch_img = numpy.concatenate((img1_rgb,img2_rgb),axis = 0)  
#filtered_img = f(minibatch_img)  
#  
#  
## plot original image and two convoluted results  
#pylab.subplot(2,3,1);pylab.axis('off');  
#pylab.imshow(img1)  
#  
#pylab.subplot(2,3,4);pylab.axis('off');  
#pylab.imshow(img2)  
#  
#pylab.gray()  
#pylab.subplot(2,3,2); pylab.axis("off")  
#pylab.imshow(filtered_img[0,0,:,:]) #0:minibatch_index; 0:1-st filter  
#  
#pylab.subplot(2,3,3); pylab.axis("off")  
#pylab.imshow(filtered_img[0,1,:,:]) #0:minibatch_index; 1:1-st filter  
#  
#pylab.subplot(2,3,5); pylab.axis("off")  
#pylab.imshow(filtered_img[1,0,:,:]) #0:minibatch_index; 0:1-st filter  
#  
#pylab.subplot(2,3,6); pylab.axis("off")  
#pylab.imshow(filtered_img[1,1,:,:]) #0:minibatch_index; 1:1-st filter  
#pylab.show()  
