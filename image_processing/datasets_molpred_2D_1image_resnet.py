# import necessary packages for data (image) loading and processing
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import MinMaxScaler
from tensorflow.keras import backend as K

import pandas as pd
import numpy as np

import glob
import cv2
import os

from PIL import Image



def load_house_attributes(inputPath): # we define this load_house_attribute function which accepts the path to the input dataset



    df = pd.read_csv(inputPath, delimiter = "\t", header = None)
    #df = pd.read_csv(inputPath, delimiter = " ", header = None)


    print(df)
    return df




def image_data_extract(imagePaths):


    images = []


    for (i, imagepath) in enumerate(imagePaths):

        image = cv2.imread(imagepath)
        image = cv2.resize(image, (224,224),interpolation=cv2.INTER_AREA)  # width followed by height

        image_to_append = image
        
        # cv2.imshow('img', image_to_append)
        # cv2.waitKey(0) 
        # cv2.destroyAllWindows() 

        # add the tiled images to the images list on which the network will be trained
        images.append(image_to_append)
        

    # return the set of images as an array
    return np.array(images)
    #return images/255.0
        


