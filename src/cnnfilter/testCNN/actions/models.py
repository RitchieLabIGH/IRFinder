import os,pickle
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow.keras import  layers

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import json
import numpy as np

from progress.bar import Bar
from utils.reader import ImageArchive
from actions.extract import getImageArrayFromRegion

def conv2d_bn(x,
              filters,
              num_row,
              num_col,
              padding='same',
              strides=(1, 1),
              name=None):
    """Utility function to apply conv + BN.
    # Arguments
        x: input tensor.
        filters: filters in `Conv2D`.
        num_row: height of the convolution kernel.
        num_col: width of the convolution kernel.
        padding: padding mode in `Conv2D`.
        strides: strides in `Conv2D`.
        name: name of the ops; will become `name + '_conv'`
            for the convolution and `name + '_bn'` for the
            batch norm layer.
    # Returns
        Output tensor after applying `Conv2D` and `BatchNormalization`.
    """
    if name is not None:
        bn_name = name + '_bn'
        conv_name = name + '_conv'
    else:
        bn_name = None
        conv_name = None
    #===========================================================================
    # if backend.image_data_format() == 'channels_first':
    #     bn_axis = 1
    # else:
    #===========================================================================
    bn_axis = 3
    x = layers.Conv2D(
        filters, (num_row, num_col),
        strides=strides,
        padding=padding,
        use_bias=False,
        name=conv_name)(x)
    x = layers.BatchNormalization(axis=bn_axis, scale=False, name=bn_name)(x)
    x = layers.Activation('relu', name=name)(x)
    return x

class IntronModeller():
    def __init__(self, verbosity=0):
        if verbosity > 3:
            verbosity=3
        os.environ['TF_CPP_MIN_LOG_LEVEL'] ="{}".format(3-verbosity)
        self.verbosity = verbosity


    def _get_model(self, json_model):
        if json_model == None:

            inputs = layers.Input(shape=(self.img_size, 2, 1))
            #rescale replace by preprocess
            x = layers.experimental.preprocessing.Rescaling(1./255,offset=-1, input_shape=(self.img_size, 2, 1))(inputs)

            x = conv2d_bn(inputs, 64, 3, 2, strides=(1, 1), padding='valid')
            #x = layers.LayerNormalization(axis=1)(x)
            x = layers.MaxPooling2D((3, 1), strides=(2, 1))(x)
            x = conv2d_bn(x, 32, 3, 1, strides=(1, 1))
            x = layers.MaxPooling2D((3, 1), strides=(2, 1))(x)
            x = conv2d_bn(x, 16, 3, 1, strides=(1, 1))
            x = layers.MaxPooling2D((3, 1), strides=(2, 1))(x)
            x = layers.Dropout(0.1)(x)
            x = layers.Flatten()(x)
            x = layers.Dense(8, activation='relu')(x)
            x = layers.Dense(self.num_classes)(x)


            return inputs,x


        else:
            print(json_model)
            with open(json_model) as json_file:
                json_savedModel= json_file.read()
                #config = json.load(json_file)

                return tf.keras.models.model_from_json(json_savedModel)
            #return tf.keras.models.model_from_config(config)

    def test(self, img_dir=None, arr_file=None, bed_file=None, model_dir="./model/",
             output_file="./predictions.tsv", colorN=3, imagemode=0):
        if (img_dir == None) == (arr_file == None):
            raise ValueError("img_dir and arr_file are mutually exclusive and one is required.")
        self._load_model(model_dir)
        self._model_dir = model_dir
        self._out_f = open(output_file, "wt")

        if img_dir != None:
            if os.path.isfile(os.path.join(img_dir, "IRFinder-IR-dir-AI.txt")):
                filesdir = "IRFinder-IR-dir"
            else:
                filesdir = "IRFinder-IR-nondir"
            arr_file = os.path.join(img_dir, filesdir + "-AI.txt")
            bed_file = os.path.join(img_dir, filesdir + ".txt")
            self._test_irfinder_result(arr_file,bed_file)

        elif arr_file != None:

            self._test_img_arr(arr_file, bed_file)

        self._out_f.close()


    def print_batsh_test_res(self, batch, batch_names, ori_res):

        batch = np.array(batch)
        pred = self.model["Model"].predict(batch)
        for i in range(len(batch_names)):
            score = tf.nn.softmax(pred[i])
            idx_max = np.argmax(score)
            pred_lab = self.model["Dataset"]["class_names"][idx_max]
            # line = "{}\t{}\t{}\n".format(batch_names[i],  np.max(score),pred_lab)
            # self._out_f.write(line)

            if pred_lab=="hIR":
                line = ori_res.readline().split("\t")
                # print(batch_names[i])
                # print(line[0]+":"+line[1]+"-"+line[2])
                while line[0] + ":" + line[1] + "-" + line[2] != batch_names[i]:
                    line = ori_res.readline().split("\t")

                #line[4] = str(np.max(score))
                line[4] = str(score[0].numpy())

                #line[-1]= pred_lab
                self._out_f.write(("\t").join(line) + "\n")
        return

    def _test_irfinder_result(self, arr_file, bed_file):

        arch = ImageArchive(None,arr_file )
        ori_res = open(bed_file , "rt")
        line = ori_res.readline().split("\t")
        line[4] = "CNN_IRscore"
        self._out_f.write(("\t").join(line))

        bar = Bar("Predicting images in {}".format(bed_file), max=len(arch))
        batch = []
        batch_names = []
        bar.start()

        for bed, arr in arch:
            if arr.is_valid:
                batch.append(getImageArrayFromRegion(arr.region, self.model["Image size"]))
                name=arr.name.split(":")
                pos=name[1].split("-")
                batch_names.append("{}:{}-{}".format(name[0],int(pos[0])+15,int(pos[1])-15))
                #batch=np.array(batch.tolist(), dtype="uint8")


                if len(batch_names) == 500:
                    self.print_batsh_test_res( batch, batch_names, ori_res)
                    batch = []
                    batch_names = []
                    bar.goto(arch.getIndex())
        if len(batch_names) !=0:
            self.print_batsh_test_res( batch, batch_names, ori_res)
        bar.finish()
        ori_res.close()
        return



    def _load_model(self, model_dir):
        print("Loading the best_model in {}".format(model_dir))
        model_info_file="{}/model_info.json".format(model_dir)
        model_file="{}/best_model.h5".format(model_dir)
        if not os.path.exists(model_info_file) or not os.path.exists(model_file):
            raise FileNotFoundError("Error! files model_info.json and best_model.h5 have to be in the model folder! ")
        with open(model_info_file, "rt") as fp:
            self.model=json.load(fp)
        self.model["Model"]=tf.keras.models.load_model(model_file)
        print("Done.")
        return


        return data






