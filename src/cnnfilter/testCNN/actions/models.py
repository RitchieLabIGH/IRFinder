import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import tflite_runtime.interpreter as tflite
from scipy.special import softmax

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import json
import numpy as np

from utils.reader import ImageArchive
from actions.extract import getImageArrayFromRegion


class IntronModeller():
    def __init__(self, verbosity=0):
        if verbosity > 3:
            verbosity=3
        os.environ['TF_CPP_MIN_LOG_LEVEL'] ="{}".format(3-verbosity)
        self.verbosity = verbosity


    def test(self, img_dir=None, model_dir="./model/", colorN=3, imagemode=0):
        if (img_dir == None):
            raise ValueError("img_dir is required.")
        self._load_model(model_dir)
        self._model_dir = model_dir
        

        if img_dir != None:
            for filesdir in ["IRFinder-IR-dir", "IRFinder-IR-nondir"]:
                if os.path.isfile(os.path.join(img_dir, filesdir+"-AI.txt")):
                    output_file=os.path.join(img_dir, filesdir+"-val.txt")
                    self._out_f = open(output_file, "wt")
                    arr_file = os.path.join(img_dir, filesdir + "-AI.txt")
                    bed_file = os.path.join(img_dir, filesdir + ".txt")
                    print("Processing "+filesdir+".txt")
                    self._test_irfinder_result(arr_file,bed_file)
                    print("Done.")
                    self._out_f.close()

    def _predict(self, arr):
        self.model["Model"].reset_all_variables()
        self.model["Model"].set_tensor(self.model["InputDetails"][0]['index'], [arr])
        self.model["Model"].invoke()
        return self.model["Model"].get_tensor(self.model["OutputDetails"][0]['index'])[0]



    def _test_irfinder_result(self, arr_file, bed_file):
        ori_res = open(bed_file , "rt")
        line = ori_res.readline().split("\t")
        line[4] = "CNN_IRscore"
        self._out_f.write(("\t").join(line))
        ori_res.close()
        arch = ImageArchive(bed_file,arr_file )
        for bed, arr in arch:
            if arr.is_valid:
                pred= self._predict(getImageArrayFromRegion(arr.region, self.model["Image size"]))
                score = softmax(pred)
                idx_max = np.argmax(score)
                pred_lab = self.model["Dataset"]["class_names"][idx_max]
                if pred_lab=="hIR":
                    line = bed
                    line[4] = str(score[0])
                    self._out_f.write(("\t").join(line)+"\n")
        return



    def _load_model(self, model_dir):
        print("Loading the best_model in {}".format(model_dir))
        model_info_file="{}/model_info.json".format(model_dir)
        model_file="{}/best_model.tflite".format(model_dir)
        if not os.path.exists(model_info_file) or not os.path.exists(model_file):
            raise FileNotFoundError("Error! files model_info.json and best_model.h5 have to be in the model folder! ")
        with open(model_info_file, "rt") as fp:
            self.model=json.load(fp)
        self.model["Model"]=tflite.Interpreter(model_path=model_file)
        self.model["Model"].allocate_tensors()
        self.model["InputDetails"]=self.model["Model"].get_input_details()
        self.model["OutputDetails"]=self.model["Model"].get_output_details()
        print("Done.")
        return


        return data






