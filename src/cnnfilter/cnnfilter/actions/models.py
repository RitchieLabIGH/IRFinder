import os,pickle
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow.keras import  layers,  Sequential
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.callbacks import ModelCheckpoint
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import matplotlib.pyplot as plt
import json
import numpy as np
import re
from sklearn.metrics import classification_report, confusion_matrix,roc_auc_score,roc_curve,auc
from sklearn.preprocessing import LabelBinarizer
from progress.bar import Bar
from cnnfilter.utils.reader import ImageArchive
from cnnfilter.actions.extract import getImageArrayFromRegion
#from intron_scanner.actions.extract import get_data_from_array
from tensorflow.keras.models import Model
from sklearn.utils import class_weight
#from tensorflow.keras.applications import inception_v3
#from tensorflow.keras.applications.inception_v3 import conv2d_bn
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
            if self.colorN == 0:
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

    def _make_graphs(self, history, out_dir,sufix=""):
        acc = history.history['accuracy']
        val_acc = history.history['val_accuracy']
        loss=history.history['loss']
        val_loss=history.history['val_loss']
        epochs_range = range(len(acc))
        plt.figure(figsize=(8, 8))
        plt.subplot(1, 2, 1)
        plt.plot(epochs_range, acc, label='Training Accuracy')
        plt.plot(epochs_range, val_acc, label='Validation Accuracy')
        plt.legend(loc='lower right')
        plt.title('Training and Validation Accuracy')


        plt.subplot(1, 2, 2)
        plt.plot(epochs_range, loss, label='Training Loss')
        plt.plot(epochs_range, val_loss, label='Validation Loss')
        plt.legend(loc='upper right')
        plt.title('Training and Validation Loss')
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig("{}/accuracyplot{}.png".format(out_dir,sufix))



    def _get_array_train_test(self,data_names,data,arrposition):
        """
        seprarate dataset into 2 dataset for test and training, proportion given by self.validation_split
        :param data_names: classification array for each IR
        :param data: IR data
        :param arrposition: IR position
        :return:
        """


        class_names=np.sort(np.unique(data_names))
        val_ds = []
        val_target = []
        train_ds = []
        train_target = []
        val_arrpos = []
        train_arrpos = []
        counts=[[0,0] for _ in class_names]
        cl_idx=0
        data_names = np.array(data_names)
        arrposition = np.array(arrposition)
        data = np.array(data, dtype="uint8")
        for cl in class_names:

            np.random.seed(self.seed)

            classdata=data[data_names==cl]
            pos=arrposition[data_names==cl]
            #np.random.shuffle(classdata)
            randomize = np.arange(len(classdata))
            np.random.shuffle(randomize)
            classdata = classdata[randomize]
            pos = pos[randomize]

            val_num=int(round( len(classdata) * self.validation_split))

            val_ds.extend(classdata[:val_num,:])
            val_target.extend([np.argwhere(class_names==cl)[0,0] for i in range(val_num) ])
            val_arrpos.extend(pos[:val_num])
            # val_target.extend(np.ones(val_num,dtype=np.int8) *class_names)
            train_ds.extend(classdata[val_num:,:])
            train_arrpos.extend(pos[val_num:])
            train_target.extend([np.argwhere(class_names==cl)[0,0] for i in range(len(classdata)-val_num) ])

            counts[class_names.tolist().index(cl)][1]+=val_num
            counts[class_names.tolist().index(cl)][0]+=len(classdata)-val_num



        return  np.array(train_ds, dtype="uint8"), train_target,train_arrpos,np.array(val_ds, dtype="uint8") , val_target,val_arrpos, { "counts" : counts , "class_names" : class_names.tolist() }


    def train_from_array(self, img_dir="./output/", out_dir="./model/",
              img_size=256, batch_size=50,
              validation_split= 0.2, seed=123, epochs=10,
              num_threads = None, model_json= None , colorN=0,earlystop=0):
        '''
        Train your model from the array resultfile . bed tsv will be used for classes.
        You can give a custom model as json file.
        '''
        if num_threads != None:
            tf.config.threading.set_inter_op_parallelism_threads(num_threads)
        info={"Image directory" : img_dir ,
              "Output directory" : out_dir,
              "Validation split" : validation_split,
              "Epochs" : epochs,
              "Batch size" : batch_size,
              "Model json" : model_json,
              "Image size" : img_size,
              "Number of colors" : colorN,
              "Seed" : seed,
              "Threads" : num_threads,
              }
        if self.verbosity > 0:
            print("Starting the creation of a tensorflow model.\n\nParameters:")
            for k in info:
                print("  {} \t{}".format(k, info[k]))
        if os.path.isdir(out_dir):
            i=0
            out_dir=re.sub(r'\/$', '', out_dir)
            while os.path.isdir("{}_{}".format(out_dir, i)):
                i+=1
            out_dir="{}_{}".format(out_dir, i)
            print("Warning! Output folder existed. Changed to: {}".format(out_dir))
        os.makedirs(out_dir)

        ### Init self variables
        self.img_size=img_size;
        self.colorN=colorN;
        self.img_dir=img_dir;
        self.batch_size=batch_size;
        self.validation_split=validation_split;
        self.seed=seed;
        self.out_dir=out_dir;
        ### 
        from cnnfilter.actions.extract import Extractor
        extractor= Extractor( self.verbosity)
        if ".npy" in img_dir :
            if "train_val_data.npy" in img_dir:
                (
                train_ds, train_target, train_arrposition, val_ds, val_target, val_arrposition, dataset_info) = np.load(
                    img_dir, allow_pickle=True)
                data = val_ds
                labelarr = [ dataset_info["class_names"][a] for a in val_target]
                arrposition = val_arrposition
            else:
                print("Loading data from {}".format(img_dir))
                (data, labelarr, arrposition) = np.load(img_dir, allow_pickle=True)

                data = data.tolist()
                labelarr = labelarr.tolist()
                arrposition = arrposition.tolist()
        else:
            listfolds = os.listdir(img_dir)
            if "DATAbin.npy" in listfolds:
                print("Loading data from {}".format(os.path.join(img_dir, "DATAbin.npy")))
                (data, labelarr,arrposition) = np.load(os.path.join(img_dir, "DATAbin.npy"), allow_pickle=True)

                data = data.tolist()
                labelarr =labelarr.tolist()
                arrposition = arrposition.tolist()
            else:
                data=[]
                labelarr=[]
                arrposition=[]
                for fold in listfolds:
                    if os.path.isdir(os.path.join(img_dir,fold)):
                        print("extraction of "+fold,end=' ')
                        if False :#os.path.isfile(os.path.join(img_dir,fold, "IRFinder-IR-dir-AI.txt")):
                            filesAI="IRFinder-IR-dir-AI.txt"
                        else:
                            filesAI="IRFinder-IR-nondir-AI.txt"
                        #filesAI="IRFinder-IR-nondir-AI.txt"
                        print("extraction of "+os.path.join(img_dir, fold, filesAI),end=' ')
                        print(os.path.join(img_dir,fold,filesAI))
                        tmparr,tmplabelarr,tmparrposition=extractor.get_data_from_array( os.path.join(img_dir,fold,filesAI),os.path.join(img_dir,"..","classes",fold+".tsv") , subset=None, labels="lIR,mIR,hIR,noIR",img_size=img_size,colorN=colorN)
                        data.extend(tmparr)
                        labelarr.extend(tmplabelarr)
                        arrposition.extend(tmparrposition)
                        print("{} regions found".format(len(data)))
                np.save(os.path.join(img_dir, "DATAbin.npy"), (data, labelarr,arrposition))

        print("Imgages size {}\nNumber of regions : {}".format(data[0].shape,len(data)))

        train_ds, train_target,train_arrposition, val_ds, val_target,val_arrposition, dataset_info=self._get_array_train_test(labelarr,data,arrposition)
        np.save(os.path.join(out_dir,"train_val_data.npy"), (train_ds, train_target,train_arrposition, val_ds, val_target,val_arrposition, dataset_info))
        #(train_ds, train_target, val_ds, val_target,dataset_info) = np.load(os.path.join(img_dir,"testsave.npy"), allow_pickle=True)



        info["Dataset"]=dataset_info
        if self.verbosity > 0:
            print("Dataset information:\n\n\t{:10s}\t{:^20s}\t{:^20s}\n\t{}".format("Label","Training","Validation",("-"*60)))
            for cl in range(len(dataset_info["counts"])):
                tot_counts=np.sum(np.array(dataset_info["counts"][cl]))
                print("\t{:10s}\t{:^10d} ({:3.0f}%)\t{:^10d} ({:3.0f}%)".format(dataset_info["class_names"][cl],dataset_info["counts"][cl][0], (dataset_info["counts"][cl][0]*100)/tot_counts, dataset_info["counts"][cl][1], (dataset_info["counts"][cl][1]*100)/tot_counts  ))
        nflod=5
        scorevalacc =np.zeros((nflod,1))
        scoreacc = np.zeros((nflod, 1))
        scoreoutter = np.zeros((nflod, 1))
        scoreoutterbest = np.zeros((nflod, 1))
        epmaxval= np.zeros((nflod, 1),dtype=int)
        for i in range(nflod) :
            self.num_classes = len(dataset_info["class_names"])
            inputsM,outputsM=self._get_model(model_json)
            model= Model(inputs = inputsM,outputs = outputsM)

            model.compile(optimizer='adam',
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
                  metrics=['accuracy'])

            initweights=model.get_weights()

            if (self.verbosity > 0) & (i==0):
                print("\n\n")
                model.summary()
                print("\n\n")
            with open("{}/model.json".format(out_dir), "wt") as f :
                json.dump(json.loads(model.to_json()), f, indent=2)
            with open("{}/model_info.json".format(out_dir), "wt") as f :
                json.dump(info, f, indent=2)

            # save best model
            mc = ModelCheckpoint("{}/best_model{}.h5".format(out_dir,i), monitor='val_accuracy', mode='max', verbose=1,
                                 save_best_only=True)
            # simple early stopping
            if (earlystop >= 0) | (earlystop == -200):
                if earlystop == -200:  # default
                    earlystop = round(epochs / 10)
                es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=earlystop)
                cb = [es, mc]
            else:
                cb = [mc]

            print("Images size {}\nNumber of regions : {}".format(data[0].shape, len(data)))
            ttrain_ds, ttrain_target, ttrain_arrposition, tval_ds, tval_target, tval_arrposition, tdataset_info = self._get_array_train_test(train_target,train_ds,  train_arrposition)
            np.save(os.path.join(out_dir, "train_val_data_cross{}.npy".format(i)),
                    (ttrain_ds, ttrain_target, ttrain_arrposition, tval_ds, tval_target, tval_arrposition, tdataset_info))
            class_weights = class_weight.compute_class_weight('balanced',
                                                              np.unique(ttrain_target),
                                                               np.array(ttrain_target))
            # weightsamples = np.array([])
            # for pos_t in ttrain_arrposition:
            #     # Replacing "," , converting to lower and then splitting
            #     a = pos_t.split("\t")[1].split(":")[1].split("-")
            #     weightsamples = np.append(weightsamples, np.log(int(a[1]) - int(a[0])))

            model.set_weights( initweights)
            # tf.keras.utils.plot_model(
            #     model,
            #     to_file="model.png",
            #     show_shapes=False,
            #     show_dtype=False,
            #     show_layer_names=True,
            #     rankdir="TB",
            #     expand_nested=False,
            #     dpi=96,
            # )
            #tf.keras.utils.plot_model(model, to_file="model-graph.png", show_shapes=True)
            history = model.fit(ttrain_ds, np.array(ttrain_target),
                                validation_data=(tval_ds, np.array(tval_target)),
                                epochs=epochs,
                                callbacks=[es, mc],
                                class_weight=dict(enumerate(class_weights)),
                                #sample_weight=weightsamples,
                                batch_size=self.batch_size
                                )


            model.save("{}/model_{}.h5".format(out_dir,i))
            self._make_graphs(history, "{}/graphs/".format(out_dir,i),i)

            epmaxval[i]=np.argmax((history.history["val_accuracy"]))

            scorevalacc[i]=history.history["val_accuracy"][epmaxval[i,0]]
            scoreacc[i] = history.history["accuracy"][epmaxval[i,0]]
            # load from best model
            #model = tf.keras.models.load_model("{}/best_model{}.h5".format(out_dir, i))
            scoreoutter[i] = model.evaluate(np.array(val_ds),np.array(val_target))[1]
            model = tf.keras.models.load_model("{}/best_model{}.h5".format(out_dir, i))
            scoreoutterbest[i] = model.evaluate(np.array(val_ds),np.array(val_target))[1]
        print("scoreoutter")
        print(scoreoutter)
        print("scoreoutterbest")
        print(scoreoutterbest)

        bestfold=np.argmax(scoreoutter)
        #os.system("cp {}/best_model{}.h5 {}/best_model.h5".format(out_dir, bestfold,out_dir))
        # test if best model is no too bias addacp all test in function
        os.system("cp {}/model_{}.h5 {}/best_model.h5".format(out_dir, bestfold,out_dir))
        os.system("cp {}/train_val_data_cross{}.npy {}/train_val_data_cross.npy".format(out_dir, bestfold, out_dir))

        out_resumeall = open("rescrossval.txt", "at")
        if out_resumeall.tell() == 0:
            out_resumeall.write( "Path\tntrial\tbestepoque\tscoreacc,scorevalacc\tscoreoutter\n")
        for i in range(nflod):
            out_resumeall.write("{}\t{}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(os.path.basename(os.path.normpath(out_dir)),
                                                                      i,epmaxval[i,0],scoreacc[i][0],scorevalacc[i][0],scoreoutter[i][0]))
        out_resumeall.close()

        print("Done!")





    def test(self, img_dir=None, arr_file=None , bed_file=None,  model_dir="./model/", output_file="./predictions.tsv", colorN=0,imagemode=0):
        if (img_dir == None) == (arr_file == None):
            raise ValueError("img_dir and arr_file are mutually exclusive and one is required.")
        self._load_model(model_dir)
        self._model_dir = model_dir
        self._out_f= open(output_file, "wt")

        self._n_classes=len(self.model["Dataset"]["class_names"])
        self._prediction_counts=[[0 for _ in range(self._n_classes) ] for _ in range(self._n_classes+1) ]
        self._pred=[]
        self._target=[]
        self.countimg=0
        self.countimgtot = 0
        if img_dir != None:
            if   ".npy" in img_dir :
                self._out_f.write("Name\tPosition\tLabel\tPredicted_Label\tScore\n")
                self._test_imgvector_dir(img_dir,colorN=colorN)
            else:
                if os.path.isfile(os.path.join(img_dir, "IRFinder-IR-dir-AI.txt")):
                    filesdir = "IRFinder-IR-dir"
                else:
                    filesdir = "IRFinder-IR-nondir"
                arr_file = os.path.join(img_dir, filesdir + "-AI.txt")
                bed_file = os.path.join(img_dir, filesdir + ".txt")
                self._test_irfinder_result(arr_file, bed_file)
        elif arr_file != None:
#             if imagemode==0:
#                 self._test_arr(arr_file, bed_file)
#             else:
            self._test_img_arr(arr_file, bed_file,colorN=colorN)

        self._out_f.close()
        if np.sum(self._prediction_counts )!=0 :
            self._print_prediction_summary();
        return


    def _multiclass_roc_auc_score(self, average="macro"):#y_test, y_pred,
        fig, c_ax = plt.subplots(1,1, figsize = (12, 8))

        y_test=self._target
        y_pred=self._pred

        lb = LabelBinarizer()
        lb.fit(y_test)
        y_test = lb.transform(y_test)
        y_pred = lb.transform(y_pred)
        try:
            for (idx, c_label) in enumerate(self.model["Dataset"]["class_names"]): # all_labels: no of the labels, for ex. ['cat', 'dog', 'rat']
                fpr, tpr, thresholds = roc_curve(y_test[:,idx].astype(int), y_pred[:,idx])
                c_ax.plot(fpr, tpr, label = '%s (AUC:%0.2f)'  % (c_label, auc(fpr, tpr)))
            c_ax.plot(fpr, fpr, 'b-', label = 'Random Guessing')
            plt.legend()
            plt.savefig("ROCcurve.png")
        except:
            print("ROCcurve plot generation failed")

        return roc_auc_score(y_test, y_pred, average=average)

    def _print_prediction_summary(self):
        header="    |{:^15}".format("True\Pred");
        for i in range(self._n_classes):
            header+="|{:^15}".format(self.model["Dataset"]["class_names"][i]);
        header+="|{:>15} |".format("Total")
        print("\n\n   "+("-"*len(header)))
        print(header)
        print("   "+("-"*len(header)))
        totals_cols=[0 for _ in range(self._n_classes) ]
        for i in range(len(self._prediction_counts)):
            line="    "
            true_label=i < len(self.model["Dataset"]["class_names"])
            if true_label:
                line+="{:>15} ".format(self.model["Dataset"]["class_names"][i])
            else:
                line+="{:>15} ".format("NA")
            for j in range(len(self._prediction_counts[i])):
                line+="{:^16}".format(self._prediction_counts[i][j])
                totals_cols[j]+=self._prediction_counts[i][j]
            tot_row=sum(self._prediction_counts[i])
            if true_label or tot_row > 0 :
                line+="|{:^15}".format(tot_row)
                print(line)
        print("    "+("-"*len(header)))
        line="    {:>15} ".format("Total")
        for i in range(self._n_classes):
            line+="{:^16}".format(totals_cols[i])
        print(line+"\n\n")
        if not isinstance(self._target[0], int):
            if (len(self._target)!=len(self._pred) ) or ( ','.join(np.sort(self.model["Dataset"]["class_names"]))!= ','.join(np.sort(np.unique(np.array(self._target)))) ):
                print("New target label(s) in the test dataset, some complementary analysis will not be able to perform. ")
                conflabel = ["hIR", "mIR", "lIR", "noIR"]  # todo : remove just  for order do : np.unique(self._target))
                confmat = confusion_matrix(self._target, [self.model["Dataset"]["class_names"][i] for i in self._pred],
                                           labels=conflabel).astype(float)
                print(("\t{}").format(conflabel))
                out_resumeall = open("../resforall.txt", "at")
                if out_resumeall.tell() == 0:
                    out_resumeall.write("Path\t{}_to_noIR\t{}\t{}\t{}%_to_noIR\n".format("_to_noIR\t".join(conflabel), "total_IR_to_noIR", "%totalIR_to_noIR", "%_to_noIR\t".join(conflabel)))

                out_resumeall.write("{}\t{}\t{:.0f}\t{:.2f}\t".format(os.path.basename(os.path.normpath(os.getcwd())),"\t".join(["%.0f" % i for i in confmat[:, 3]]), sum(confmat[0:3, 3]),sum(confmat[0:3, 3])/sum(sum(confmat[0:3, ]))*100))
                for i in range(0, confmat.shape[0]):
                    confmat[i] /= sum(confmat[i]).astype(float) / 100
                    print(("{}\t{}").format(conflabel[i], confmat[i]))

                out_resumeall.write("{}\n".format("\t".join(["%.2f" % i for i in confmat[:, 3]])))
                out_resumeall.close()
        else:
            report=classification_report(self._target, self._pred, target_names=self.model["Dataset"]["class_names"])
            print(report)
            reportdict = classification_report(self._target, self._pred, target_names=self.model["Dataset"]["class_names"],
                                  output_dict=True)


            out_resume = open("../resvalidationforall.txt", "at")
            if out_resume.tell() == 0:
                out_resume.write("Path\tf1weighted score\t txt report flat\n")

            out_resume.write("{}\t{}\t{}\n".format(os.path.basename(os.path.normpath(os.getcwd())),
                                                               reportdict["weighted avg"]["f1-score"],
                                                               report.replace("\n","\t")))

            out_resume.close()




    def _print_prediction(self, names, pred):
        targ=True
        if ','.join(np.sort(self.model["Dataset"]["class_names"]))!= ','.join(np.sort(np.unique(np.array(names)[:,1]))):
            targ=False



        for i in range(len(names)):
            score = tf.nn.softmax(pred[i])
            idx_max=np.argmax(score)
            pred_lab=self.model["Dataset"]["class_names"][idx_max]
            line="{}\t{}\t{}\t{}\n".format(names[i][0], names[i][1], pred_lab , np.max(score))
            if names[i][1] in self.model["Dataset"]["class_names"]:
                self._prediction_counts[self.model["Dataset"]["class_names"].index(names[i][1])][idx_max]+=1
            else:
                self._prediction_counts[-1][idx_max]+=1
            self._pred.append(idx_max)
            #print("target {}".format(names[i][1]))
            #self._target.append(self.model["Dataset"]["class_names"].index(names[i][1]))
            if targ:
                self._target.append(self.model["Dataset"]["class_names"].index(names[i][1]))
            else:
                self._target.append(names[i][1])


            self._out_f.write(line)



        return

    def print_batsh_test_res(self, batch, batch_names, ori_res):
        """
        test a batch of data and write outputfile as subset of IRfinderresult but only IR confirm by the network

        """

        batch = np.array(batch)
        pred = self.model["Model"].predict(batch)
        for i in range(len(batch_names)):
            score = tf.nn.softmax(pred[i])
            idx_max = np.argmax(score)
            pred_lab = self.model["Dataset"]["class_names"][idx_max]
            if pred_lab=="hIR":
                line = ori_res.readline().split("\t")
                # print(batch_names[i])
                # print(line[0]+":"+line[1]+"-"+line[2])
                while line[0] + ":" + line[1] + "-" + line[2] != batch_names[i]:
                    line = ori_res.readline().split("\t")

                line[4] = str(score[0].numpy())
                #line[-1]= pred_lab
                self._out_f.write(("\t").join(line))
        return

    def _test_irfinder_result(self,arr_file, bed_file, colorN=0):
        """
        test IR regarding  IRfinder result file
        :param arr_file: array of data from the IRs generate by IRfinder
        :param bed_file: bed .tsv file result of IRfinder

        """
        arch = ImageArchive(None, arr_file)
        ori_res = open(bed_file, "rt")
        line = ori_res.readline().split("\t")
        line[4] = "CNN_IRscore"
        self._out_f.write(("\t").join(line))

        bar = Bar("Predicting images in {}".format(bed_file), max=len(arch))
        batch = []
        batch_names = []
        bar.start()

        for bed, arr in arch:
            if arr.is_valid:
                batch.append(getImageArrayFromRegion(arr.region, self.model["Image size"], colorN))
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

    def _test_imgvector_dir(self, img_dir):
        """
        test saved data from cnnfilter
        :param img_dir: *.npy
        :return:
        """

        batch=[]
        batch_names=[]
        print("Loading ",img_dir)


        if "train_val_data.npy" in img_dir :
            (train_ds, train_target,train_arrposition, val_ds, val_target,val_arrposition, dataset_info) = np.load(img_dir, allow_pickle=True)
            batch=val_ds
            names=val_target
            pos=val_arrposition


        else :
            (batch, names,pos) = np.load(img_dir,allow_pickle=True)
            batch = np.array(batch.tolist(), dtype="uint8")

            #batch = self.preprocess_data(batch)
        print("Done.")


        targ=True
        names=list(names)
        print(type(names[0]))
        if isinstance(names[0],str):
            if ','.join(np.sort(self.model["Dataset"]["class_names"]))!= ','.join(np.sort(np.unique(np.array(names)[:]))):
                targ=False


        #self._print_prediction( batch_names, self.model["Model"].predict(batch) )
        pred=self.model["Model"].predict(batch)


        for i in range(len(names)):
            score = tf.nn.softmax(pred[i])
            idx_max=np.argmax(score)
            pred_lab=self.model["Dataset"]["class_names"][idx_max]
            if type(names[i])==np.int64 :
                val_lab=self.model["Dataset"]["class_names"][names[i]]
            else :
                val_lab=names[i]
            line="{}\t{}\t{}\t{}\n".format(pos[i], val_lab, pred_lab , np.max(score))
            if val_lab in self.model["Dataset"]["class_names"]:
                self._prediction_counts[self.model["Dataset"]["class_names"].index(val_lab)][idx_max]+=1
            else:
                self._prediction_counts[-1][idx_max]+=1
            # if self.saveimg:
            #     self.savefileimgarray(pos[i],  val_lab, pred_lab, np.max(score), batch[i])
            self._pred.append(idx_max)
            #print("target {}".format(names[i][1]))
            #self._target.append(self.model["Dataset"]["class_names"].index(names[i][1]))
            if targ:
                self._target.append(self.model["Dataset"]["class_names"].index(val_lab))
            else:
                self._target.append(names[i])
            self._out_f.write(line)

    def savefileimgarray(self, pos,  val_lab, pred_lab, score,batch):

        if  True : #((val_lab == "hIR") & (pred_lab == "noIR")) | ((pred_lab == "hIR") & (val_lab == "noIR")):
            big=""
            # if (np.diff(np.fromstring( pos.split(':')[1], dtype=int, sep='-'))> 3000):
            #     big="large"
            # if ((val_lab == "noIR") & (pred_lab != "noIR")):
            #     self.countimg+=1
            #     if (self.countimg>200) :
            #         return
            # else :
            #
            #
            #  if big=="" and ((self.countimg>200) or (self.countimgtot > 700)):
            #     return
            # self.countimgtot+=1
            plt.figure(figsize=(8, 8))
            plt.plot(batch[:, 0, 0])
            plt.plot(batch[:, 1, 0])
            plt.ylim(0,256 )
            name="{}_{}_{}_{:.2f}".format(pos.replace('\t',"_"),  val_lab, pred_lab, score)
            plt.title("{} {} predict as {}, score {:.2f}".format(pos.replace("\t"," "),  val_lab, pred_lab, score))

            out_dir=os.path.join(self._model_dir,"imageval")
            os.makedirs(os.path.join(out_dir,"imageval"), exist_ok=True)
            extrapath=os.path.join(out_dir,val_lab+"_"+pred_lab)
            os.makedirs(extrapath, exist_ok=True)
            plt.savefig("{}/{}{}.jpg".format(extrapath,big,name.replace('/',"_")))
            plt.close()




    def _test_img_arr(self, arr_file, bed_file,colorN=0,onlylabelmatchingclass=False):
        arch = ImageArchive(bed_file, arr_file)
        bar  = Bar("Predicting images in {}".format(arr_file), max=len(arch))
        batch=[]
        batch_names=[]
        bar.start()
        for bed, arr in arch:
            if arr.is_valid:
                if (not onlylabelmatchingclass) | ( bed[-1]in self.model["Dataset"]["class_names"]):
                    batch_names.append([arr.name , bed[-1] ])
                    #print(batch_names)
                    #print(bed[-1])
                    batch.append(getImageArrayFromRegion(arr.region, self.model["Image size"] ,colorN))
                    if len(batch_names) == 50:
                        self._print_prediction(batch_names, self.model["Model"].predict(np.array(batch)) )
                        batch=[]
                        batch_names=[]
                        bar.goto(arch.getIndex())
        bar.finish()
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






