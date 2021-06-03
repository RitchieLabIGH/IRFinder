import numpy as np
import os
import time
import json
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
from tensorflow.keras.preprocessing.image import array_to_img, img_to_array, load_img
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
from pysam import AlignmentFile
from pybedtools import BedTool, Interval
from progress.bar import Bar
from math import floor
from scipy import interpolate
import gzip
from scipy import signal

#for debug
import matplotlib.pyplot as plt
import time


from cnnfilter.utils.reader import ImageArchive


class Extractor:
    def __init__(self, verbosity=0):
        self.verbosity = verbosity;
        self.times={};
        return
    def _init_annotation(self, annotation_file):
        self.annotation_file=annotation_file
        if not os.path.exists(annotation_file):
            raise  Exception("File {} doesn't exists.".format(annotation_file))
        self.annot = BedTool(annotation_file)
        if not self.annotation_file.endswith("gz"):
            self.annot = self.annot.tabix()
        self.type_ann = self.annot[0][0].startswith("chr")
    
    def _tick(self, key, force=False):
        if self.verbosity > 1 or force:
            self.times[key]=time.perf_counter() #time.clock()
        
    def _tock(self, key, force=False):
        if self.verbosity > 1 or force:
            self.times[key]=time.perf_counter()-self.times[key] #time.clock()
    
    def _save_log_time(self, file=None):
        if self.verbosity > 1:
            if file == None:
                file="./{}.times.tsv".format(time.perf_counter()())
            with open(file, "w") as f:
                f.write("Process\ttime\n")
                for k in self.times:
                    f.write("{}\t{}\n".format(k, self.times[k]))
            
             
    

    def _add_flank(self, b, flank_size):
        '''
            Add a flank to the given region
        '''
        b.start -= flank_size-1
        b.end+= flank_size
        if b.start < 0 :
            b.start = 0
        return b

    def _filter_annot(self, b):
        '''
        Filter the annotation file to consider only exon annotation from protein coding genes
        '''
        if b[2] == 'exon':
            return str(b.attrs).find("protein_coding") != -1;
        return False

    def _message(self, message, verbosity=0):
        ''' 
        Print a message based on the verbosity
        '''
        if ( verbosity <= self.verbosity):
            if type(message) == type(""):
                print(' '+message, end='', flush=True)
            else:
                print(message, end='', flush=True)

    def _check_files(self, fnames):
        for fname in fnames:
            if not os.path.exists(fname):
                raise FileNotFoundError(fname)
        

    def _step(self, message):
        '''
        Start a new step printing the given message
        '''
        hbar= '-' * 40
        white = ' ' * round((40-len(message))/2)
        white_b= ' ' *10
        message="\n\n "+ white_b + hbar+ "\n"+white_b+ "|"+ white + message + white + "|\n "+white_b + hbar+"\n\n";
        self._message(message, 0)
    def get_data_from_array(self, array_file, bed_file, subset=None, labels=None ,img_size=256,colorN=3, output_type="t"):
        '''
        Generate images from an array file, output of the extract command, and bed file giving the 
        labels of the corresponding image in the last column.
        The labels that you want to save in the training folder as link has to be comma separated.
        The output type is a string containing the letters "t" and/or "o", where
        "t" -> training set ( the lines having the label in 'labels' )
        "o" -> others ( the lines not having the label in 'labels' )
        The images names will be: [bed_file_name]_[index].png
        
        '''
        self._tick("global", True)
        self._check_files([array_file, bed_file])
        experiment_name=os.path.split(bed_file)[1].split('.')[0]
        if labels != None:
            
            labels=labels.split(",")
            labels.sort()
        else:
            labels=[]
        has_subset=False
        names=set()
        arrays=[]
        arrlabel=[]
        arrposition=[]
        if subset != None:
            has_subset=True
            self._check_files([subset])
            with open(subset, "r") as sub:
                for line in sub:
                    names.add(line.split("\t")[0])
            
        do_training = "t" in output_type.lower()
        do_others = "o" in output_type.lower()
        if not (do_others or do_training):
            raise ValueError("The output type has to contain at least one of the letters 't' or 'o'")

        self._step("Parameters")
        self._message("Array file:\t{}\n".format(array_file))    
        self._message("Bed file:\t{}\n".format(bed_file))
        self._message("Image size:\t{}\n".format(img_size))        
        self._message("Produce target images:\t{}\n".format(do_training))
        self._message("Produce other images:\t{}\n".format(do_others))
        if has_subset:
            self._message("Subset:\t{} names\n".format(len(names)))
        if len(labels) > 0:
            self._message("Labels:\n")
            for lab in labels:
                self._message("\t\t{}\n".format(lab))
            self._message("\n")
        img_archive= ImageArchive(bed_file, array_file)
        s_name=bed_file.split("/")[-1].replace(".bed", "").replace(".tsv", "")
        bar=Bar("Generating the images for {}...".format(s_name), max=int(len(img_archive)))
        processed=0

        for bed, arr in img_archive:
            processed+=1
            if processed % 1000 == 0 :
                bar.goto(processed)
            label=bed[-1]
            is_training=labels.count(label) != 0
            if has_subset and not arr.name in names:
                continue
            if (is_training and do_training) :
                
                if arr.is_valid :
                    try:
                        # if "10:3779" in arr.name:
                        #     print("ln164")
                        #     1
                        arr.region=getImageArrayFromRegion(arr.region,img_size=img_size,colorN=colorN)
                        if arr.strand=="-" :
                            #arr.region=arr.region[::-1]
                            arr.region[:,:,0]=np.flip(arr.region[:,:,0],axis=0)

                    except Exception as e:
                        bar.clearln();
                        self._message("\n",0)
                        self._message("Line {} ".format(img_archive.getIndex()), 0)
                        self._message(e, 0)
                        self._message("\n",0)

                    arrays.append(arr.region)
                    arrlabel.append(label)
                    arrposition.append("{}\t{}".format(experiment_name,arr.name))#.replace(":","\t").replace("-","\t")))#labels.index(
                    #link_name="{}/{}/{}_{}_{}.png".format(img_dirs["training"],label,s_name, arr.name,  img_archive.getIndex() ) 

                    #link_name="{}/{}_{}_{}.png".format(img_dirs["other"], s_name, arr.name, img_archive.getIndex())

                
        bar.goto(processed)
        bar.finish()
        img_archive.close()
        self._tock("global", True)
        self._message("\n\nDone in {} seconds.\n".format(self.times["global"]))
        self._save_log_time()
        #print(arrays.length)
        return arrays,arrlabel,arrposition
        



def generateImagesArrayGreyFromRegion(region, img_size=None):


    '''
    Return the arrays composing an image from a given region squize or extend to img_size
    '''

    region_size = len(region)

    depth = max([sum(i) for i in region])
    if depth == 0:
        raise ArithmeticError("Error! trying to generate an image with zero depth.")
    reads_img = (np.array(region)[:, :] / depth) * 255


    if region_size < img_size:
        kindinterp = "nearest"
    else:
        kindinterp = "zero"  #"linear"


    f0 = interpolate.interp1d(np.arange(0, region_size-30), reads_img[15:-15,0], kind=kindinterp)
    f1 = interpolate.interp1d(np.arange(0, region_size-30), reads_img[15:-15,1], kind=kindinterp)

    reads_imgd1 = np.array([np.array(reads_img[0:15, 0])])
    reads_imgd1 = np.append(reads_imgd1, f0(np.arange(0, (img_size - 30)) * ((region_size - 31) / (img_size - 30))))
    reads_imgd1 = np.append(reads_imgd1, reads_img[-15:, 0])

    reads_imgd2 = np.array([np.array(reads_img[0:15, 1]+reads_img[0:15, 0])])


    reads_imgd2 = np.append(reads_imgd2, f1(np.arange(0, (img_size - 30)) * ((region_size - 31) / (img_size - 30)))+reads_imgd1[15:-15])

    reads_imgd2 = np.append(reads_imgd2, reads_img[ -15:,1]+reads_img[-15:, 0])

    reads_img2 = np.array([reads_imgd1,reads_imgd2])

    reads_img2 = np.expand_dims(np.rot90(np.round(reads_img2).astype("uint8"), k=3), axis=2)

    return reads_img2, None



def getImageArrayFromRegion(region, img_size=None,colorN=3):
    '''
    Return the numpy array representing the image from a given region
    '''
    img_sizey=img_size
    if colorN == 0:
        read_img, ann_img = generateImagesArrayGreyFromRegion(region,img_size)
        img_sizey=2

    if img_size==None:
        if ann_img==None:
            return read_img
        else:
            return np.concatenate((read_img, ann_img))
    else:
        if ann_img==None:
            if colorN ==0:
                return read_img

        else:
            
            return np.concatenate((img_to_array(array_to_img(read_img).resize((img_sizey,img_size-1)).transpose(1)), img_to_array(array_to_img(ann_img).resize((img_size,1))))) 


def getImageFromRegion(region, img_size=None,colorN=3):
    '''
    Return the image from a given region
    '''
    return array_to_img(getImageArrayFromRegion(region, img_size,colorN))

def saveImageFromRegion(region, file_name, img_size=None ,colorN=3):
    '''
    Save an image from a region. If img_size is None, the dimensions will be the ones of the region array.
    '''
    getImageFromRegion(region, img_size,colorN).save(file_name)


def generateArrayGenomicRegion(r_chr, r_start, r_end, bam, annot, target_size, max_rescale, exact):
        '''
        Generate an array representing a genomic region. 
        The region will be rescaled in according to the target_size, 
        but won't be higher than the max rescale.
        The output array contains
        [[ overlap_read, overlap_splicing, overlap_exon_annotation ] ,.. ] 
        
        '''
        
        region_size= r_end-r_start
        rescale=1
        if region_size > target_size:
            rescale = floor(region_size / target_size  )
            if rescale > max_rescale:
                rescale = max_rescale
            region_size= floor(region_size / rescale )  
        
        region  = [ [0,0,0] for i in range(region_size) ]
        # Encode the reads informations
        if bam.header.references[0].startswith("chr"):
            reads = bam.fetch(contig = "chr"+r_chr, start= r_start, end = r_end)
        else:
            reads = bam.fetch(contig = r_chr, start= r_start, end = r_end)
        
        for read in reads:
            if exact:
                if read.aend > r_end or read.pos < r_start:
                    continue
            blocks=read.get_blocks()
            for i in range(len(blocks)):
                if i > 0:
                    for j in range(floor((blocks[i-1][1]-r_start)/rescale)+1, floor((blocks[i][0]-r_start)/rescale)):
                        if j >= 0 and j < region_size:
                            region[j][1]+=1
                for j in range(floor((blocks[i][0]-r_start)/rescale), floor((blocks[i][1]-r_start)/rescale)+1):
                    if j >= 0 and j < region_size:
                            region[j][0]+=1
        if annot != None:                    
            if annot[0][0].startswith("chr"):
                bed_chr= "chr"+r_chr
            else:
                bed_chr=r_chr
            exons = annot.all_hits(Interval(bed_chr,r_start, r_end ))
            for ex in exons:
                if ex[2] == 'exon':
                    for i in range(floor((ex.start - r_start)/rescale),floor(( ex.end - r_start ) / rescale)):
                        if i >= 0 and i < region_size:
                            region[i][2]=1
        return region


