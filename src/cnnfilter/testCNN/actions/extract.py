import numpy as np
from scipy.interpolate import interp1d





def getImageArrayFromRegion(region, img_size=None):
    '''
    Return the numpy array representing the image from a given region
    '''


    read_img, ann_img = generateImagesArrayGreyFromRegion(region,img_size)

    return read_img



def generateImagesArrayGreyFromRegion(region, img_size=None):


    '''
    Return the arrays composing an image from a given region
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


    f0 = interp1d(np.arange(0, region_size-30), reads_img[15:-15,0], kind=kindinterp)
    f1 = interp1d(np.arange(0, region_size-30), reads_img[15:-15,1], kind=kindinterp)

    reads_imgd1 = np.array([np.array(reads_img[0:15, 0])])
    reads_imgd1 = np.append(reads_imgd1, f0(np.arange(0, (img_size - 30)) * ((region_size - 31) / (img_size - 30))))
    reads_imgd1 = np.append(reads_imgd1, reads_img[-15:, 0])

    reads_imgd2 = np.array([np.array(reads_img[0:15, 1]+reads_img[0:15, 0])])


    reads_imgd2 = np.append(reads_imgd2, f1(np.arange(0, (img_size - 30)) * ((region_size - 31) / (img_size - 30)))+reads_imgd1[15:-15])

    reads_imgd2 = np.append(reads_imgd2, reads_img[ -15:,1]+reads_img[-15:, 0])

    reads_img2 = np.array([reads_imgd1,reads_imgd2])

    reads_img2 = np.expand_dims(np.rot90(np.round(reads_img2).astype("float32"), k=3), axis=2)

    return reads_img2, None

