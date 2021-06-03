#!/usr/bin/env python3
# encoding: utf-8
'''
intron_scanner.intron_scanner -- shortdesc

intron_scanner.intron_scanner is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2020 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os

LIBRARY_LOCATION = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(LIBRARY_LOCATION)

from optparse import OptionParser


def title(title):
    bar = '-' * 50
    white = ' ' * round((50 - len(title)) / 2)
    white_b = ' ' * 5
    title = "\n " + white_b + bar + "\n" + white_b + "|" + white + title + white + "|\n " + white_b + bar + "\n\n"
    print(title)


def print_help():
    print("Usage: intron_scanner action [options]\n\nPossible actions:")
    print("\t- extract:\t exctract genomic regions (bed file) from a bam file as images.")
    print(
        "\t- a2i:    \t exctract images from an array produced by the extract command and organize in according to given labels")
    print("\t- train:  \t train a tensorflow model on a given set of images")
    print("\t- test:  \t use a trained model to predict the class of a given set of images")
    print("\n")


def main(argv=None):
    '''Command line options.'''
    program_name = os.path.basename(sys.argv[0])
    title("Intron scanner")
    action = "help"
    if len(sys.argv) > 1:
        action = sys.argv[1]

    if argv is None:
        argv = sys.argv[2:]
    try:

        if action == "train":
            program_name = program_name + " " + 'training process'
            parser = OptionParser(usage="Usage: %prog train ",
                                  description="Train a neural network model on a given set of images")
            parser.add_option("-d", "--img-dir", dest="dir",
                              help="The directory containing the images, divided in subdirectories in according to the classes. Example: ./training/ -> ./training/labelA/ ./training/labelB/ ",
                              metavar="DIR", type="string")
            parser.add_option("-o", "--out", dest="outdir", help="Output directory [default: %default]", metavar="DIR",
                              type="string")
            parser.add_option("-b", "--batch", dest="batch", help="Number of images to load. [default: %default]",
                              metavar="INT", type="int")
            parser.add_option("-s", "--image-size", dest="size", help="Images size [default: %default]", metavar="INT",
                              type="int")
            # parser.add_option("-c", "--color-number", dest="colorN",
            #                   help="Number of color's dimension : 1 -> grey (read)\n\t\t\t 2 ->(read and annotation ??),\n\t\t\t3 -> 3colors  [default: %default]",
            #                   metavar="INT", type="int")
            parser.add_option("-S", "--seed", dest="seed", help="Seed for the validation split [default: %default]",
                              metavar="INT", type="int")
            parser.add_option("-t", "--threads", dest="threads", help="Number of threads to use. [default: %default]",
                              metavar="INT", type="int")
            parser.add_option("-V", "--validation-split", dest="vsplit", metavar="FLOAT",
                              help="Fraction of the dataset to use for the validation [default: %default]",
                              type="float")
            parser.add_option("-e", "--epochs", dest="epochs", metavar="INT",
                              help="Number of training epoch [default: %default]", type="int")
            parser.add_option("-E", "--earlystop", dest="earlystop", metavar="INT",
                              help="Number of patience epoch for earlystop , -1 for no earlystop  [default: 0.1*epochs (10 percent of the total number of epochs]",
                              type="int")
            parser.add_option("-m", "--json-model", dest="model", metavar="FILE",
                              help="Load the tensorflow model from a json file [default: %default]", type="string")
            parser.add_option("-v", "--verbose", dest="verbose", action="count", help="Set tensorflow verbosity level")
            # set defaults
            parser.set_defaults(outdir="./model/", size=256, colorN=0, verbose=0, epoch=10, earlystop=-200, ext="png",
                                model=None, vsplit=0.20, batch=50, seed=123, threads=None, epochs=10)
            # process options
            (opts, _) = parser.parse_args(argv)
            required = "dir ".split()
            for r in required:
                if opts.__dict__[r] is None:
                    parser.error("Parameter %s required\n\nUse --help to get more information\n" % r)
            from cnnfilter.actions.models import IntronModeller
            modeller = IntronModeller(opts.verbose)
            # modeller.train(opts.dir, opts.outdir, opts.size, opts.batch, opts.vsplit, opts.seed, opts.epochs, opts.threads, opts.model, opts.colorN, opts.earlystop)
            modeller.train_from_array(opts.dir, opts.outdir, opts.size, opts.batch, opts.vsplit, opts.seed, opts.epochs,
                                      opts.threads, opts.model, opts.colorN, opts.earlystop)

        elif action == "test":
            program_name = program_name + " " + 'test process'
            parser = OptionParser(usage="Usage: %prog test ",
                                  description="Test a neural network model on a given set of images")
            parser.add_option("-d", "--img-dir", dest="dir",
                              help="The directory containing the images to predict. If they are in subdirectories, the subfolder name is used as true label.",
                              metavar="DIR", type="string")
            parser.add_option("-a", "--array-file", dest="array", metavar="FILE",
                              help="Use a file conaining the image information, produced by the extract process",
                              type="string")
            parser.add_option("-b", "--bed-file", dest="bed", metavar="FILE",
                              help="bed file associated to the array (-a). Can be a general tsv file. The last column is used as true label",
                              type="string")
            parser.add_option("-m", "--model-dir", dest="model", metavar="DIR",
                              help="Folder containing the model. It has to contain the files best_model.h5 and model_info.json [default: %default]",
                              type="string")
            # parser.add_option("-c", "--color-number", dest="colorN",
            #                   help="Number of color's dimension : 1 -> grey (read)\n\t\t\t 2 ->(read and annotation ??),\n\t\t\t3 -> 3colors  [default: %default]",
            #                   metavar="INT", type="int")
            parser.add_option("-o", "--out", dest="out", help="Output file [default: %default]", metavar="FILE",
                              type="string")
            parser.add_option("-v", "--verbose", dest="verbose", action="count", help="Set tensorflow verbosity level")
            # set defaults
            parser.set_defaults(out="./predictions.tsv", verbose=0, dir=None, array=None, bed=None,
                                model="./model/")
            # process options
            (opts, _) = parser.parse_args(argv)
            if (opts.dir != None) == (opts.array != None):
                parser.error(
                    "Parameters -a and -d are mutual exclusive and at least one is required\n\nUse --help to get more information\n")
            from cnnfilter.actions.models import IntronModeller
            modeller = IntronModeller(opts.verbose)
            modeller.test(opts.dir, opts.array, opts.bed, opts.model, opts.out)
        elif action == "help" or action == "-h" or action == "--help":
            print_help()
        else:
            raise ValueError("Action %s not recognized." % action)

    except Exception as e:
        print(program_name + ": " + repr(e) + "\n")
        print("\n\nFor help use --help\n\n")
        print(e)
        if __debug__:
            raise e
        return 2


if __name__ == "__main__":
    sys.exit(main())


#json,gzip,time,tensorflow,matplotlib,numpy,sklearn,re,progress