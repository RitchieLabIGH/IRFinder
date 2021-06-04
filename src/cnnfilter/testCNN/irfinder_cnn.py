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
import json
import gzip

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

    print("\t- test:  \t use a trained model to predict the class of a given set of images")
    print("\n")


def main(argv=None):
    '''Command line options.'''
    program_name = os.path.basename(sys.argv[0])
    title("CNN filter")
    if argv is None:
        argv = sys.argv
    try:
        program_name = program_name + " " + 'test process'
        parser = OptionParser(usage="Usage: %prog test ",
                              description="Test a neural network model on a given set of images")
        parser.add_option("-d", "--img-dir", dest="dir",
                          help="The directory containing the IRFinder results to predict. ",
                          metavar="DIR", type="string")
        
        parser.add_option("-m", "--model-dir", dest="model", metavar="DIR",
                          help="Folder containing the model. It has to contain the files best_model.h5 and model_info.json [default: %default]",
                          type="string")
        parser.add_option("-o", "--out", dest="out", help="Output file [default: %default]", metavar="FILE",
                          type="string")
        parser.add_option("-v", "--verbose", dest="verbose", action="count", help="Set tensorflow verbosity level")
        # set defaults
        parser.set_defaults(out="./predictions.tsv", verbose=0, colorN=3, dir=None, array=None, bed=None,
                            model="./model/")
        # process options
        (opts, _) = parser.parse_args(argv)
        if (opts.dir != None) == (opts.array != None):
            parser.error(
                "Parameter -d is required\n\nUse --help to get more information\n")
        from actions.models import IntronModeller
        modeller = IntronModeller(opts.verbose)
        modeller.test(opts.dir, opts.model, opts.out)
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