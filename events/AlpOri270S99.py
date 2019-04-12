#!/usr/bin/python3

import numpy as np
occname  = "AlpOri270S99"
# What are the wavelengths?
lambdas  = np.loadtxt('/home/foster/data/vimsSaturnOccs/data/lambda.txt')[-256:]
# Location of Files
cubdir   = "cubfiles/"
cubfiles = open(cubdir+'/cubs.txt').read().splitlines()
ncubs    = len(cubfiles)
# SAMPLE CHANNEL FOR MOVIE
continua = (np.arange(10), np.arange(60,78), np.arange(100,120), np.arange(146,159), np.arange(189,220))
# STAR LOCATIONS TODO: automate by making aperture around brightest pixel
starpixx = (1,6)
starpixy = (1,5)
skipcol1 = True
# clipping and binning
normclip = 250 # number of frames to clip for normalization of spectrum (if first few are bad)
binning  =  10 # temporal binning for binned spectra
stretch  =   3 # stretch for center finding
# CUBES TO ZOOM IN ON TODO: automate by finding region of high stddev?
zoomin   = 500
zoomax   = 750
# MOVIES?
movies	 = False
gamma    = 0.4
