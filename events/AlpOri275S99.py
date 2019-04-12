#!/usr/bin/python3

import numpy as np
occname  = "AlpOri275S99"
# SAMPLE CHANNEL FOR MOVIE
continuum = np.arange(98,120) #np.concatenate((np.arange(10), np.arange(98,120), np.arange(189,220)))
# STAR LOCATIONS TODO: automate by making aperture around brightest pixel
starpixx = (4,7)
starpixy = (0,3)
# clipping and binning
normclip =    0 # number of frames to clip for normalization of spectrum (if first few are bad)
binning  =  100 # temporal binning for binned spectra
# CUBES TO ZOOM IN ON TODO: automate by finding region of high stddev?
zoomin   =  800
zoomax   = 1000
# MOVIES?
movies	 = True
gamma    = 0.4
