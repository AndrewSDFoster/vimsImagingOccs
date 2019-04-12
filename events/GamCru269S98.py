#!/usr/bin/python3
# Filename: analysis.py
# Produces plots of saturn occultations from .cub files

# prep:
# put .cub files in a directory called cubfiles, with a cubs.txt file listing all desired files
# mkdir -p ./figs/{spectra,frames}
# copy animationmaker.sh into ./figs
# ./analysis.py && cd figs/frames/ && ../animationmaker.sh && cd ../spectra && ../animationmaker.sh && cd ../.. && xviewer figs

# futz with parameters at the top of this file (caps-locked comments) to get good/desired plots
#TODO: use times instead of cube numbers

import pysis as ps
import numpy as np
import matplotlib.pyplot as plt

occname  = "GamCru269S98"
# SAMPLE CHANNEL FOR MOVIE
continuum = np.concatenate((np.arange(10), np.arange(98,120), np.arange(189,220)))
# STAR LOCATIONS TODO: automate by making aperture around brightest pixel
starpixx = (0,-1)
starpixy = (0,-1)
# CUBES TO ZOOM IN ON TODO: automate by finding region of high stddev?
zoomin   = 0
zoomax   = -1
# MOVIES?
movies	 = True
