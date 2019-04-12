#!/usr/bin/python3

import numpy as np
occname  = "AlpOri268S98"
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

