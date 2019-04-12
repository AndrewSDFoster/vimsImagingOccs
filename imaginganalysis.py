#!/usr/bin/python3
# Filename: analysis.py
# Produces plots of saturn occultations from .cub files

# prep:
# put .cub files in a directory called cubfiles, with a cubs.txt file listing all desired files
# mkdir -p ./figs/{spectra,frames} NOTE: there are more directories needed now, grep this file for savefigs to find them all

# futz with parameters at the top of this file (caps-locked comments) to get good/desired plots
#TODO: use times instead of cube numbers in plots

import pysis as ps
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from config import *

  #####################
  #     READ IN DATA  #
  #####################

# Get cube-size and allocate space
shape    = ps.CubeFile(cubdir+cubfiles[0]).shape
nspec    = shape[0]
height   = shape[1]
width    = shape[2]
maxpix   = np.max((height,width))
cubdata  = np.zeros((ncubs, nspec, height, width))

# Read in data

try:
  cubdata = np.load(occname+"data.npy")
except:
  print("Reading in data from %d cubes"%ncubs)
  for i in range(ncubs):
    cubdata[i] = ps.CubeFile(cubdir+cubfiles[i]).data
  np.save(occname+"data.npy", cubdata)

  #####################
  #  CREATE APERTURE  #
  #####################

# Square aperture set by starpix variable in config file
aper 	= np.zeros(cubdata.shape) 
aper[:,:,starpixy[0]:starpixy[1], starpixx[0]:starpixx[1]] += 1

  #####################
  #  CENTERING ALLOC  #
  #####################

# Allocate space for centering frames
indices = np.indices((height,width))
nconts  = len(continua)
col	= np.zeros((ncubs, 2, nconts))
bincols = np.zeros((ncubs, 2, nconts))
stretchcol	= np.zeros((ncubs, 2, nconts))
binstretchcols  = np.zeros((ncubs, 2, nconts))

  #####################
  #  GENERATE SPECTRA #
  #####################

#Create Spectrum Array
print("Generating spectra of pixels where x={%d,%d}, y={%d,%d}"%(starpixx[0], starpixx[1], starpixy[0], starpixy[1]))

# Background gradient / flatfield is corrected first
flatfield  = np.mean(cubdata[zoomax:], axis=0)
flatfield /= np.nanmedian(flatfield, axis=(1,2)).reshape((nspec,1,1))
flatcor    = cubdata / flatfield

# Now subtract out the sky level
specsky = flatcor.copy()
# Sky value is calculated for each frame and channel individually as the median of the out-of-aperture pixels
specsky[:,:,starpixy[0]:starpixy[1],starpixx[0]:starpixx[1]] = np.nan
specsky = np.nanmedian(specsky, axis=(2,3))
# Sky-value is subtracted from the frame
corspec = cubdata - specsky.reshape(ncubs,nspec,1,1)
# Spectra is the sum of pixels in the aperture
spectra = (corspec * aper).sum(axis=(2,3))
# Normalize each channel
normspectra = spectra / np.nanmedian(spectra[normclip:zoomin], axis=0)
# TODO use binned spectra for something?
#binspec     = np.zeros(binning, ncubs/binning, nspec)
#binspec     = normspectra.resize(binning, ncubs/binning, nspec).mean(axis=0)

  #####################
  # FOR EACH CONTINUA #
  #####################

# FOR LOOP OVER CONTINUA SEGMENTS
for i in range(nconts):
  print("continuum band %d"%i)
  # Set the continuum array of channels to the first entry in the continua tuple from config
  continuum = continua[i]

  #####################
  #     MAKE FRAMES   #
  #####################

  # Create frame, which is cubdata summed over desired wavelength channels
  frames  = cubdata[:,continuum].sum(axis=1)
  # Sum skyspectra over relevant wavelengths to get sky for each frame
  sky     = specsky[:,continuum].sum(axis=1)
  # Subtract sky/background for the centering frames
  frames -= sky.reshape((ncubs,1,1))

  #####################
  # MAKE SUMMED FRAME #
  #####################
  if backgroundcheck == True:
    # Sum frames and plt.imshow
    summedframe = frames.sum(axis=0)
    sfmin       = summedframe.min()
    sfmax       = summedframe.max()
    print("Summed frame for continuum %d ranges from %d to %d DN"%(i,sfmin,sfmax))
    plt.imshow(summedframe, interpolation='none', vmin=sfmin/2, vmax=-sfmin/2, cmap='copper')
    plt.colorbar(extend="both")
    plt.xlim(-0.5,width-0.5)
    plt.ylim(-0.5,height-0.5)
    plt.title("%s Summed Frame, %f-%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
    plt.savefig('figs/summedframecont%d.png'%(i))
    plt.clf()
    plt.close()


  #####################
  #CENTERING ALGORITHM#
  #####################

  # CURRENT CENTERING ALGORITHM: CENTER OF LIGHT

  # Find centers of frame using a center of light algorithm
  col[:,0,i]   = (frames*indices[0]*aper[:,0]).sum(axis=(1,2))/((frames*aper[:,0]).sum(axis=(1,2)))
  col[:,1,i]   = (frames*indices[1]*aper[:,0]).sum(axis=(1,2))/((frames*aper[:,0]).sum(axis=(1,2)))
  # Do it again where the light is stretched by some power
  stretchcol[:,0,i]   = (frames**stretch*indices[0]*aper[:,0]).sum(axis=(1,2))/((frames**stretch*aper[:,0]).sum(axis=(1,2)))
  stretchcol[:,1,i]   = (frames**stretch*indices[1]*aper[:,0]).sum(axis=(1,2))/((frames**stretch*aper[:,0]).sum(axis=(1,2)))
  # Make arrays of binned/smoothed centers to see if that reduces noise
  for bin in range(binning):
    bincols[:,:,i]          += np.roll(col[:,:,i]       , bin-np.int(binning/2), axis=0) / binning
    binstretchcols[:,:,i]   += np.roll(stretchcol[:,:,i], bin-np.int(binning/2), axis=0) / binning
 
#####################
#       PLOTS       #
#####################
  
  #####################
  #  LIGHTCURVE PLOTS #
  #####################

  # Make Lightcurves
  print("Making Lightcurve and Centering plots for Each Channel")
  plt.figure(num=None, figsize=(20,20), dpi=300, facecolor='w', edgecolor='k')

  # Unzoomed lightcurve
  plt.subplot(2,2,1)
  plt.plot(spectra[:,continuum].sum(axis=1))
  plt.title("%s %f:%f micron Lightcurve"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
  plt.xlabel("Cube Number")
  plt.ylabel("Summed Brightness")

  # zoomed lightcurve
  plt.subplot(2,2,2)
  plt.plot(np.arange(zoomin,zoomax), spectra[zoomin:zoomax,continuum].sum(axis=1))
  plt.title("%s %f:%f micron Lightcurve"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
  plt.xlabel("Cube Number")
  plt.ylabel("Summed Brightness")

  #####################
  # CENTERING PLOTS   #
  #####################
  
  # Plot centering
  plt.subplot(2,2,3)
  plt.plot(col[:,0,i], 'b.', label="Z Position")
  plt.plot(col[:,1,i], 'g.', label="X Position")
  plt.plot(bincols[:,0,i], 'c,', label="Z Binned")
  plt.plot(bincols[:,1,i], 'r,', label="X Binned")
  plt.legend(loc="upper left")
  plt.ylim(0,maxpix)
  plt.title("%s Center v Time, %f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
  plt.xlabel("Frame Number")
  plt.ylabel("Star's Central Pixel")

  # zoomed centering
  plt.subplot(2,2,4)
  plt.plot(np.arange(zoomin,zoomax), col[zoomin:zoomax,0, i], 'b.', label="Z Position")
  plt.plot(np.arange(zoomin,zoomax), col[zoomin:zoomax,1, i], 'g.', label="X Position")
  plt.plot(np.arange(zoomin,zoomax), bincols[zoomin:zoomax,0,i], 'c,', label="Z Binned")
  plt.plot(np.arange(zoomin,zoomax), bincols[zoomin:zoomax,1,i], 'r,', label="X Binned")
  plt.legend(loc="upper left")
  plt.ylim(0,(np.max((height,width))))
  plt.title("%s Center v Time - Zoomed, %f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
  plt.xlabel("Frame Number")
  plt.ylabel("Star's Central Pixel")

  plt.savefig("figs/centeringphtomcont%d-zoom.png"%i)
  plt.clf()
  plt.close()

  #####################
  # DIFF CENTER PLOTS #
  #####################

  # Plot differential centering
  print("Differential Centering Plots")
  for k in range(nconts):
    if i>k:
      plt.plot(col[:,0,i] - col[:,0,k], 'b.', label="Z Position")
      plt.plot(col[:,1,i] - col[:,1,k], 'g.', label="X Position")
      #plt.plot(np.arange(len(col)), np.mean(col[:zoomin,0,i] - col[:zoomin,0,k])*np.ones(len(col)), 'k')
      #plt.plot(np.arange(len(col)), np.mean(col[:zoomin,1,i] - col[:zoomin,1,k])*np.ones(len(col)), 'y')
      plt.legend(loc="upper left")
      plt.ylim(-5,5)
      plt.title("%s Center v Time %f:%f-%f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]], lambdas[continua[k][0]], lambdas[continua[k][-1]]))
      plt.xlabel("Frame Number")
      plt.ylabel("Differential Centering Position")
      plt.savefig("figs/differentialcentering%d-%d.png"%(i,k))
      plt.clf()
      plt.close()
      plt.plot(bincols[:,0,i] - bincols[:,0,k], 'b.', label="Z Binned")
      plt.plot(bincols[:,1,i] - bincols[:,1,k], 'g.', label="X Binned")
      plt.legend(loc="upper left")
      plt.ylim(-5,5)
      plt.title("%s Center v Time %f:%f-%f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]], lambdas[continua[k][0]], lambdas[continua[k][-1]]))
      plt.xlabel("Frame Number")
      plt.ylabel("Differential Binned Centering Position")
      plt.savefig("figs/binneddifferentialcentering%d-%d.png"%(i,k))
      plt.clf()
      plt.close()

  if stretch != 0:
    ######################
    #STRETCH CENTER PLOTS#
    ######################

    # Redo it for the stretched ones
    print("Making Stretched Centering Plots")
    plt.plot(stretchcol[:,0,i], 'b.', label="Z Position")
    plt.plot(stretchcol[:,1,i], 'g.', label="X Position")
    plt.plot(binstretchcols[:,0,i], 'c,', label="Z Binned")
    plt.plot(binstretchcols[:,1,i], 'r,', label="X Binned")
    plt.legend(loc="upper left")
    plt.ylim(0,maxpix)
    plt.title("%s Center v Time, %f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
    plt.xlabel("Frame Number")
    plt.ylabel("Star's Central Pixel")
    plt.savefig("figs/stretchedcenteringcont%d.png"%i)
    plt.clf()
    plt.close()
    plt.plot(np.arange(zoomin,zoomax), stretchcol[zoomin:zoomax,0,i], 'b.', label="Z Position")
    plt.plot(np.arange(zoomin,zoomax), stretchcol[zoomin:zoomax,1,i], 'g.', label="X Position")
    plt.plot(np.arange(zoomin,zoomax), binstretchcols[zoomin:zoomax,0,i], 'c,', label="Z Binned")
    plt.plot(np.arange(zoomin,zoomax), binstretchcols[zoomin:zoomax,1,i], 'r,', label="X Binned")
    plt.legend(loc="upper left")
    plt.ylim(0,(np.max((height,width))))
    plt.title("%s Center v Time - Zoomed, %f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]]))
    plt.xlabel("Frame Number")
    plt.ylabel("Star's Central Pixel")
    plt.savefig("figs/stretchedcenteringcont%d-zoom.png"%i)
    plt.clf()
    plt.close()
  
    ######################
    #STRETCH DIFF CENTER #
    ######################
  
    # Plot differential centering
    print("Stretched Differential Centering Plots")
    for k in range(nconts):
      if i>k:
        plt.plot(stretchcol[:,0,i] - stretchcol[:,0,k], 'b.', label="Z Position")
        plt.plot(stretchcol[:,1,i] - stretchcol[:,1,k], 'g.', label="X Position")
          #plt.plot(np.arange(len(stretchcol)), np.mean(stretchcol[:zoomin,0,i] - stretchcol[:zoomin,0,k])*np.ones(len(stretchcol)), 'k')
        #plt.plot(np.arange(len(stretchcol)), np.mean(stretchcol[:zoomin,1,i] - stretchcol[:zoomin,1,k])*np.ones(len(stretchcol)), 'y')
        plt.legend(loc="upper left")
        plt.ylim(-5,5)
        plt.title("%s Center v Time %f:%f-%f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]], lambdas[continua[k][0]], lambdas[continua[k][-1]]))
        plt.xlabel("Frame Number")
        plt.ylabel("Differential Centering Position")
        plt.savefig("figs/stretcheddifferentialcentering%d-%d.png"%(i,k))
        plt.clf()
        plt.close()
        plt.plot(binstretchcols[:,0,i] - binstretchcols[:,0,k], 'b.', label="Z Binned")
        plt.plot(binstretchcols[:,1,i] - binstretchcols[:,1,k], 'g.', label="X Binned")
        plt.legend(loc="upper left")
        plt.ylim(-5,5)
        plt.title("%s Center v Time %f:%f-%f:%f microns"%(occname, lambdas[continua[i][0]], lambdas[continua[i][-1]], lambdas[continua[k][0]], lambdas[continua[k][-1]]))
        plt.xlabel("Frame Number")
        plt.ylabel("Differential Centering Position")
        plt.savefig("figs/stretchedbinneddifferentialcentering%d-%d.png"%(i,k))
        plt.clf()
        plt.close()
  #####################
  # MOVIE             #
  #####################

  if movies:
    # Make all of the frames all-positive STRICTLY for plotting and stretching purposes
    # NOTE that this is the same shift for every frame to ensure that the animation stays faithful
    frames         -= frames.min()
    #gamma=1
    # Find minimum and maximum values for consistent colorscaling
    framin = frames.min()
    framax = frames.max()
    # Make Image movie
    print("Creating Movie of Continuum Wavelengths, value range %d, %d"%(framin, framax))
    for j in range(ncubs):
      plt.imshow(frames[j]**gamma, interpolation='None', vmin=framin**gamma, vmax=framax**gamma, cmap='copper')
      plt.colorbar()
      plt.plot(col[j][1][i],col[j][0][i], 'r+')
      plt.xlim(-0.5,width-0.5)
      plt.ylim(-0.5,height-0.5)
      plt.title("%s Frame: %04d"%(occname,j))
      plt.savefig('figs/framescont%d/frame%04d.png'%(i,j))
      plt.clf()
      plt.close()

  #####################
  # SPECTRUM PLOT     #
  #####################

# Make Spectrum plot
print("Making Spectrum Plots")
plt.imshow(normspectra.transpose(), vmin=0, vmax=2, interpolation='none', extent=[0,ncubs, lambdas[-1], lambdas[0]], aspect='auto', cmap='copper')
plt.colorbar()
for q in range(nconts):
  plt.plot((0,ncubs), (lambdas[continua[q][ 0]],lambdas[continua[q][ 0]]), 'c')
  plt.plot((0,ncubs), (lambdas[continua[q][-1]],lambdas[continua[q][-1]]), 'y')
plt.title("Spectrum %s"%occname)
plt.ylabel("Wavelength, microns")
plt.xlabel("Cube Number")
plt.savefig('figs/occultation.png')
plt.clf()
plt.close()
plt.imshow(normspectra[zoomin:zoomax].transpose(), vmin=0, vmax=2, interpolation='none', extent=[zoomin, zoomax, lambdas[-1], lambdas[0]], aspect='auto', cmap='copper')
plt.colorbar()
for q in range(nconts):
  plt.plot((zoomin, zoomax), (lambdas[continua[q][ 0]],lambdas[continua[q][ 0]]), 'c')
  plt.plot((zoomin, zoomax), (lambdas[continua[q][-1]],lambdas[continua[q][-1]]), 'y')
plt.title("Spectrum - Zoomed %s"%occname)
plt.ylabel("Wavelength, microns")
plt.xlabel("Cube Number")
plt.savefig('figs/occultation-zoom.png')
plt.clf()
plt.close()

  #####################
  # SPECTRUM MOVIE    #
  #####################

if movies:
  # Make Spectrum movie
  print("Making Spectrum Movie")
  for j in range(ncubs):
    plt.plot(lambdas, normspectra[j])
    for q in range(nconts):
      plt.plot((lambdas[continua[q][ 0]],lambdas[continua[q][ 0]]), (0,2), 'c')
      plt.plot((lambdas[continua[q][-1]],lambdas[continua[q][-1]]), (0,2), 'y')
    plt.title("%s Frame: %04d"%(occname,j))
    plt.ylim(0,2)
    plt.xlim(lambdas.min(),lambdas.max())
    plt.savefig("figs/spectra/spectra%04d.png"%(j))
    plt.clf()
    plt.close()

print("./analysis.py complete")
