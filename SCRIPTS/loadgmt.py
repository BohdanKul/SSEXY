''' loadgmt - load a GMT colormap from a file.

Given a directory and file name, returns a colormap suitable for plotting in pylab which can
be broken up in various ways.
'''
from __future__ import division
#from pylab import *
import matplotlib.colors as col

# =======================================================================

def getcmap(dirName,fileName,GMTPath = None):
     ''' Load a user defined colormap in GMT format from a file.

	 Given a directory and colormap name, returns a pylab colormap object.'''

     import colorsys
     import numpy as N
     if type(GMTPath) == type(None):
         filePath = '/usr/local/share/cpt-city/' + dirName + '/' + fileName  + '.cpt'
     else:
         filePath = GMTPath + '/' + dirName + '/' + fileName +'.cpt'
     try:
         f = open(filePath)
     except:
         print 'file ',filePath, 'not found'
         return None

     #print filePath 
     lines = f.readlines()
     f.close()

     x = []
     r = []
     g = []
     b = []
     colorModel = 'RGB'
     for l in lines:
         ls = l.split()
         if l[0] == '#':
            if ls[-1] == 'HSV':
                colorModel = 'HSV'
                continue
            else:
                continue
         if ls[0] == 'B' or ls[0] == 'F' or ls[0] == 'N':
            pass
         else:
             x.append(float(ls[0]))
             r.append(float(ls[1]))
             g.append(float(ls[2]))
             b.append(float(ls[3]))
             xtemp = float(ls[4])
             rtemp = float(ls[5])
             gtemp = float(ls[6])
             btemp = float(ls[7])

     x.append(xtemp)
     r.append(rtemp)
     g.append(gtemp)
     b.append(btemp)

     nTable = len(r)
     x = N.array( x , N.float)
     r = N.array( r , N.float)
     g = N.array( g , N.float)
     b = N.array( b , N.float)
     if colorModel == 'HSV':
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
     if colorModel == 'HSV':
        for i in range(r.shape[0]):
            rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
     if colorModel == 'RGB':
         r = r/255.
         g = g/255.
         b = b/255.
     xNorm = (x - x[0])/(x[-1] - x[0])

     red = []
     blue = []
     green = []
     for i in range(len(x)):
         red.append([xNorm[i],r[i],r[i]])
         green.append([xNorm[i],g[i],g[i]])
         blue.append([xNorm[i],b[i],b[i]])
     colorDict = {'red':red, 'green':green, 'blue':blue}
     #cmap = cm.colors.LinearSegmentedColormap('cname',
     #         colorDict,rcParams['image.lut'])
     cmap = col.LinearSegmentedColormap('cname', colorDict) 
     return (cmap)

# =======================================================================
def getColorList(dirName,fileName,numColors):
	''' Return a list of colors in RGB tuple form. '''

        cmap = getcmap(dirName,fileName)
	colorlist = []
	for n in range(0,numColors):
		colorlist.append(cmap(1.0*n/(1.0*(numColors)))[0:3])

	return colorlist

# =======================================================================
def getMarkerList(numMarkers=0):
	''' Returns a list of markers used for plotting. '''
	markerList = ['o','s','^','v','>','<','p','D','*','H','1','2','3','4']
	#markerList = ['o','s','p','^','v','>','<','D','*','H','1','2','3','4']
	return markerList

# =======================================================================
def getMap(cmap,A):
	''' Return a scalar mappable instance useful for creating a colorbar. '''

	import matplotlib.cm as cm
	map = cm.ScalarMappable(cmap=cmap)
	map.set_array(A)

	return map
