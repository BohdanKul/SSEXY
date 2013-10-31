import matplotlib

import os,sys
from pylab import *
import argparse
import numpy as np

import matplotlib.pyplot as plt; plt.rcdefaults()
import kevent
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib

w = 'white'
b = 'black'
vtypes = np.array([ [w,w,w,w],
                    [b,b,b,b],
                    [b,w,w,b],
                    [w,b,b,w],
                    [w,b,w,b],
                    [b,w,b,w]
                 ]);

offRcolor = 'black'
dRcolor   = 'none'
linkKwarg = {'color': 'black','linestyle':'--','linewidth':4}
loopKwarg = {'color': 'red'  ,'linestyle':'-'}

# Vertex constants
scale   = 1.0
rwidth  = 0.15*scale
rheight = rwidth/5.0*scale
cradius = rwidth/5.0
offset  = rwidth*0.5

# Algorithmic constants
Nspins     = 5
Noperators = 6

#Figure constants
xSpacing = rwidth/2.0*0
xL = 0
cellWidth = rwidth + xSpacing
xW = Nspins*cellWidth
xH = xL + xW

ySpacing = rwidth*1.25
yL = 0
cellHeight = 2*(2*cradius+offset)+rheight+ySpacing
yW = Noperators*cellHeight
yH = yL + yW

#cellHeight = yW/float(Noperators - 1.0)
#cellWidth  = xW/float(Nspins - 1.0)

def TransformGrid(cgrid):
    cgrid = cgrid[[1,2,3,0]]
    t = np.copy(cgrid[2])
    cgrid[2] = cgrid[1]
    cgrid[1] = t
    return cgrid

def GetOperXYCoords(coords):
    (leg,oper,spin) = coords
    if spin == 0: offset = np.array([0,0])
    else:         offset = np.array([-2*cradius,0])
    if (leg == 1) or (leg == 2):
       spin -=1
    xyCoords = ((spin+0.5)*cellWidth,((oper+0.5)*cellHeight)) + offset*0 
    return xyCoords

def GetLegsXYCoords(coords):
    (leg,oper,spin) = coords
    offsetgrid = TransformGrid(np.mgrid[-rwidth/2.0:rwidth/2.0:2j, 
                          rheight+offset:-rheight-offset:2j].reshape(2,-1).T)
    xyCoords = GetOperXYCoords(coords)+offsetgrid 
    return xyCoords[leg]

def DrawLink(ax,From,To,kwarg):
    xyFrom = GetLegsXYCoords(From) + (xL,yL)
    xyTo   = GetLegsXYCoords(To)   + (xL,yL)
    line = plt.Line2D((xyFrom[0],xyTo[0]),(xyFrom[1],xyTo[1]),**kwarg)
    ax.add_line(line)

def DrawPeriodicLink(ax,From,To,kwargs):
    (legF,operF,spinF) = From
    (legT,operT,spinT) = To
#    if not(((legF == 0) and (legT = 3)) or ((legF == 1) and (legT == 2))):
#       print "Wrongly connected periodic legs: from=%i to=%i" %(legF,legT)
    if (legF == 0) or (legF == 3): shift = -rwidth/2.5
    else:                          shift =  rwidth/2.5
    xyFrom = GetLegsXYCoords(From) + (xL,yL)
    xyTo   = GetLegsXYCoords(To)   + (xL,yL)
    line = plt.Line2D((xyFrom[0],xyFrom[0]+shift),(xyFrom[1],xyFrom[1]),**kwargs)
    ax.add_line(line)
    line = plt.Line2D((xyFrom[0]+shift,xyTo[0]+shift),(xyFrom[1],xyTo[1]),**kwargs)
    ax.add_line(line)
    line = plt.Line2D((xyTo[0],xyTo[0]+shift),(xyTo[1],xyTo[1]),**kwargs)
    ax.add_line(line)

def DrawSpin(ax,spin,state):
    rcenter = GetOperXYCoords((0,0,spin))
    cgrid   = np.mgrid[rcenter[0]-rwidth/2.0:rcenter[0]+rwidth/2.0:2j, 
                       rcenter[1]+rheight+offset:rcenter[1]-rheight-offset:2j].reshape(2,-1).T
    # add legs
    cgrid = TransformGrid(cgrid)
    xy = cgrid[0]
    if state == 1: fcol = 'black'
    else         : fcol = 'white'
    circle = mpatches.Circle(xy, cradius,ec="black", facecolor=fcol)
    ax.add_patch(circle)

    

def DrawVertex(ax,rcenter,vtype):
    #rcenter += (cellWidth/2.0,cellHeight/2.0)
    cgrid   = np.mgrid[rcenter[0]-rwidth/2.0:rcenter[0]+rwidth/2.0:2j, 
                       rcenter[1]+rheight+offset:rcenter[1]-rheight-offset:2j].reshape(2,-1).T
    # add legs
    cgrid = TransformGrid(cgrid)
    for (i,xy) in enumerate(cgrid):
        circle = mpatches.Circle(xy, cradius,ec="black",
                                facecolor=vtypes[vtype][i])
        ax.add_patch(circle)

    # add operator
    if vtype < 4 : rcolor = dRcolor
    else:          rcolor = offRcolor
    fancybox = mpatches.FancyBboxPatch(
            np.array(rcenter) - np.array((rwidth/2.0, rheight/2.0)), rwidth, rheight,
            boxstyle=mpatches.BoxStyle("Square", pad=0.02),fc=rcolor)
    ax.add_patch(fancybox)


def GetFileParams(fileName):
    [Nx,Ny,T] = fileName.split('-')[1:4] 
    return [int(Nx),int(Ny),float(T)]

def GetSpin(leg,bspins):
    '''
        Get the spin number based on the leg #
        and the spins associated with the operator
        attached to that leg
    '''    
    if (leg == 0) or (leg == 3): return bspins[0]
    else:                        return bspins[1]
    
def TypeLinks(bonds,links,nonIoper):
    '''
        Build a map LinksTypes which stores for a given leg #
        what type of link it is connected by. The choice is 
        periodic (True) or not periodic (False)
    '''    
    flinks = np.copy(links)
    mlinks = np.copy(links)
    LinksTypes = {}
    for i,link in enumerate(mlinks):
#        print '(i,link)', (i,link)
        if link == -1: continue
        # Get the associated spin
        # Below: spin = GetSpin(operator,leg,bond_sites[bond#])
        sFrom   = GetSpin(i%4,bonds[:,nonIoper[i//4]//2])
        # Set helper variables to their defaults
        tlinks  = np.copy(flinks)
        (tlinks[i],tlinks[link]) = (-1,-1)
        periodic = True
        # Check whether the spin is being acted at more than once
        for j,tlink in enumerate(tlinks):
            #Dont bother if the relevant link had been already visited
            if (tlink == -1): continue
            # Get the assosiated spin
            # Below: spin = GetSpin(operator,leg,bond_sites[bond#])
            
            tsFrom = GetSpin(tlink%4,bonds[:,nonIoper[tlink//4]//2])
            # In case we find another link acting on the same spin
            if  tsFrom == sFrom:
                # Check whether the new-link's legs have a higher p-coordinate
                if  max(max(link,i),max(tlink,j)) == max(tlink,i): 
                    # This signales that the current link is not a periodic one
                    periodic = False
                    break
        if (tlink != -1) or (link != -1):
           LinksTypes[i]    = periodic
           LinksTypes[link] = periodic
           # Mark as visited
           mlinks[link]= -1
           mlinks[i] = -1

    return LinksTypes

###########################################################################
def main():
    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('--output',  '-o', help= 'SSEXY output files', nargs='+')
    parser.add_argument('--bond','-b', help= 'Bond file', type=str)
    parser.add_argument('--meas','-m', help= 'Measurement # to display', type=int)
    

    args = parser.parse_args()
    meas = args.meas 
    
###########################################################################
    # Load data
    FTypes = ['bond', 'vertex','link', 'operator', 'loop','spin']
    FData  = (bonds,vertices,links,opers,loops,spins) = ([],[],[],[],[],[])

    for fname in args.output:
        for i,fType in enumerate(FTypes):
            if  fname.find(fType) != -1:
                ffile   = open(fname,'r')
                RawData = ffile.readlines()
                for line in RawData:
                    FData[i].append(np.fromstring(line,dtype=int,sep=' '))

    # Add a 0 element
    bonds  = np.array([np.insert(bonds[0][0::2],0,0),np.insert(bonds[0][1::2],0,0)])

   
    # Load simulation parameters
    (Nx,Ny,T) = GetFileParams(fname)
    Nspins     = Nx
    Noperators = len(vertices[meas])

    # Iniatiate graphics
    fig, ax = plt.subplots()
    plt.connect('key_press_event',kevent.press)

###########################################################################
    print
    print '----------------------------'
    print 'Spins'
    print spins[meas]
    # Display current spins state
    for (i,spin) in enumerate(spins[meas]):
        DrawSpin(ax,i,spin)

    
###########################################################################
    print
    print '----------------------------'
    print opers[meas]
    # Operator's coordiantes
    (p,s) = (0,0)
    # List of non-identity operators
    nonIoper = []
    # Display current operator state
    for oper in opers[meas]:
        if oper == 0: continue 
        nonIoper.append(oper)
        s = bonds[0][oper//2]
        xy = GetOperXYCoords((0,p,s))
        vtx = vertices[meas][p]
        DrawVertex(ax,xy,vtx-1)
        print '---------------\n', '(s,p) = ', (s,p), '\n(x,y) =', xy
        p += 1

    # Load links type map
    LinksTypes = TypeLinks(bonds,np.copy(links[meas]),nonIoper)
    
###########################################################################
    print
    print "-----------Links-----------"
    print links[meas]
    print LinksTypes
    for i,link in enumerate(links[meas]):
        # Skip if the relevant link has been already visit
        if link == -1: continue

        # Get the coordiantes of the link's origin
        (pFrom,legFrom)  = (i//4,i%4)
        sFrom   = GetSpin(legFrom,bonds[:,nonIoper[pFrom]//2])

        # Get the coordiantes of the link's end 
        (pTo,legTo)   = (link//4,link%4)
        sTo   = GetSpin(legTo,bonds[:,nonIoper[pTo]//2])

        # Draw link
        if LinksTypes[i]:   DrawPeriodicLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),linkKwarg)
        else:               DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),linkKwarg)

        # Mark link as "visited"
        links[meas][link]= -1
        links[meas][i] = -1

        # Print link's coordiantes as well
        print '===============\n', 'From (leg,p,s) ',(legFrom,pFrom,sFrom), '\nTo   (leg,p,s) ',(legTo,pTo,sTo)
        

###########################################################################
    print
    print "-----------Loops-----------"
    Llinks = zip( loops[meas], 
                  np.array(loops[meas])[range( 1, len(loops[meas]) ) + [0]]
                )
                
    Llinks.pop()
    print LinksTypes
    print Llinks 
    for i,link in enumerate(Llinks):

        # Get legs coordiantes 
        (pFrom,legFrom)  = (link[0]//4,link[0]%4)
        sFrom   = GetSpin(legFrom,bonds[:,nonIoper[pFrom]//2])

        (pTo,legTo)   = (link[1]//4,link[1]%4)
        sTo   = GetSpin(legTo,bonds[:,nonIoper[pTo]//2])

        # Print link's coordiantes as well
        print '===============\n', 'From (leg,p,s) ',(legFrom,pFrom,sFrom), '\nTo   (leg,p,s) ',(legTo,pTo,sTo)
        
        # Draw link
        Switch = i%2 == 0
        if Switch:
            DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),loopKwarg) 
            print 'Odd switch'
        else:
            if (pFrom == pTo) and (sTo != sFrom):  
                DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),loopKwarg) 
                print 'Switch'
            else:
                if LinksTypes[pFrom]:   
                    DrawPeriodicLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),loopKwarg)    
                    print 'Periodic'
                else:                   
                    DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo), loopKwarg)           
                    print 'Link'


###########################################################################
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.axis('equal')
    #plt.axis('off')

    plt.show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


