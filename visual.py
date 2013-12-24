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
linkKwarg = {'color': 'black','linestyle':'--','linewidth':2}
loopKwarg = {'color': 'red'  ,'linestyle':'-'}

# Vertex constants
scale   = 1.0
rwidth  = 0.15*scale
rheight = rwidth/5.0*scale
cradius = rwidth/5.0
offset  = rwidth*0.5

# Algorithmic constants
Nspins     = 10
Noperators = 10

#Figure constants
xSpacing = rwidth/2.0*0
xL = 0
xW = 0
cellWidth = rwidth + xSpacing

ySpacing = rwidth*1.25
yL = 0
yW = 0
cellHeight = 2*(2*cradius+offset)+rheight+ySpacing

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
    if (min(legF,legT) == 0) or (min(legF,legT) == 3): shift = -rwidth/3.0
    else:                          shift =  rwidth/3.0
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
    #circle = plt.Circle(xy, cradius,ec="black", facecolor=fcol)
    circle = mpatches.Circle(xy, cradius,ec="black", facecolor=fcol)
    ax.add_patch(circle)

    

def DrawVertex(ax,rcenter,vtype,veroper):
    if veroper:
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
    
def DrawSpins(ax,spins):
    # Display current spins state
    for (i,spin) in enumerate(spins):
        DrawSpin(ax,i,spin)

def DrawOperators(ax,opers,bonds,vertices,veroper):
    # Operator's coordiantes
    (p,s) = (0,0)
    # Display current operator state
    for oper in opers:
        if oper == 0: continue 
        s = bonds[0][oper//2]
        xy = GetOperXYCoords((0,p,s))
        if veroper:
            vtx = vertices[p]
            DrawVertex(ax,xy,vtx-1,veroper)
        else:
            if oper%2 == 1: vtx = 5
            else:            vtx = 2
            DrawVertex(ax,xy,vtx,veroper)
        #print '---------------\n', '(s,p) = ', (s,p), '\n(x,y) =', xy
        p += 1

def DrawLinks(ax,links,bonds,opers,LinksTypes):
    for i,link in enumerate(links):
        # Skip if the relevant link has been already visit
        if link == -1: continue

        # Get the coordiantes of the link's origin
        (pFrom,legFrom)  = (i//4,i%4)
        sFrom   = GetSpin(legFrom,bonds[:,opers[pFrom]//2])

        # Get the coordiantes of the link's end 
        (pTo,legTo)   = (link//4,link%4)
        sTo   = GetSpin(legTo,bonds[:,opers[pTo]//2])

        # Draw link
        if LinksTypes[i]:   DrawPeriodicLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),linkKwarg)
        else:               DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),linkKwarg)

        # Mark link as "visited"
        links[link]= -1
        links[i] = -1

        # Print link's coordiantes as well
        #print '===============\n', 'From (leg,p,s) ',(legFrom,pFrom,sFrom), '\nTo   (leg,p,s) ',(legTo,pTo,sTo)
     
def DrawLoop(ax,Llinks,bonds,opers,LinksTypes):
    for i,link in enumerate(Llinks):

        # Get legs coordiantes 
        (pFrom,legFrom)  = (link[0]//4,link[0]%4)
        sFrom   = GetSpin(legFrom,bonds[:,opers[pFrom]//2])

        (pTo,legTo)   = (link[1]//4,link[1]%4)
        sTo   = GetSpin(legTo,bonds[:,opers[pTo]//2])

        # Print link's coordiantes as well
        #print '===============\n', 'From (leg,p,s) ',(legFrom,pFrom,sFrom),' leg# ',link[0], '\nTo   (leg,p,s) ',(legTo,pTo,sTo),' leg# ',link[1]
        
        # Draw link
        Switch = i%2 == 0
        if Switch:
            DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),loopKwarg) 
            #print 'Odd switch'
        else:
            if (pFrom == pTo) and (sTo != sFrom):  
                DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),loopKwarg) 
                #print 'Switch'
            else:
                if LinksTypes[link[0]]:   
                    DrawPeriodicLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo),loopKwarg)    
                    #print 'Periodic'
                else:                   
                    DrawLink(ax,(legFrom,pFrom,sFrom),(legTo,pTo,sTo), loopKwarg)           
                    #print 'Link'

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

meas = 0
###########################################################################
def main():
    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('--output',  '-o', help= 'SSEXY output files', nargs='+')
    parser.add_argument('--meas','-m', help= 'Measurement # to display', type=int)
    

    args = parser.parse_args()
    meas = args.meas 
    
###########################################################################
    # Load data
    titles = ['Initial state','After diagonal move','Loop itinerary','After off-diagonal move']
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
    xW = Nspins*cellWidth
    xH = xL + xW

    
    # Iniatiate graphics
    #plt.connect('key_press_even',kevent.press)
    bdone = False
    for meas in range(args.meas,args.meas+5,1):#range(len(spins)):
        fig, axes = plt.subplots(1,4,sharey = True,sharex=True,squeeze=True)
        meast = meas*3
        Noperators = 0
        for i,ax in enumerate(axes):
            # Generate a list of non-identity operators
            ax.set_title(titles[i])
            nonIoper = []
            for oper in opers[meast]:
                if oper == 0: continue 
                nonIoper.append(oper)
            
            if len(nonIoper)>Noperators: Noperators=len(nonIoper)
            yW = Noperators*cellHeight
            yH = yL + yW
        
            if i<2:
                    #Draw Spins
                    print '\n----------------------------\nSpins\n',spins[meas]
                    DrawSpins(ax,spins[meas])
                    
                    #Draw Operators
                    print '\n----------------------------\nOperators\n',nonIoper
                    print '\n----------------------------\nVertices\n', vertices[meas*2]
                    DrawOperators(ax,nonIoper,bonds,[],False)
            if i==2:
                    #Draw Spins
                    print '\n----------------------------\nSpins\n',spins[meas]
                    DrawSpins(ax,spins[meas])
                    
                    #Draw Operators
                    print '\n----------------------------\nOperators\n',nonIoper
                    DrawOperators(ax,nonIoper,bonds,vertices[meas*2],True)
                    print '\n----------------------------\nVertices\n', vertices[meas*2]
                    
                    # Load links type map
                    print links[meas]
                    print nonIoper
                    LinksTypes = TypeLinks(bonds,np.copy(links[meas]),nonIoper)
                    print LinksTypes

                    # Build loop tuple list
                    Llinks = zip( loops[meas], 
                                  np.array(loops[meas])[range( 1, len(loops[meas]) ) + [0]]
                                )
                    Llinks.pop()

                    #Draw Links
                    #print "\n-----------Links-----------\n",links[meas]
                    DrawLinks(ax,np.copy(links[meas]),bonds,nonIoper,LinksTypes)
                       
                    #Draw Loops
                    #print "\n-----------Loops-----------\n",Llinks
                    DrawLoop(ax,Llinks,bonds,nonIoper,LinksTypes)
            if i==3:
                    #Draw Spins
                    print '\n----------------------------\nSpins\n',spins[meas+1]
                    DrawSpins(ax,spins[meas+1])
                    
                    #Draw Operators
                    print '\n----------------------------\nOperators\n',nonIoper
                    DrawOperators(ax,nonIoper,bonds,vertices[meas*2+1],True)
                    print '\n----------------------------\nVertices\n', vertices[meas*2+1]
     
            #plt.axis('equal')
            ax.set_autoscaley_on(False)
            ax.set_autoscalex_on(False)
            ax.set_xlim([xL-cradius-rwidth/3.0,xL-cradius-rwidth/3.0+xW])#max(xW,yW)])
            ax.set_ylim([yL,yL+yW])#max(xW,yW)])
            meast += 1
            #plt.gca().set(adjustable='box')
            #if (prev != 0) and (max(xW,yW)>prev):
            #    axes[0].set_xlim([xL-cradius-rwidth/3.0,xL-cradius-rwidth/3.0+max(xW,yW)])
            #    axes[0].set_ylim([yL,yL+max(xW,yW)])
            #prev = 1    
            #break 


        #plt.tight_layout()
        #plt.subplots_adjust(left=0, right=1, bottom=0, top=1,wspace=0)
        #plt.ylim(yL,yH)
        #plt.xlim(yL,yH)
        #plt.ylim(xL-cradius-rwidth/3.0,xH)
        #plt.xlim(xL-cradius-rwidth/3.0,xH)
        #plt.axis('off')

        #plt.xlim(xL+xH,xH/2.0)
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        fig.set_size_inches(30,16)
        plt.savefig('Config/ssexy_'+str(meas))
        #plt.show()
        meas += 1
    print 'haha'
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()


