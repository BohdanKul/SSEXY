import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas
import matplotlib.pyplot as window
import time
import sys
import kevent

class Annotate(object):
    def __init__(self,Nx,Ny,ofname,sfname,filled):
        self.colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546"]
        if  ofname:
            print "Loading from file: ", ofname
            self.State=np.loadtxt(ofname,dtype=bool)
            (self.Ny, self.Nx) =self.State.shape  
        else:
            (self.Nx, self.Ny) = (Nx,Ny)
            if not(filled):   self.State = np.zeros((self.Ny,self.Nx),dtype=bool) 
            else:             self.State = np.ones((self.Ny,self.Nx),dtype=bool) 

        (self.ofname,self.sfname) = (ofname,sfname)
        (self.xlim, self.ylim)    = (1,1)
        (self.width, self.height) = (self.xlim/float(self.Nx),self.ylim/float(self.Ny))

        self.ax = plt.gca()
        self.ax.set_xlim(0,self.xlim) 
        self.ax.set_ylim(0,self.ylim) 
        
        self.Drag   = False
        self.start=0;

        self.rect  = Rectangle((0,0), self.xlim, self.ylim,alpha=0.1)
        if filled: self.Rects = np.ones((self.Ny,self.Nx),dtype=Rectangle)
        else:      self.Rects = np.zeros((self.Ny,self.Nx),dtype=Rectangle)
        
        self.rColors = {False: dict(facecolor="white"), 
                        True:  dict(facecolor=self.colors[1])}

        self.DrawState()
        self.ax.figure.canvas.draw()

        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        plt.connect('key_press_event',self.on_keypress)

    def on_keypress(self,event):
        
        #Save region A
        if  event.key == 'a':
            footer = "\n"
            for i in range(self.Ny):
                for j in range(self.Nx):
                    if  self.State[i,j]:
                        footer += str(j+i*self.Nx)+" "
            if (self.ofname or self.sfname):
                np.savetxt(self.sfname, self.State, fmt="%i",header="Nx=%i\nNy=%i\n" %(self.Nx,self.Ny),footer=footer) 
                print "State saved to: ", self.sfname
            else:
                print "No storage file specified"

        #Close window
        if  event.key == 'q': 
	    window.close()

        #Save regions A by increment 
        if  event.key == 't': 
            tfiles = []
            (lastx,lasty) = (-1,-1)
            counter   = 0 
            increment = 1
            tState  = np.zeros((self.Ny,self.Nx),dtype=bool)
            tfooter = "\n" 
            mainPart = self.sfname.split('.')[0]
            saveas = mainPart+('_a-%03d' %counter)+'.dat'
            tfiles.append(saveas)
            np.savetxt(saveas, tState, fmt="%i",header="Nx=%i\nNy=%i\n" %(self.Nx,self.Ny),footer=tfooter) 
            print "Saved to: ",saveas 
            for m in range(self.Ny):
                for n in range(self.Nx):
                    for i in range(m,self.Ny):
                        for j in range(self.Nx)[::-1]:
                            if  self.State[i,j]:
                                if  (j<lastx) or (i>lasty):
                                    (lastx,lasty) = (j,i)
                                    counter += 1
                                    tfooter += str(j+i*self.Nx)+" "
                                    tState[i,j]   = True
                                    if  (counter%increment == 0):
                                        mainPart = self.sfname.split('.')[0]
                                        saveas = mainPart+('_a-%03d' %counter)+'.dat'
                                        np.savetxt(saveas, tState, fmt="%i",header="Nx=%i\nNy=%i\n" %(self.Nx,self.Ny),footer=tfooter) 
                                        print "Saved to: ",saveas 
                                        tfiles.append(saveas)
            print "Generated files: "
            print ' '.join(tfiles[:-1])
            print ' '.join(tfiles[1:])
    def DrawState(self):
        for y in range(self.Ny):
            for x in range(self.Nx):
                self.Rects[y,x] = Rectangle((x*self.width,y*self.height),self.width,self.height,edgecolor=self.colors[3],**self.rColors[self.State[y,x]])
                self.ax.add_patch(self.Rects[y,x])
        row_labels = range(self.Nx)
        col_labels = range(self.Ny) 
        plt.xticks(np.arange(self.Nx)*self.width+self.width/2.0, row_labels)
        plt.yticks(np.arange(self.Ny)*self.height+self.height/2.0, col_labels)

    def XY2Rectangle(self,X,Y): 
        Xr = int(X//self.width)
        Yr = int(Y//self.height)
        return (Yr,Xr)
   

    def on_press(self, event):
        self.Drag = True
        self.x1 = event.xdata
        self.y1 = event.ydata
        if event.button == 1: self.newStatus = 1
        if event.button == 3: self.newStatus = 0
        #print 'button=%d, x=%d, y=%d, (ydata,xdata)=%s'%(
        #       event.button, event.x, event.y, self.XY2Rectangle(event.xdata, event.ydata),)
        self.SwitchStatus() 
        self.ax.figure.canvas.draw()

    def on_release(self, event):
        self.Drag = False
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.ax.figure.canvas.draw()

    def SwitchStatus(self):
        coords = (Yr,Xr) = self.XY2Rectangle(self.x1,self.y1) 
        if  (self.State[coords]!=self.newStatus):
            self.State[coords] = self.newStatus
            self.Rects[coords].update(self.rColors[self.newStatus])



    def on_motion(self, event):
        if self.Drag: 
            self.x1 = event.xdata
            self.y1 = event.ydata
            self.SwitchStatus() 
            if (time.clock() - self.start>0.15):
               self.ax.figure.canvas.draw()
               self.start=time.clock() 



# setup the command line parser options 
parser = argparse.ArgumentParser(description='Graphical editor A region')
parser.add_argument('--config','-c', help='Path to region A config file',type=str)
parser.add_argument('--save', '-n',help='New config file name', type=str)
parser.add_argument('--filled','-f',help='Should initial lattice state be filled', action="store_true")
parser.add_argument('--Lx',  '-x',help='Lattice width', type=int)
parser.add_argument('--Ly',  '-y',help='Lattice height',type=int)
parser.set_defaults(filled = False)
args = parser.parse_args()

if  (not(args.config)) and ((not(args.Lx)) or (not(args.Ly))):
    print "Specify lattice dimensions via Lx and Ly parameters or a config file"
    sys.exit(0) 

a = Annotate(args.Lx, args.Ly, args.config, args.save, args.filled)
plt.tight_layout()
plt.show()
