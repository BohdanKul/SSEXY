import os,sys,glob
import loadgmt,kevent
import ssexyhelp
import MCstat
from optparse import OptionParser
import argparse
from scipy import integrate
from pylab import *
import mplrc
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy import interpolate
from sets import Set
from simpserror import simps_error
from uncertainties   import unumpy
from utils import *
from scipy.optimize import curve_fit
from collections import OrderedDict
def GetDicKey(dic, val=''):

    keys = []
    for key,value in dic.items():
            if  value==val: 
                keys.append(key)

    return keys

def GetZ(Beta,E,dE, offset):        
    
    #spline = interpolate.splrep(Beta,E, w=1./dE)
    #Beta  = np.arange(min(Beta),max(Beta),0.001)
    

    ##E = np.random.normal(E,dE)
    #spline = interpolate.splrep(Beta,E)
    #Beta  = np.arange(0,max(Beta),0.01)
    #E   = interpolate.splev(Beta,spline)
    
    Sleng = len(Beta) - offset 
    Z  = np.zeros(Sleng)
    dZ = np.zeros(Sleng)
    for i in range(offset,len(Beta)):
        #print i, E[i], Beta[i]
        
        Z[i-offset]  = integrate.simps((-1)*E[:i],Beta[:i], even='first')
        if  i<3:
            dZ[i-offset] = 0    
        else:
            dZ[i-offset] = simps_error(dE[:i], Beta[:i], even='first')
    Beta = Beta[offset:]
    return unumpy.uarray(Z,dZ),Beta

def CheckUniq(l):
    l = list(l)
    size = len(l)
    lo = []
    for i in range(size):
        elem = l.pop()
        if  elem in l:
            #print "Not unique temperature ", elem
            lo.append(size - i -1)

    return lo

def CheckEqual(l1,l2):
    
    l1 = list(l1)
    l2 = list(l2)
    if  len(l1) > len(l2): 
        bigger  = Set(l1)
        smaller = Set(l2)                                    
    else:                 
        bigger  = Set(l2)
        smaller = Set(l1)                                    
    print bigger - smaller



# -----------------------------------------------------------------------------
# Renyi finite-beta, finite-L correction  
# -----------------------------------------------------------------------------
def RenyiCorrection(betas, A,B,C):
    return A*np.exp(-B*betas)-C

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    # define the mapping between short names and label names 
    parMap = {'x': r'L_x',
              'y': r'L_y',
              'b': r'\beta',
              'T': r'T',
              'r': r'r'}


    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    parser.add_argument('--energy',    action='store_true', default=False, help = "Show energy plots")
    parser.add_argument('--renyi',     action='store_true', default=False, help = "Show energy plots")
    parser.add_argument('--partition', action='store_true', default=False, help = "Show Renyi entropy plots plots")
    parser.add_argument('--mutual', action='store_true', default=False, help = "Show Renyi entropy plots plots")
    args = parser.parse_args() 

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    figure(1,(7,6))
    connect('key_press_event',kevent.press)
    print args.energy
    nplots = args.renyi + args.mutual + args.energy
    iplot  = nplots
    if  args.mutual:
        ax  = subplot(nplots,1,iplot)
        iplot -= 1
        ax.set_ylabel(r'$\mathrm{I/L}$')
        ax.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
    if  args.renyi:    
        ax1  = subplot(nplots,1,iplot)
        iplot -= 1
        ax1.set_ylabel(r'$\mathrm{(S_2-S_2[\beta_{max}])/L}$')
        ax1.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
        #ax1.set_yscale('log')
    if  args.energy:    
        ax2  = subplot(nplots,1,iplot)
        iplot -= 1
        ax2.set_xlabel(r'$\mathrm{%s[K^{-1}]}$' %parMap['b'])
        ax2.set_ylabel(r'$\mathrm{E/N [K^{-1}]}$')
       


    order = ['Lx','Ly','T','b']
    
    Lxs = {}
    offset = 1

    for fileName in args.fileNames:
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        Lx = int(parmap['Lx'])
        if  not Lx in Lxs:
            Lxs[Lx] = {}

        geom = GetDicKey(parmap)[0]
        r    = int(parmap['r'])
        Lxs[Lx][(geom,r)] = scalarhelp

    if  (args.renyi or args.mutual): 
        for Lx,dic in Lxs.items():
            if  not( ('full',1) in dic):
                print "One replica reduce file was not found for ", Lx
                sys.exit(0)
            else:
                t    = dic[('full',1)]
                del dic[('full',1)]
                Lxs[Lx] = OrderedDict(dic)
                Lxs[Lx][('full',1)] = t
    
    data = {}
    i = 0
    j = 0
    equalA = False
    shift = 0
    for Lx, dic in Lxs.items():
        for geom, rscalar in reversed(dic.items()):
            print geom
            
            rscalar.loadData()
            rb, Beta = rscalar.getrParams()
            E, dE    = rscalar.getAverages('ET')

            #if (1.0/Beta[0]) > 20:
            #    print "Inserting value at infinity"
            #    Beta = np.insert(Beta,0,0)
            #    E    = np.insert(E,0,0)
            #    dE   = np.insert(dE,0,0)
   
            mask = CheckUniq(Beta)
            if len(mask)>0:
                print "Deleting repeating betas: ", Beta[mask]
                Beta = np.delete(Beta,mask)
                E    = np.delete(E,mask)
                dE   = np.delete(dE,mask)
            
            Beta = Beta[shift:]
            if j>0: CheckEqual(Beta2,Beta)
            Beta2 = Beta
            
            Norm = int(Lx)**2 
            E    = E[shift:]*Norm
            dE   = dE[shift:]*Norm
            
            Z,BetaZ = GetZ(Beta,E,dE,offset)
        
            if args.energy:
                        if i==0: ax2.errorbar(Beta, E, dE, linestyle = 'None', color = colors[i], marker='o', label=r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))
                        else:    ax2.errorbar(Beta, E, dE, linestyle = 'None', color = colors[i], marker='o')
            if args.partition:
                        if i==0: ax2.plot(BetaZ, Z, linestyle = '-', color = colors[i], label=r'$\mathrm{r=%0.2d;\, A=%s}$' %(geom[1],geom[0]))
                        else:    ax2.plot(BetaZ, Z, linestyle = '-', color = colors[i])
            
            if  args.renyi:
                if  geom == ('full',1): Zf1 = Z 
                else:
                    S2  = -Z + 2*Zf1  
                    if   geom == ('half',2): S2A  = S2B = S2
                    elif geom == ('A',2):    S2A  = S2   #+ np.log(2)*(int(Lx)**2)*0.5  
                    elif geom == ('B',2):    S2B  = S2   #+ np.log(2)*(int(Lx)**2)*0.5
                    elif geom == ('full',2): S2AB = S2   #+ np.log(2)*(int(Lx)**2)

                    #FancyErrorbar(ax1,np.array(BetaZ),((S2-S2[-1])*1.0)/float(Lx),\
                    #      colors[i],((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))
                    S2n = unumpy.nominal_values(((S2-S2[-1])*1.0)/float(Lx))
                    S2d = unumpy.std_devs(((S2-S2[-1])*1.0)/float(Lx))
                    ax1.errorbar(np.array(BetaZ),S2n,S2d,
                          color = colors[i],label = ((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ (r'$\mathrm{S_{2%s}}$' %geom[0] )))

            j += 1
        #T, MI = loadtxt("MI_ED.dat", unpack=True)
        #ax.plot(1.0/np.array(T), MI/float(Lx),\
        #        marker='',color=colors[-1],\
        #        label=r'$\mathrm{ED}$')
            
            
       
        if args.renyi:
            (A,B,C) = (0.4,1,0.4)
            flshift = 50
            fhshift = -1# 895
            coeff, var_matrix = curve_fit(RenyiCorrection,np.array(BetaZ)[flshift:fhshift],S2n[flshift:fhshift],p0=(A,B,C))
            (A,B,C) = coeff
            S2pred  = RenyiCorrection(np.array(BetaZ),A,B,C)[flshift:fhshift]
            ax1.plot(BetaZ[flshift:fhshift], S2pred, linewidth = 2, label = "Fit")
            
        if args.mutual:
            I = S2A+S2B-S2AB
            FancyErrorbar(ax,np.array(BetaZ),(I*1.0)/float(Lx),\
            colors[i],r'$\mathrm{Integration}$')#Lx=%2.0d}$' %Lx)

        i += 1

    if  args.mutual:
        ax.legend(loc='best',frameon=False)
    if args.renyi:
       ax1.legend(loc='best',frameon=False)
    if args.energy:
       ax2.legend(loc='best',frameon=False)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

