import os,sys,glob
import loadgmt,kevent
import ssexyhelp
import MCstat
from optparse import OptionParser
import argparse
from pylab import *
import mplrc
import numpy as np
from matplotlib.ticker import MultipleLocator
from numpy import array
from uncertainties       import unumpy
from utils import *

def GetDicKey(dic, val=''):

    keys = []
    for key,value in dic.items():
            if  value==val: 
                keys.append(key)

    return keys

def ComputeS2(Zr):
    b  = -np.sum(unumpy.log(Zr))
    return b 

def ComputeMI(Zr):
    mid = len(Zr)//2
    b  = np.sum(unumpy.log(Zr[:mid]))
    t  = np.sum(unumpy.log(Zr[mid:]))
    return t-b
# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    # define the mapping between short names and label names 
    parMap = {'x': r'L_x',
              'y': r'L_y',
              'b': r'\beta',
              'T': r'T',
              'r': r'r',
              'a': r'N_A',
              'ran': r'ran',
              'det': r'det'}


    parser = argparse.ArgumentParser(description='Calculate entropy based on ratio of partition')
    parser.add_argument('fileNames', help='Scalar reduce files', nargs='+')
    parser.add_argument('--renyi',   action='store_true', default=False, help='Compute Renyi entropy')
    parser.add_argument('--mutual',  action='store_true', default=False,  help='Compute mutual information')
    parser.add_argument('--raw',     action='store_true', default=False,  help='Show raw reduce data')
    args = parser.parse_args() 

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    figure(1,(8,6))
    connect('key_press_event',kevent.press)
    nplots = args.renyi + args.raw + args.mutual
    iplot  = nplots
    if args.raw:
        ax  = subplot(nplots,1,iplot)
        iplot -= 1
        ax.set_xlabel(r'$\mathrm{%s}$' %parMap['a'])
        ax.set_ylabel(r'$\mathrm{\sigma_{[Z_{A_{i+1}}/Z_{A_{i}}]}}$')
        ax.set_xlim(-1,16)
    if  args.renyi:
        ax1 = subplot(nplots,1,iplot)
        iplot -= 1
        ax1.set_xlabel(r'$\mathrm{\Beta [K^{-1}]}$')
        ax1.set_ylabel(r'$\mathrm{S_2/L}$')
    if args.mutual:
        ax2 = subplot(nplots,1,iplot)
        ax2.set_xlabel(r'$\mathrm{\Beta [K^{-1}]}$')
        ax2.set_ylabel(r'$\mathrm{I_2/L}$')

    Betas = {}
    for fileName in args.fileNames:
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        b = parmap['b']
        if  not b in Betas:
            Betas[b] = {}

        if 'ran' in parmap: geom = 'ran'
        if 'det' in parmap: geom = 'det'
        #if 'A' in parmap: geom = parmap['A']
        #else:             geom = 'half'
        
        Betas[b][geom] = scalarhelp

    markers = {'ran': 's', 'det': 'o', 'full': 's', 'half': 'o', 'A': '>', 'B': '<'}
    dS2A, dS2Ae = {},{}
    dS2B, dS2Be = {},{}
    dS2F, dS2Fe = {},{}
    dMI,  dMIe  = {},{}
    i = 0
    j = 0
    for b, dic in Betas.items():
        for A, rscalar in dic.items():
        
            parmap  = rscalar.getParmap()
            Lx      = parmap['Lx']
        
            rscalar.loadData()
            rb, As = rscalar.getrParams()
            
            # Get Z ratios-----------------------------------------------

            # Simple one-sided ratio 
            nAr, dnAr  = rscalar.getAverages('nAred')
            nAe, dnAe  = rscalar.getAverages('nAext')
            inds_Al    = np.where(nAr==1.0)    
            inds_Ah    = np.where(nAe==1.0)

            print As[inds_Al]
            nnAr = unumpy.uarray(nAe[inds_Al],dnAe[inds_Al])
            nnAe = unumpy.uarray(nAr[inds_Ah],dnAr[inds_Ah]) 
            Zr_SN  = nnAr/nnAe 

            # Loop one-sided normalized Zr (intermediate version)
            LRatio,dLRatio  = rscalar.getAverages('LRatio')
            LRr   = unumpy.uarray(LRatio[inds_Al],dLRatio[inds_Al])
            LRe   = unumpy.uarray(LRatio[inds_Ah],dLRatio[inds_Ah])
            Zr_LN   = LRr/LRe
            
            # Loop two-sided normalized Zr (advanced version)
            ALRatio,dALRatio  = rscalar.getAverages('ALRatio')
            Zr_ALN = unumpy.uarray(ALRatio[inds_Al],dALRatio[inds_Al])

            if  j==0:
                norm = np.minimum(unumpy.std_devs(Zr_LN),unumpy.std_devs(Zr_ALN))
            # Plot raw data if needed -----------------------------------
            if  args.raw:
                
                #ax.plot(As[inds_Al],  unumpy.std_devs(Zr_SN)/norm, ls = '', marker='s',color=colors[(j)%len(colors)], label=r'$\mathrm{\beta=%0.3f \, SRT \, %s}$' %(b,parMap[A])) 
                #ax.plot(As[inds_Al],  unumpy.std_devs(Zr_LN)/norm,  ls = '', marker='o',color=colors[(j)%len(colors)], label=r'$\mathrm{\beta=%0.3f \, LRT \, %s}$' %(b,parMap[A])) 
                #ax.plot(As[inds_Al],  unumpy.std_devs(Zr_ALN)/norm, ls = '', marker='>',color=colors[(j)%len(colors)], label=r'$\mathrm{\beta=%0.3f \, ALRT \, %s}$' %(b,parMap[A])) 
                ax.errorbar(As[inds_Al],  unumpy.nominal_values(Zr_SN),  unumpy.std_devs(Zr_SN),  ls = '', marker='s',color=colors[(j)%len(colors)], label=r'$\mathrm{\beta=%0.3f \, SRT \, %s}$' %(b,parMap[A])) 
                ax.errorbar(As[inds_Al],  unumpy.nominal_values(Zr_LN),  unumpy.std_devs(Zr_LN),  ls = '', marker='o',color=colors[(j)%len(colors)], label=r'$\mathrm{\beta=%0.3f \, LRT \, %s}$' %(b,parMap[A])) 
                ax.errorbar(As[inds_Al],  unumpy.nominal_values(Zr_ALN), unumpy.std_devs(Zr_ALN), ls = '', marker='>',color=colors[(j)%len(colors)], label=r'$\mathrm{\beta=%0.3f \, ALRT \, %s}$' %(b,parMap[A])) 
                #if i==0: ax.errorbar(As[:],  unumpy.nominal_values(Zr), unumpy.std_devs(Zr), ls = '', marker=markers[A],color=colors[i%len(colors)], label=r'$\mathrm{\beta=%0.3f \, A=%s}$' %(b,A)) 
                #else:   
                #    if  (A=='A') or (A=='half'): ax.errorbar(As[:],  unumpy.nominal_values(Zr), unumpy.std_devs(Zr), ls = '', marker=markers[A],color=colors[i%len(colors)], label=r'$\mathrm{\beta=%0.3f}$' %b) 
                #    else:                        ax.errorbar(As[:],  unumpy.nominal_values(Zr), unumpy.std_devs(Zr), ls = '', marker=markers[A],color=colors[i%len(colors)]) 
            # Compute Renyi enetropy for subregions ---------------------
            if  A=='A'   : S2A = ComputeS2(Zr)/float(Lx)
            if  A=='B'   : S2B = ComputeS2(Zr)/float(Lx)
            if  A=='full': S2F = ComputeS2(Zr)/float(Lx)
            if  A=='half': 
                           S2A = ComputeS2(Zr[:len(Zr)//2])/float(Lx)
                           S2F = ComputeS2(Zr)/float(Lx)
                           S2B = S2A 
            j += 1
        print 'Lx=%2.0d b=%0.2f' %(int(Lx),float(b))
        if  args.renyi:
            if i==0:
                ax1.errorbar(b,S2A.n,S2A.s, marker='s',color= colors[0], label = r'$\mathrm{S_2A}$')
                if  not(A=='half'):
                    ax1.errorbar(b,S2B.n,S2B.s, marker='o',color= colors[1], label = r'$\mathrm{S_2B}$')
                ax1.errorbar(b,S2F.n,S2F.s, marker='>',color= colors[2], label =r'$\mathrm{S_{2A \bigcup B}}$')
            else:
                ax1.errorbar(b,S2A.n,S2A.s, marker='s',color= colors[0])
                if  not(A=='half'):
                    ax1.errorbar(b,S2B.n,S2B.s, marker='o',color= colors[1])
                ax1.errorbar(b,S2F.n,S2F.s, marker='>',color= colors[2])

            (dS2A[b],dS2Ae[b]) = (S2A.n, S2A.s)
            (dS2B[b],dS2Be[b]) = (S2B.n, S2B.s)
            (dS2F[b],dS2Fe[b]) = (S2F.n, S2F.s)

        if  args.mutual:
            MI = S2A+S2B-S2F #ComputeMI(Zr)/float(Lx)
            ax2.errorbar(b,MI.n,MI.s, marker='s',color=colors[0])
            (dMI[b],dMIe[b]) = (MI.nominal_value, MI.std_dev)
            print '     I/L  = ', MI 
        
        
        i += 1



    if  args.mutual:
        print '--------MI:'
        print 'dic04x04  = ', dMI, '\n'
        print 'dic04x04e = ', dMIe
    if  args.renyi:
        print '--------S2A:'
        print 'dicS2A04x04  = ', dS2A, '\n'
        print 'dicS2A04x04e = ', dS2Ae
        print '--------S2F:'
        print 'dicS2F04x04  = ', dS2F, '\n'
        print 'dicS2F04x04e = ', dS2Fe
    tight_layout()
    if args.raw:
        ax.legend(loc='best',frameon=False)
        #ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
        #      fancybox=True, shadow=True, ncol=5)
    if args.renyi:
        ax1.legend(loc='best',frameon=False)
        
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

