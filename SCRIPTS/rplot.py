import argparse
import os,sys,glob
import loadgmt,kevent
import ssexyhelp
import MCstat
from optparse import OptionParser
from pylab import *
import mplrc
import numpy as np
from matplotlib.ticker import MultipleLocator
from uncertainties       import unumpy


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main():

    order = ['det', 'ran','Set','Lx','Ly','T','b','r', 'full','half', 'Ar', 'Ae']

    parser = OptionParser() 
    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('fileNames', help='Scalar estimator files', nargs='+')
    parser.add_argument('--estimator','-e', help='A list of estimator names that \
                        are to be plotted.', type=str)
    parser.add_argument('--xl', help='Lower x-axis limit', type=float)
    parser.add_argument('--xh', help='Higher x-axis limit', type=float)
    parser.add_argument('--yl', help='Lower y-axis limit', type=float)
    parser.add_argument('--yh', help='Higher y-axis limit', type=float)

    args = parser.parse_args()

    fileNames = args.fileNames
    
    if not fileNames:
        print "No files detected"
        sys.exit()

    headers = ssexyhelp.ScalarReduce(fileNames[0]).getHeaders()

    if (not args.estimator) or not(args.estimator in headers):
       print headers
       print "Specify a correct estimator"
       sys.exit()
        
    rvariable = headers[0]

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']

    figure(1,(8,6))
    connect('key_press_event',kevent.press)
    ax = subplot(111)
    for i,fileName in enumerate(fileNames):
        sReduce   = ssexyhelp.ScalarReduce(fileName)
        sReduce.loadData()
        t,x   = sReduce.getrParams()
        y,dy  = sReduce.getAverages(args.estimator)
        #x     = 1.0/np.array(x)
        print sReduce.paramMap
        if  args.estimator=='SS':
            x = 1.0/np.array(x)
        print sReduce.getTupleIdstr(order)
        if  i==0: y1 = unumpy.uarray(y,dy)
        else:     y2 = unumpy.uarray(y,dy)

        errorbar(x, y, dy,\
                marker='s',mec=colors[i],mfc=colors[i],color=colors[i],\
                ls='',capsize=4,label=r'$\mathrm{%s}$' %sReduce.getTupleIdstr(order))

    deltay = y2-y1  
    y  = np.absolute(unumpy.nominal_values(deltay), unumpy.std_devs(deltay))
    #ax.plot(x, y,  marker='s',mec=colors[i],mfc=colors[i],color=colors[i],\
    #     ls='',label=r'$\mathrm{%s}$' %sReduce.getTupleIdstr(order))
    if  args.estimator=='SS':
        plot([0,1],[2.0/np.pi*0,2.0/np.pi*1], 
             label=r'$\mathrm{Magic \, line}$')
    xlabel(r'$\mathrm{%s}$' %rvariable)
    ylabel(r'$\mathrm{%s}$' %args.estimator)
    #ylabel(r'$\mathrm{|\triangle E}|/\sigma_{\triangle E}$')
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if args.xl: xmin = args.xl
    if args.xh: xmax = args.xh
    if args.yl: ymin = args.yl
    if args.yh: ymax = args.yh
    xlim(xmin,xmax)
    ylim(ymin,ymax)
    #legend(loc='best',frameon=False)
    legend(loc='upper right',frameon=False)
    tight_layout()
    show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

