import os,sys,glob
import loadgmt,kevent
import ssexyhelp
import MCstat
from optparse import OptionParser
from pylab import *
import mplrc
import numpy as np
from matplotlib.ticker import MultipleLocator

# ----------------------------------------------------------------------
def getStats(data,SpanJobAverage=-1,waverage,dim=0):
    ''' Get the average and error of all columns in the data matrix. '''

    #SpanJobAverage tells whether the estimator file contains the Bins
    #column added by the merge scripts during merge of a spanned job
    if ndim(data) > dim:
        numBins  = size(data,dim) 
        if SpanJobAverage != -1:
           dataAve  = sum(data[:,:-1]*data[:,-1][:,newaxis],dim)/(1.0*sum(data[:,-1])) 
           dataErr = std(data[:,:-1],0)/sqrt(len(data[:,0])-1.0) 
        elif waverage:
            vals = data[0::2]
            errs = data[1::2]
            weights = 1.0/(errs*errs)
            
            # Little hack of the average function to compute the weighted std
            # http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
            dataErr  = np.sqrt(1.0/numpy.average(1.0/weights, weights=weights))
            dataAve  = np.average(vals,weights=weights,dim) 
        else:
            dataAve  = average(data,dim) 
            bins = MCstat.bin(data) 
            dataErr = amax(bins,axis=0)
        
        
    return dataAve,dataErr

# -----------------------------------------------------------------------------
def getScalarEst(type,ssexy,outName,reduceFlag, skip=0):
    ''' Return the arrays containing the reduced averaged scalar
        estimators in question.'''

    fileNames = ssexy.getFileList(type)
    lAveraged = []
    oldFormat = False
    small = 100
    for i, fname in enumerate(fileNames):
        headers   = ssexyhelp.getHeadersFromFile(fname)
        #if len(headers)<small: small = len(headers)
        if '+/-' in headers:
            headers = headers[::2]
            print "Weigthed average merged files detected"
        if 'dnT' in headers:
            headers = headers[::2]
            oldFormat = True
            print "Old estimator file format detected"
        if 'Bins' in headers: 
           SpanJobAverage = headers.index('Bins')
           headers.pop(SpanJobAverage) 
           print "Normal average merged files detected"
        else: 
            SpanJobAverage = -1
        lAveraged += [SpanJobAverage]     
        small = min(small,len(headers))

    ave = zeros([len(fileNames),small],float)
    err = zeros([len(fileNames),small],float)
    for i,fname in enumerate(fileNames):
        # Compute the averages and error
        waverage = False
        headers   = ssexyhelp.getHeadersFromFile(fname)
        #if 'ZRatio' in headers:
        #    headers.pop(3)
        #    data = loadtxt(fname,ndmin=2,usecols=(0,1,2,4,5,6,7))[skip:,:]
        #else:
        data = loadtxt(fname,ndmin=2)[skip:,:]
        if 'dnT' in headers:
            headers = headers[::2]
        if '+/-' in headers:
            headers = headers[::2]
            waverage = True
        #if 'ZRatio' in headers:
        #    headers = np.array(headers)
        #    headers = np.delete(headers,3,axis=0)
        #    data = np.delete(data,3,axis=1)
        if oldFormat:
            data = data[:,::2]
        print fname
        print headers
        ave[i,:],err[i,:] = getStats(data,lAveraged[i],waverage)
    
    # output the estimator data to disk
    outFile = open('reduce-%s%s.dat' % (type,outName),'w');

    # the headers
    outFile.write('#%15s' % reduceFlag[0])
    for head in headers:
        outFile.write('%16s%16s' % (head,'+/-'))
    outFile.write('\n')

    # sort based on the reduced parameter
    Vals = []
    for i,fname in enumerate(fileNames):
        Vals.append(float(ssexy.params[ssexy.id[i]][reduceFlag]))
   
    print Vals
    lorder = range(len(fileNames))
    print lorder
    lorder = zip(*sorted(zip(Vals,lorder)))[1]  
    print lorder
    
    
    # record into a file
    for i in lorder:
        outFile.write('%16.8E' % float(ssexy.params[ssexy.id[i]][reduceFlag]))
        for j,h in enumerate(headers):
            outFile.write('%16.8E%16.8E' % (ave[i,j],err[i,j]))
        outFile.write('\n')
    outFile.close()

    return headers,ave,err;


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
              'a': r'N_A'}


    parser = OptionParser() 
    parser.add_option("-T", "--temperature", dest="T", type="float",
                      help="simulation temperature in Kelvin") 
    parser.add_option("-b", "--beta", dest="B", type="float",
                      help="number of particles") 
    parser.add_option("-v", "--reduce", dest="reduce",
                      choices=['r','x','y','T','b','a'], 
                      help="variable name for reduction [r,x,y,T,b,a]") 
    parser.add_option("-r", "--replica", dest = "r",type="int",
                      help="number of replica copies") 
    parser.add_option("-x", "--Lx", dest="x", type="int",
                      help="lattice width") 
    parser.add_option("-y", "--Ly", dest="y", type="int",
                      help="lattice height") 
    parser.add_option("-s", "--skip", dest="skip", type="int",
                      help="number of measurements to skip") 
    parser.add_option("-p", "--plot", help="plot a particular estimator", type="str")
    parser.set_defaults(plot=False)
    parser.set_defaults(skip=0)
    (options, args) = parser.parse_args() 
    if (not options.reduce):
        parser.error("need a correct reduce flag (-r,--reduce): [r,x,y,T,b,a]")
    # parse the command line options and get the reduce flag


    # Check that we are in the correct ensemble
    dataName,outName = ssexyhelp.getWildCardString(options) 
    
    ssexy = ssexyhelp.SSEXYHelp(options)
    ssexy.getSimulationParameters()
    if  (not ssexy.id): 
        print "No filenames detected satisfying those criteria"
        sys.exit()

    # We first reduce the scalar estimators and output them to disk
    head1,scAve1,scErr1 = getScalarEst('estimator',ssexy,outName,options.reduce, options.skip)
    if options.plot:
        if  not(options.plot in head1):
            print "Incorrect estimator to plot. Choose among: "
            print head1
        else:
            pindex =  head1.index(options.plot)
        rcParams.update(mplrc.aps['params'])
        # Get the changing parameter that we are plotting against
        param = []
        for ID in ssexy.id:
            param.append(float(ssexy.params[ID][options.reduce]))
        lab = ''
        options_dic = vars(options)
        for item in ['x','y','r','T','B']:
            if options_dic[item]:
               lab += r'%s=%s \,' %(parMap[item],options_dic[item])

        colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']

        figure(1,(8,6))
        connect('key_press_event',kevent.press)
        ax = subplot(111)
        
        errorbar(param, scAve1[:,pindex], yerr=scErr1[:,pindex],\
                marker='s',mec=colors[1],mfc=colors[1],\
                ls='None',capsize=4,label=r'$\mathrm{}%s$' %lab)
        xlabel(r'$\mathrm{%s}$' %parMap[options.reduce])
        ylabel(r'$\mathrm{E}$')
        legend(loc='best',frameon=False)
        #tight_layout()
        show()

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()

