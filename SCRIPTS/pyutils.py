"""
pyutils.py

Useful python routines that we have come up with over the years.

JBLux, AGDelma 2009
"""

import re
import gzip
import math
import sys

#TODO - add a dump file method that is complementary

# ----------------------------------------------------------------------------
def now():
    ''' Print the current date and time.'''
    from datetime import datetime
    print datetime.now().strftime("Last Updated %Y/%m/%d at %H:%M:%S")

# ----------------------------------------------------------------------------
def loadFile(filename,format='float',comment='#'):
    """ Loads filename into an array of data with same dimensions 
    as rows and columns of file and ignores comments.  
    Reads in all numbers as the specified format (string) 
    (int : integer, float : float, scientific : scientific notation).
    """
    
    data_array = []
    
    if format == 'int':
        type_converter = lambda x: int(x)
    elif (format == 'float' or format == 'scientific'):
        type_converter = lambda x: float(x) 
    
    if re.match('.*\.gz',filename):
        f = gzip.GzipFile(filename,'r')
    else:
        f = file(filename,'r')
       
    n=0 
    for line in f.readlines():
        if (line[0] != comment):# and (n<5000)):
            numbers = line.split()
            if len(numbers) == 1:
                data_array.append([type_converter(numbers[0])])
            else:
                data_array.append([type_converter(num) for num in numbers])
        n = n+1
    f.close()
    return data_array

# ----------------------------------------------------------------------------
def getDimensions(data):
    """ Get the dimensions of a matrix. """
    N = []
    N.append(len(data))
    try:
        N.append(len(data[0]))
    except:
        N.append(0)

    return N

# ----------------------------------------------------------------------------
def vector2Matrix(vec):
    """ Convert a numpy array to a matrix. """
    import numpy as n
    mat = n.zeros([1,len(vec)],float)
    mat[0,:] = vec
    return mat

# ----------------------------------------------------------------------------
def extrema(x, max = True, min = True, strict = False, withend = False):
    """
    This function will index the extrema of a given array x.
    
    Options:
        max         If true, will index maxima
        min         If true, will index minima
        strict      If true, will not index changes to zero gradient
        withend     If true, always include x[0] and x[-1]
    
    This function will return a tuple of extrema indexies and values
    """
    
    # This is the gradient
    from numpy import zeros
    dx = zeros(len(x))
    from numpy import diff
    dx[1:] = diff(x)
    dx[0] = dx[1]
    
    # Clean up the gradient in order to pick out any change of sign
    from numpy import sign
    dx = sign(dx)
    
    # define the threshold for whether to pick out changes to zero gradient
    threshold = 0
    if strict:
        threshold = 1
        
    # Second order diff to pick out the spikes
    d2x = diff(dx)
    
    if max and min:
        d2x = abs(d2x)
    elif max:
        d2x = -d2x
    
    # Take care of the two ends
    if withend:
        d2x[0] = 2
        d2x[-1] = 2
    
    # Sift out the list of extremas
    from numpy import nonzero
    ind = nonzero(d2x > threshold)[0]
    
    return ind, x[ind]


# ----------------------------------------------------------------------------
def average(data,dim=1):

    """ Average a data matrix over the dimension given in dim. 
    
    Only works with 2d lists. dim=0 means average over rows, while dim=1 means
    average over columns (default)."""

    # We first determine whether or not we are dealing with a 1d or 2d python list
    N = getDimensions(data)

    if N[dim] == 0:
        raise NameError('Dimension to average over does not match with dimension of array!')

    ave = []
    if dim == 0:
        for i in range(N[0]):
            temp = 0.0
            for j in range(N[1]):
                temp += data[i][j]
            ave.append(temp/(1.0*N[1]))
    elif dim == 1:
        for j in range(N[1]):
            temp = 0.0
            for i in range(N[0]):
                temp += data[i][j]
            ave.append(temp/(1.0*N[0]))

    return ave

# ----------------------------------------------------------------------------
def error(data,dim=1):
    """ Computes the error in a data matrix over the dimension given in dim. 
    
    Only works with 2d lists. dim=0 means average over rows, while dim=1 means
    average over columns (default)."""

    # we get the dimensions of our data matrix
    N = getDimensions(data)

    if N[dim] == 0:
        raise NameError('Dimension to average over does not match with dimension of array!')

    # Create the data^2 matrix
    data2 = []
    for i in range(N[0]):
        data2.append([d*d for d in data[i]])

    # now we get the averages of data and data2
    ave  = average(data,dim)
    ave2 = average(data2,dim)

    # try and get the error matrix
    norm = N[not dim] - 1
    if norm == 0:
        err = [0.0 for n in range(len(ave))]
    else:
        err = [math.sqrt(abs(ave2[n]-ave[n]*ave[n])/(1.0*norm)) for n in range(len(ave))]

    return err

# ----------------------------------------------------------------------------
def bootstrap(data,dim=1):
    """ Computes the bootstrap error in a data matrix over the dimension given in dim. 
    
    Only works with 2d lists."""

    import random

    # the number of bootstrap samples
    numBoot = 500

    # we get the dimensions of our data matrix
    N = getDimensions(data)

    # make sure we have at least 2 elements of data
    if N[not dim] > 1:
        sampleAve = []

        # we need to compute things differently based on which dimension
        # the average is being performed over

        # Average over rows (non-standard)
        if dim == 0:
            for i in range(N[0]):
                ave = []
                for b in range(numBoot):
                    temp = 0.0
                    for j in range(N[1]):
                        k = random.randint(0,N[1]-1)
                        temp += data[i][k]
                    ave.append(temp/(1.0*N[1]))
                sampleAve.append(ave)

        # average over columns (standard)
        elif dim == 1:
            for j in range(N[1]):
                ave = []
                for b in range(numBoot):
                    temp = 0.0
                    for i in range(N[0]):
                        k = random.randint(0,N[0]-1)
                        temp += data[k][j]
                    ave.append(temp/(1.0*N[0]))
                sampleAve.append(ave)

        sampleAve2 = []
        for i in range(len(sampleAve)):
            sampleAve2.append([d*d for d in sampleAve[i]])

        bootAve = average(sampleAve,0)
        bootAve2 = average(sampleAve2,0)

        norm = 1.0*numBoot/(1.0*(numBoot-1))
        bootErr = [math.sqrt(norm*abs(bootAve2[n]-bootAve[n]*bootAve[n])) for n in range(N[dim])]

    else:
        bootErr= [0.0 for n in range(N[dim])]

    return bootErr

# ----------------------------------------------------------------------------
def isList(x):
    try:
        len(x)
        return True
    except:
        return False
#    return hasattr(x,'__iter__')

# ----------------------------------------------------------------------------
def DedekindEta(tau):
    ''' Compute the Dedekind Eta function for imaginary argument tau. 
        See: http://mathworld.wolfram.com/DedekindEtaFunction.html. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    qbar = mp.exp(-2.0*mp.pi*tau)
    
    return mp.nsum(lambda n: ((-1)**n)*(qbar**((6.0*n-1.0)*(6.0*n-1.0)/24)),[-mp.inf,mp.inf])

    
# ----------------------------------------------------------------------------
def dDedekindEta(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        See: http://mathworld.wolfram.com/DedekindEtaFunction.html. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    qbar = mp.exp(-2.0*mp.pi*tau)

    res = mp.nsum(lambda n: ((-1)**n)*((6.0*n-1.0)**2)*(qbar**((6.0*n-1.0)*(6.0*n-1.0)/24.0)),[-mp.inf,mp.inf])
    return -(mp.pi/ 12.0)*res
## ----------------------------------------------------------------------------
def DedekindEtaA(tau):
    ''' Alternative definition of the Dedekind Eta function in terms 
        of the Jacobi theta function. 10x faster computation. 
        http://functions.wolfram.com/EllipticFunctions/DedekindEta/27/01/02/'''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    
    return mp.cbrt(0.5*mp.jtheta(1,0,mp.exp(-mp.pi*tau),1)) 

 #----------------------------------------------------------------------------
def dDedekindEtaA(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        Numerically. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    
    return mp.diff(lambda tau:DedekindEtaA(tau),tau,1)

# ----------------------------------------------------------------------------
def DedekindEtaA1(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        Numerically. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False
    return 1.0/mp.sqrt(3.0)*mp.jtheta(2,mp.pi/6,mp.exp(-(mp.pi/3.0)*tau))

# ----------------------------------------------------------------------------
def DedekindEtaA2(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        Numerically. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    
    return mp.exp(-mp.pi/12.0)*mp.jtheta(3,mp.pi*(mp.j*tau+1.0)/2.0,mp.exp(-3.0*mp.pi))


# ----------------------------------------------------------------------------
def DedekindEtaA4(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        Numerically. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    
    return mp.cbrt(0.5*mp.jtheta(2,0,mp.exp(-mp.pi*tau))*mp.jtheta(3,0,mp.exp(-mp.pi*tau))*mp.jtheta(4,0,mp.exp(-mp.pi*tau)))

# ----------------------------------------------------------------------------
def d2DedekindEta(tau):
    ''' Compute the 2nd derivative of the Dedekind Eta function for imaginary argument tau. 
        See: http://mathworld.wolfram.com/DedekindEtaFunction.html. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    qbar = mp.exp(-2.0*mp.pi*tau)

    res = mp.nsum(lambda n: ((-1)**n)*((6.0*n-1.0)**4)*(qbar**((6.0*n-1.0)*(6.0*n-1.0)/24.0)),[-mp.inf,mp.inf])
    return ((mp.pi/12.0)**2)*res 

def d2DedekindEtaA(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        Numerically.'''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    
    return mp.diff(lambda tau:dDedekindEtaA(tau),tau,1)


# ----------------------------------------------------------------------------
def dThetaRatio(z,q,n):
    ''' Compute the ratio of the nth derivative of the third Jacobi theta function
        with itself.'''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
        
    #return (mp.djtheta(3,z,q,n)/mp.jtheta(3,z,q))
    return (mp.jtheta(3,z,q,n)/mp.jtheta(3,z,q))

# ----------------------------------------------------------------------------
class PyFit:
    ''' A class which performs non-linear regression fits to 1D data. '''

    # Test that we can load the desired modules
    try:
        import numpy as n
        numpy_loaded = True
    except ImportError:
        numpy_loaded = False 

    try:
        import scipy.optimize.minpack as min
        scipy_loaded = True
    except ImportError:
        scipy = False 


    # ------------------------------------------------------------------------
    def __init__(self,initialFitPar,fitFunc,inX,inY,err=None,funcPar=None,dFunc=None):
        ''' We setup all local class variables, which may be overwritten
            usint initFit. '''

        # Setup all class variables
        self.params    = funcPar
        self.function  = fitFunc
        self.dfunction = dFunc
        self.x         = inX
        self.y         = inY
        self.err       = err 
        self.aFit      = initialFitPar

        self.dataLength = len(inX)
        self.fitLength  = len(self.aFit)

        # Initialize and update the weights
        self.weight = self.n.ones([self.dataLength],float)
        self.getWeights()

    # ------------------------------------------------------------------------
    def objective(self,a): 
        '''Calculate the residual deviation with the fitting function. '''

        obj = self.n.zeros([self.dataLength],float)
        print a

        for i in range(self.dataLength):
#           print self.weight[i],self.y[i],self.function(self.x[i],a,self.params)
            #obj[i] = self.weight[i]*(self.y[i] - self.function(self.x[i],a,self.params))
            obj[i] = self.weight[i]*(self.y[i] - self.function(self.x[i],a))

        return obj

    # ------------------------------------------------------------------------
    def dobjective(self,a): 
        '''Calculate the residual deviation with the fitting function. '''

        dobj = self.n.zeros([self.dataLength,len(a)],float)

        for i in range(self.dataLength):
            dobj[i,:] = -self.weight[i]*self.dfunction(self.x[i],a,self.params)

        return dobj

    # ------------------------------------------------------------------------
    def getWeights(self):
        ''' Get the appropriate weights. ''' 
        if self.err is None:
            self.weight = self.n.ones([self.dataLength],float)
        else:
            for i in range(self.dataLength):
                if (self.err[i] == 0.0):
                    self.weight[i] = 1.0
                else:
                    self.weight[i] = 1.0/self.err[i]

    # ------------------------------------------------------------------------
    def updateFit(self,initx,inity,inita,err=None,funcPar=None): 
        ''' Update the data that will be fit. '''

        self.params = funcPar
        self.x = initx
        self.y = inity
        self.err = err
        self.aFit = inita

        # Update the weights
        self.getWeights()

    # ------------------------------------------------------------------------
    def getFit(self): 
        ''' Perform the actual least-squares fit. '''
        
        # Have we supplied a Jacobian for our fit function?
        if self.dfunction is not None:
            df = self.dobjective
        else:
            df = None

        # Perform the fit
        self.aFit,covMat,infodict,mesg,ierr = self.min.leastsq(self.objective, self.aFit,\
                Dfun=df, epsfcn=1.0E-9,full_output=True) 
        print 'Object fit: ',self.aFit
        # Get chi^2
        #self.yFit = self.n.zeros([self.dataLength],float)
        #self.chi2 = 0.0
        #dof = 0
        #for i in range(self.dataLength):
        #    #self.yFit[i] = self.function(self.x[i],self.aFit,self.params)
        #    self.yFit[i] = self.function(self.x[i],self.aFit)
        #    #if self.err[i] > 1.0E-12:
        #    #    dof += 1
        #    #    self.chi2 += (self.y[i]-self.yFit[i])**2/(self.err[i]**2)
        #    #self.chi2 += (self.y[i]-self.yFit[i])**2/abs(self.yFit[i])
        #self.chi2 /= 1.0*(dof - self.fitLength - 1)

        ## this is the adjusted covariance matrix
        #s_sq = (self.objective(self.aFit)**2).sum()/(self.dataLength-self.fitLength) 
        #covMat = covMat * s_sq

        ## Find the error in the fit parameters
        #if not isList(self.aFit):
        #    self.aErr = self.n.sqrt(covMat[0,0])
        #else:
        #    self.aErr = []
        #    for i in range(self.fitLength):
        #        self.aErr.append(self.n.sqrt(covMat[i,i]))
