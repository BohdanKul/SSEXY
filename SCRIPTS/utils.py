from uncertainties import unumpy
from simpserror    import simps_error
from scipy         import interpolate
from scipy         import integrate
import numpy as np

#==============================================================================
#==============================================================================
def MinMax(data):
    values = unumpy.nominal_values(data)
    errors = unumpy.std_devs(data)
    mins = values - errors
    maxs = values + errors
    return values,mins,maxs
        

#==============================================================================
#==============================================================================
def FancyErrorbar(ax,x,data,**kwargs):
    values, mins, maxs = MinMax(data)    
    ax.fill_between(x, mins,maxs, facecolor = kwargs['color'], alpha =0.25, edgecolor='None')
    ax.plot(x, values+0.0, marker='',ls='-',ms=6,mec='None', **kwargs)


#==============================================================================
#==============================================================================
def Bootstrap(func,x,y,dy,ns):
    ''' Generates Gaussian noise on y normalized by dy. The newly generated data 
        is analyzed with func() ns-times which is required to produce an array 
        of the same shape as x. The mean and std of accumulated func() output 
        are returned in the end '''

    datan = len(x)
    accum = np.zeros((ns,datan))
    for n in range(ns):
        noise = np.random.normal(0,1.0,datan)*dy # Generate noise
        newy  = y + noise                        # Update the y variable
        accum[n,:] = func(x,newy,dy)             # Analiyze the new data-set
    
    mean = np.average(accum,0)
    std  = np.std(accum,0)

    return mean, std



#==============================================================================
#==============================================================================
def SplineIntegrate(xs,ys,dys):
    ''' Integrates a spline fit to the data at 
        all x's in the interval [x[0],x[-1]] '''

    ospline = interpolate.UnivariateSpline(xs,ys, w=1./dys,k=3)
    ispline = ospline.antiderivative()
    C       = ispline(xs[0])

    return ispline(xs)-C

#==============================================================================
#==============================================================================
def SimpsonIntegrate(xs,ys,dys):        
    ''' Use Simpson rules to integrate a numberical data set''' 
    integral    = np.zeros(len(xs))
    integral[0] = 0
    for i in range(1,len(xs)):
        integral[i]  = integrate.simps(ys[:i],xs[:i], even='first')
    return integral

#==============================================================================
#==============================================================================
def SimpsonIntegrateError(xs,ys,dys):        
    ''' Use Simpson rules to integrate a numerical data set
        and estimate its error via uncertainty propagation''' 

    integral  = np.zeros(len(xs))
    error     = np.zeros(len(xs))
    integral[0] = 0
    for i in range(1,len(xs)):
        integral[i]  = integrate.simps(ys[:i],xs[:i], even='first')
        if  (i>2):  error[i] = simps_error(dys[:i], xs[:i], even='first')
    return integral,error


#==============================================================================
#==============================================================================
def SplineGenerate(x,y,dy,nx):
    ''' Generates a new dataset based on spline fits to the original y
        and dy data sets. '''
        
    # Fit a smooth function to y and an "exact" one to dy 
    spline  = interpolate.UnivariateSpline(x, y, w=1./dy,k=3)
    espline = interpolate.UnivariateSpline(x,dy,         k=3, s = 0.1)

    return spline(nx), espline(nx)

