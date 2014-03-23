from uncertainties   import unumpy

def MinMax(data):
    values = unumpy.nominal_values(data)
    errors = unumpy.std_devs(data)
    mins = values - errors
    maxs = values + errors
    return values,mins,maxs
        

def FancyErrorbar(ax,x,data,col,lab,lines='-'):
    values, mins, maxs = MinMax(data)    
    ax.fill_between(x, mins,maxs, facecolor = col, alpha =0.25, edgecolor='None')
    ax.plot(x, values+0.0, marker='',color=col,ls=lines,ms=6,mec='None', label=lab)

