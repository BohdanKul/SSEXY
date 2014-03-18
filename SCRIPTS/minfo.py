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
    for i in range(size):
        elem = l.pop()
        if  elem in l:
            print "Not unique temperature ", elem

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
    args = parser.parse_args() 

    rcParams.update(mplrc.aps['params'])
    colors = ["#66CAAE", "#CF6BDD", "#E27844", "#7ACF57", "#92A1D6", "#E17597", "#C1B546",'b']
    figure(1,(7,6))
    connect('key_press_event',kevent.press)
    print args.energy
    if args.renyi:    
        if not(args.energy) and not(args.partition):
            ax  = subplot(211)
            ax1 = subplot(212)
        else:
            ax  = subplot(311)
            ax1 = subplot(312)
            ax2 = subplot(313)
            ax2.set_xlabel(r'$\mathrm{%s[K^{-1}]}$' %parMap['b'])
            ax2.set_ylabel(r'$\mathrm{E/N [K^{-1}]}$')
        ax1.set_ylabel(r'$\mathrm{S_2/L}$')
        ax1.set_xlabel(r'$\mathrm{%s}[K^{-1}]$' %parMap['b'])
    else:
        if not(args.energy) and not(args.partition):
            ax  = subplot(111)
        else:
            ax  = subplot(211)
            ax2 = subplot(212)
            ax2.set_xlabel(r'$\mathrm{%s[K^{-1}]}$' %parMap['b'])
            ax2.set_ylabel(r'$\mathrm{E/N [K^{-1}]}$')
        
    ax.set_ylabel(r'$\mathrm{I/L}$')


    order = ['Lx','Ly','T','b']
    
    Lxs = {}
    offset = 1
    for fileName in args.fileNames:
        scalarhelp = ssexyhelp.ScalarReduce(fileName)
        parmap = scalarhelp.getParmap()
        Lx = parmap['Lx']
        if  not Lx in Lxs:
            Lxs[Lx] = {}

        geom = GetDicKey(parmap)[0]
        r    = parmap['r']
        Lxs[Lx][(geom,r)] = scalarhelp


    data = {}
    i = 0
    j = 0
    equalA = False
    shift = 0
    for Lx, dic in Lxs.items():
        for geom, rscalar in dic.items():
            rscalar.loadData()
            rb, Beta = rscalar.getrParams()
            if (1.0/Beta[0]) > 20:
                print "Inserting value at infinity"
                Beta = np.insert(Beta,0,0)
            
            Beta = Beta[shift:]
            Norm = int(Lx)**2 
            print geom
            CheckUniq(Beta)
            if j>0: CheckEqual(Beta2,Beta)
            Beta2 = Beta
           
            if   geom == ('half',2):
                 equalA = True
                 Eh2, dEh2    = rscalar.getAverages('ET')
                 if (1.0/Beta[1]) > 20:
                    Eh2  = np.insert(Eh2,0,0)
                    dEh2 = np.insert(dEh2,0,0)
                 Eh2  = Eh2[shift:]
                 dEh2 = dEh2[shift:]
                 Eh2  = Eh2*Norm
                 dEh2 = dEh2*Norm
                 Zh2,BetaZ = GetZ(Beta,Eh2,dEh2,offset)
                 if args.energy:
                     if i==0: ax2.errorbar(Beta, Eh2, dEh2, linestyle = 'None', color = colors[i], marker='o', label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.errorbar(Beta, Eh2, dEh2, linestyle = 'None', color = colors[i], marker='o')
                 if args.partition:
                     if i==0: ax2.plot(BetaZ, Zh2, linestyle = '-', color = colors[i], label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.plot(BetaZ, Zh2, linestyle = '-', color = colors[i])
            elif geom == ('A',2):
                 EA2, dEA2    = rscalar.getAverages('ET')
                 if (1.0/Beta[1]) > 20:
                    EA2  = np.insert(EA2,0,0)
                    dEA2 = np.insert(dEA2,0,0)
                 EA2  = EA2*Norm
                 dEA2 = dEA2*Norm
                 ZA2,BetaZ = GetZ(Beta,EA2,dEA2,offset)
                 if args.energy:
                     if i==0: ax2.errorbar(Beta, EA2, dEA2, linestyle = 'None', color = colors[i], marker='o', label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.errorbar(Beta, EA2, dEA2, linestyle = 'None', color = colors[i], marker='o')
                 if args.partition:
                     if i==0: ax2.plot(BetaZ, ZA2, linestyle = '-', color = colors[i], label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.plot(BetaZ, ZA2, linestyle = '-', color = colors[i])
            elif geom == ('B',2):
                 EB2, dEB2    = rscalar.getAverages('ET')
                 if (1.0/Beta[1]) > 20:
                    EB2  = np.insert(EB2,0,0)
                    dEB2 = np.insert(dEB2,0,0)
                 EB2  = EB2*Norm
                 dEB2 = dEB2*Norm
                 ZB2,BetaZ = GetZ(Beta,EB2,dEB2,offset)
                 if args.energy:
                     if i==0: ax2.errorbar(Beta, EB2, dEB2, linestyle = 'None', color = colors[i], marker='o', label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.errorbar(Beta, EB2, dEB2, linestyle = 'None', color = colors[i], marker='o')
                 if args.partition:
                     if i==0: ax2.plot(BetaZ, ZB2, linestyle = '-', color = colors[i], label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.plot(BetaZ, ZB2, linestyle = '-', color = colors[i])
            elif geom == ('full',2):
                 Ef2, dEf2    = rscalar.getAverages('ET')
                 if (1.0/Beta[1]) > 20:
                    Ef2  = np.insert(Ef2,0,0)
                    dEf2 = np.insert(dEf2,0,0)
                 Ef2  = Ef2[shift:]
                 dEf2 = dEf2[shift:]
                 Ef2  = Ef2*Norm
                 dEf2 = dEf2*Norm
                 Zf2,BetaZ = GetZ(Beta,Ef2,dEf2, offset)
                 if args.energy:
                     if i==0: ax2.errorbar(Beta, Ef2, dEf2, ls='None', color = colors[i], marker='s', label=r'$\mathrm{r=2;\, A=full}$')
                     else:    ax2.errorbar(Beta, Ef2, dEf2, ls='None', color = colors[i], marker='s')
                 if args.partition:
                     if i==0: ax2.plot(BetaZ, Zf2, linestyle = '-', color = colors[i], label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.plot(BetaZ, Zf2, linestyle = '-', color = colors[i])
            
            elif geom == ('full',1):
                 Ef1, dEf1    = rscalar.getAverages('ET')
                 if (1.0/Beta[1]) > 20:
                    Ef1  = np.insert(Ef1,0,0)
                    dEf1 = np.insert(dEf1,0,0)
                 Ef1  = Ef1[shift:]
                 dEf1 = dEf1[shift:]
                 Ef1  = Ef1*Norm
                 dEf1 = dEf1*Norm
                 Zf1,BetaZ = GetZ(Beta,Ef1,dEf1, offset)
                 if args.energy:
                     if i==0: ax2.errorbar(np.array(Beta), Ef1, dEf1,ls='None', color = colors[i], marker='>', label=r'$\mathrm{r=1;\, A=full}$')
                     else:    ax2.errorbar(np.array(Beta), Ef1, dEf1,ls='None', color = colors[i], marker='>')
                 if args.partition:
                     if i==0: ax2.plot(BetaZ, Zf1, linestyle = '-', color = colors[i], label=r'$\mathrm{r=2;\, A=half}$')
                     else:    ax2.plot(BetaZ, Zf1, linestyle = '-', color = colors[i])
            j += 1
        T, MI = loadtxt("MI_ED.dat", unpack=True)
        ax.plot(1.0/np.array(T), MI/float(Lx),\
                marker='',color=colors[-1],\
                label=r'$\mathrm{ED}$')
        dic16x01 = {0.5: 0.0013439740551892676, 1.0: 0.004540510233360606, 2.0: 0.013825090870731693, 3.0: 0.023659939402554397, 4.0: 0.03329452300571986, 5.0: 0.041639084124480946, 0.1: 5.3647025572756135e-05, 7.0: 0.053496549803275645, 8.0: 0.05783982372091184, 9.0: 0.0616724129728039, 10.0: 0.06505550659201656, 6.0: 0.04789194046993728, 0.01: -0.0, 1.5: 0.008898557058858003} 
        #dic04x04 = {0.5: 0.022500418440747794, 1.0: 0.08979659261091766, 2.0: 0.19522255433112523, 4.0: 0.26849339636811287, 0.1: 0.0009457518941242693, 1.5: 0.1563080540879589, 7.0: 0.33251641242350705, 8.0: 0.35013882821562525, 9.0: 0.3664919686257604, 10.0: 0.3826062726961008, 2.75: 0.2313212484936219, 2.25: 0.20777695945606411, 6.0: 0.3135947068146292, 5.0: 0.29299056784373767, 1.75: 0.17724769701174903, 0.01: -0.0, 1.25: 0.12785877158962244} 
        dic04x04  =  {0.5: 0.07020309184159479, 0.75: 0.17339358634204705, 2.75: 0.331474482134048, 4.0: 0.3549492763332615, 1.1: 0.30820887423118226, 0.1: 0.002400609577131174, 1.0: 0.2818067715190119, 8.0: 0.4229487601228128, 1.75: 0.3266143562914207, 0.9: 0.24371549496221712, 2.25: 0.32212659675366573, 6.0: 0.39000909529160377, 2.5: 0.32758479927105966, 1.5: 0.3297757245965327, 0.01: 0.00016144356014557992, 1.25: 0.330461291461975}
        dic04x04e =  {0.5: 0.0003831777430893546, 0.75: 0.0004675204681600513, 2.75: 0.0008917989478065976, 4.0: 0.0008956688447896057, 1.1: 0.0008952625741067421, 0.1: 0.000302110811727677, 1.0: 0.0005032714435446223, 8.0: 0.0009852341788363702, 1.75: 0.000636725768302311, 0.9: 0.0008710658144792192, 2.25: 0.0006451020052804455, 6.0: 0.0009645559487371858, 2.5: 0.0006282633150395536, 1.5: 0.0008893248430503498, 0.01: 0.0003161381258478166, 1.25: 0.0006137569649840289}
        dic04x04A  =  {0.5: 0.10151677825047312, 1.0: 0.3629985489039014, 1.25: 0.41212133372654125, 4.0: 0.40733595847152876, 0.1: 0.0033365649569256917, 1.1: 0.39561309309807613, 8.0: 0.47970763147821227, 1.5: 0.4065193074056275, 0.9: 0.3185515971303039, 2.25: 0.3776005950606687, 6.0: 0.44819954052931366, 1.75: 0.388433356790134, 2.5: 0.37995684759477166, 2.75: 0.38282874713989407} 
        dic04x04Ae =  {0.5: 0.0007450409825502224, 1.0: 0.0009481296290678376, 1.25: 0.0010187346802027124, 4.0: 0.0010978116222615142, 0.1: 0.0006106006697038598, 1.1: 0.0010203347930939602, 8.0: 0.0011901492541913123, 1.5: 0.0010339790609518711, 0.9: 0.0009323954510706496, 2.25: 0.00104793668898419, 6.0: 0.001148244412638139, 1.75: 0.0010281105861974981, 2.5: 0.001085870766010387, 2.75: 0.0010694040757170682}
        dic  = dic04x04A
        dice = dic04x04Ae
        ax.errorbar(dic.keys(),np.array(dic.values()),np.array(dice.values()),\
                    ls='', marker='s',color=colors[-2],\
                    label =r'$\mathrm{ratio}$')

        if  equalA:
            I  = Zf2 + 2*Zf1 - 2*Zh2
            S2A  = -Zh2 + 2*Zf1 + np.log(2)*(int(Lx)**2)*0.5
            S2B  = -Zh2 + 2*Zf1 + np.log(2)*(int(Lx)**2)*0.5
            S2AB = -Zf2 + 2*Zf1 + np.log(2)*(int(Lx)**2)
        else:
            I = Zf2 + 2*Zf1 - ZA2 - ZB2
            S2A  = -ZA2 + 2*Zf1  
            S2B  = -ZB2 + 2*Zf1  
            S2AB = -Zf2 + 2*Zf1  
            #for j in range(50):
            #    ZA2,Beta1 = GetZ(Beta,EA2,dEA2,offset)
            #    ZB2,Beta1 = GetZ(Beta,EB2,dEB2,offset)
            #    Zf2,Beta1 = GetZ(Beta,Ef2,dEf2, offset)
            #    Zf1,Beta1 = GetZ(Beta,Ef1,dEf1, offset)
        
        if args.renyi:
            FancyErrorbar(ax1,np.array(BetaZ),(S2A*1.0)/float(Lx),\
                  colors[i],((r'$\mathrm{Lx=%2.0d}$' %Lx) +'\n'+ r'$\mathrm{S_2A}$'))
            if not equalA:
                FancyErrorbar(ax1,np.array(BetaZ),(S2B*1.0)/float(Lx),\
                      colors[i],r'$\mathrm{S_{2B}}$')
            FancyErrorbar(ax1,np.array(BetaZ),(S2AB*1.0)/float(Lx),\
                  colors[i],r'$\mathrm{S_{2A \bigcup B}}$', lines='--')
            dicS2F04x04  =  {0.5: 2.510650769895426, 0.75: 2.1498839339301865, 2.75: 0.38014521192140444, 4.0: 0.29546525280083485, 1.1: 1.5044519736841986, 0.1: 2.7622900523274256, 1.0: 1.686407094237605, 8.0: 0.1980272884529148, 1.75: 0.7072501116689618, 0.9: 1.875704705740789, 2.25: 0.47846276847212843, 6.0: 0.23846972985730638, 2.5: 0.42244783542458075, 1.5: 0.9288488298974135, 0.01: 2.772333451026834, 1.25: 1.2571854190686764} 
            dicS2F04x04e =  {0.5: 0.0003831777430893546, 0.75: 0.0004675204681600513, 2.75: 0.0008917989478065976, 4.0: 0.0008956688447896057, 1.1: 0.0008952625741067421, 0.1: 0.00030211081172767696, 1.0: 0.0005032714435446223, 8.0: 0.0009852341788363702, 1.75: 0.000636725768302311, 0.9: 0.0008710658144792192, 2.25: 0.0006451020052804455, 6.0: 0.0009645559487371858, 2.5: 0.0006282633150395537, 1.5: 0.0008893248430503497, 0.01: 0.0003161381258478166, 1.25: 0.0006137569649840289}
            dicS2A04x04  =  {0.5: 1.2904269308685103, 0.75: 1.1616387601361167, 2.75: 0.3558098470277262, 4.0: 0.3252072645670482, 1.1: 0.9063304239576904, 0.1: 1.3823453309522784, 1.0: 0.9841069328783084, 8.0: 0.3104880242878638, 1.75: 0.5169322339801913, 0.9: 1.059710100351503, 2.25: 0.4002946826128971, 6.0: 0.3142394125744551, 2.5: 0.3750163173478202, 1.5: 0.6293122772469731, 0.01: 1.3862474472934898, 1.25: 0.7938233552653258} 
            dicS2A04x04e =  {0.5: 0.0002589116626442326, 0.75: 0.0003201733812410127, 2.75: 0.0006010387765504837, 4.0: 0.0006303926164098645, 1.1: 0.0006061105056466294, 0.1: 0.00020855636580185746, 1.0: 0.0003406438541255847, 8.0: 0.0006552075585454876, 1.75: 0.00046048335584552896, 0.9: 0.0005928122941938723, 2.25: 0.0004628809993232214, 6.0: 0.0006753574613881269, 2.5: 0.00042803429740009575, 1.5: 0.0006460616230140156, 0.01: 0.0002189592780460406, 1.25: 0.00042978216760510745}
            ax1.errorbar(dicS2F04x04.keys(),np.array(dicS2F04x04.values()),np.array(dicS2F04x04e.values()),\
                        ls='', marker='s',color=colors[-2],\
                        label =r'$\mathrm{S_{2F}}$')
            ax1.errorbar(dicS2A04x04.keys(),np.array(dicS2A04x04.values()),np.array(dicS2A04x04e.values()),\
                        ls='', marker='s',color=colors[-3],\
                        label =r'$\mathrm{S_{2A}}$')
        FancyErrorbar(ax,np.array(BetaZ),(I*1.0)/float(Lx),\
              colors[i],r'$\mathrm{Integration}$')#Lx=%2.0d}$' %Lx)

        i += 1

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

