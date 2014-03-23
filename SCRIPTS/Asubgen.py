from numpy import arange
from numpy import ones_like
from numpy import array 
import argparse

def main(): 
    parser = argparse.ArgumentParser(description='Plot Raw MC Equilibration Data for Scalar Estimators.')
    parser.add_argument('-L', help='Length', type=int)

    args = parser.parse_args()
    
    if  not args.L: 
        parser.error("Specify the dimension")
    else:
         maxL = args.L 
    midL = maxL % 2
    maxN = maxL**2 
    midN = maxL**2 // 2 
    remainder = (midN) % 8
    Lrange = [0]
    if   remainder ==  7: Lrange = [0, 7]
    elif remainder ==  6: Lrange = [0, 6]
    elif remainder ==  5: Lrange = [0, 6, 12]
    elif remainder ==  4: Lrange = [0, 6, 11]
    elif remainder ==  3: Lrange = [0, 10, 18]
    elif remainder ==  2: Lrange = [0, 10]
    elif remainder ==  1: Lrange = [0, 9]
    aLrange = array(Lrange)
    Hrange = ones_like(aLrange)*maxN - aLrange
    Hrange = list(Hrange)[::-1]
    Lrange = Lrange + range(midN,max(Lrange),-8)[::-1]
    Hrange = range(midN,min(Hrange),8)[1:] + Hrange

    Arange = Lrange + Hrange
    Astep  = list(array(Arange)[1:] - array(Arange)[0:-1]) + [0]
    print 'a', " ".join(map(str,Arange))
    print 't', " ".join(map(str,Astep))
    print 'p', "1:"+str(len(Arange))+":1"


if __name__ == "__main__": 
    main()
