#!/usr/bin/python 
#
# gensubmit.py
# Adrian Del Maestro
# 06.03.2009
#
# Generate a torque submit script for the pimc code which creates
# a pbs file for various sets of parameters.

import os,sys,clusters
from argparse import ArgumentParser


# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

    # setup the command line parser options 
    parser = ArgumentParser(description="Build submission scripts for various clusters") 
    parser.add_argument("file", help='configuration file')
    parser.add_argument('-b','--batch', type=int, help='Turn on batch mode. Sets the number of RG seeds.')
    parser.add_argument("cluster", metavar="cluster", choices=['clumeq','westgrid','sharcnet','scinet','bluemoon'],\
            help="target cluster: [westgrid,sharcnet,scinet,clumeq,bluemoon]") 
    parser.add_argument("-r",type=str, dest='run', default="",help="optional JobId number that will be added to the scripts name")
    # parse the command line options
    args = parser.parse_args() 
    inFileName = args.file

    if (not args.cluster):
        parser.error("need to specify a cluster")

    # We open up the input file, and read in all lines.
    inFile = open(inFileName,'r')
    inLines = inFile.readlines();

    # The first line of the file contains all static pimc options
    staticPIMCOps = inLines[0].rstrip('\n')

    # The next lines contains the short-name of the pimc option
    # and a list of values.  
    optionValue = {}
    numOptions = {}
    for line in inLines[1:]:
        option = line.split()
        flag = option.pop(0)
        # Now we determine if we have a range or if we have actually included the values
        if option[0].find(':') != -1:

            # determine if we have floats or ints
            if option[0].find('.') != -1:
                type_converter = lambda x: float(x)
            else:
                type_converter = lambda x: int(x)

            dataRange = option[0].split(':')
            print dataRange
            # Parse the range string
            dataMin = type_converter(dataRange[0])
            dataMax = type_converter(dataRange[1])
            dataStep = type_converter(dataRange[2])

            # construct the option list
            vals = [dataMin]
            while vals[-1] < dataMax:
                vals.append(vals[-1]+dataStep)

            # assign the option list
            if flag in optionValue.keys():
	       optionValue[flag] += vals
	       numOptions[flag]  +=len(vals)
            else:
		optionValue[flag] = vals
                numOptions[flag]  = len(vals)
       	else:
            # store the typed out values
            if flag in optionValue.keys():
	       optionValue[flag] +=option
	       numOptions[flag]  +=len(option)
       	    else:
            	optionValue[flag] = option
            	numOptions[flag]  = len(option)

    # We make sure all option strings have the same length 
    if numOptions:
        num = numOptions[flag]
        for (flag,n) in numOptions.iteritems():
           print (flag,n)        
	   if n != num:
                print flag, ' ', n
		print 'Not all parameters have the same number of values!'
                sys.exit()
                
    commandLines = []
    print optionValue
    index = 0
    for i in range(0,num):
        commandLine = './ssexy.e '
        for flag,val in optionValue.iteritems():
            commandLine += '-%s %s ' % (flag,val[i])
        commandLine += staticPIMCOps
        commandLines.append(commandLine)

    ncommandLines = []
    if  args.batch!=None:
        for j in range(args.batch):
            for commandLine in commandLines:
                ncommandLines.append(commandLine + ' -p %d ' %j)
        commandLines=ncommandLines
    
    if args.cluster == 'westgrid':
        clusters.westgrid(commandLines,args.run)

    if args.cluster == 'sharcnet':
        clusters.sharcnet(commandLines,args.run)

    if args.cluster == 'scinet':
        clusters.scinet(commandLines,args.run)

    if args.cluster == 'clumeq':
        clusters.clumeq(commandLines,args.run)

    if args.cluster == 'bluemoon':
        clusters.bluemoon(commandLines,args.run)

# ----------------------------------------------------------------------
if __name__ == "__main__": 
    main()
