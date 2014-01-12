'''savefig.py

A set of subroutines to store pylab plots in a zip file. 
Also, it includes functions to create and update a text-file
containing code that can be used to display those figures on 
MoinMoin. This txt-file is automatically saved in the zip file 
and gets updated everytime the zip-archive is changed. 
'''

import os
import zipfile

# -----------------------------------------------------------------------------
# MoinMoin related functions
# -----------------------------------------------------------------------------
def createMoinMoin(fnameMoinMoin):
    '''create a new MoinMoin txt-file'''

    File = open(fnameMoinMoin,'w')
    File.write(' * Reduce-one\n')
    File.write(' * Pimcplot\n')
    File.close()
   
    return

# -----------------------------------------------------------------------------
def appendCommandMoinMoin(filename,command,TYPE):
    '''appends a new command to the MoinMoin filename 
       txt-file while respecting its structure'''
    File = open(filename,'r')
    dataLines = File.readlines()
    index = dataLines.index(' * %s\n' %TYPE)
    dataLines.insert(index+1,command)
    File.close()
    os.remove(filename)
    File = open(filename,'w')
    File.writelines(dataLines)
    File.close()
   
    return     
   
# -----------------------------------------------------------------------------
def getCommandMoinMoin(estimator,filename,TYPE):
    '''Create MoinMoin plot displsy command string '''
 
    Command = '  * {{{%s}}}\n\
    {{attachment:%s|%s|width="1000"}}\n' % (estimator,filename,TYPE)       
   
    return Command 
  
# -----------------------------------------------------------------------------
#  Zip archive related functions
# -----------------------------------------------------------------------------
def removeZipedFile(ZipName,filename):
    '''remove filename from ZipName archive'''  
    zf = zipfile.ZipFile(ZipName, 'r') 
    tzf = zipfile.ZipFile('temp_'+ZipName, 'w')
    for zipedFile in zf.namelist():
        if zipedFile != filename:
           data = zf.read(zipedFile)
           tzf.writestr(zipedFile, data)
    zf.close()
    tzf.close()
    os.rename('temp_'+ZipName,ZipName) 

# -----------------------------------------------------------------------------
def saveToZip(ZipName,FileName,run,estimator,TYPE,period=''):
    '''save FileName to ZipName archieve. 
       Same instances of FileName get overwritten. 
       ${run}_MoinMoin.txt contains code that can be used to 
       display the figures, once unzipped, on MoinMoin'''

    fnameMoinMoin = 'run-%s_MoinMoin.txt' % run
    Filexist = False 
    if not(os.path.exists(ZipName)): 
       zf = zipfile.ZipFile(ZipName, mode='w')
       createMoinMoin(fnameMoinMoin)
    else:
        zf = zipfile.ZipFile(ZipName, mode='r')   
        zf.extract(fnameMoinMoin)
        removeZipedFile(ZipName,fnameMoinMoin) 
        if FileName in zf.namelist():
            Filexist = True  
            removeZipedFile(ZipName,FileName)  
        zf = zipfile.ZipFile(ZipName, mode='a')

    if not(Filexist):
       if TYPE == 'Pimcplot':
          Command = getCommandMoinMoin(estimator,FileName,'Period = %s' % period)
          appendCommandMoinMoin(fnameMoinMoin,Command,TYPE)
          #Command = getCommandMoinMoin(estimator,FileName,'Cumulative average')
          #appendCommandMoinMoin(fnameMoinMoin,Command,TYPE)
       if TYPE == 'Reduce-one':
          Command = getCommandMoinMoin(estimator,FileName,'Reduce-one')
          appendCommandMoinMoin(fnameMoinMoin,Command,TYPE)
    zf.write(fnameMoinMoin)
    zf.write(FileName)
    zf.close()  
    os.remove(fnameMoinMoin)
    os.remove(FileName)

# -----------------------------------------------------------------------------
def getPmcName(args,ParamValue):
    '''produce the names of a png file'''

    if (args.ParamName != None) and (len(ParamValue) != 0) : 
        minParamValue = sorted(ParamValue)[0]
        maxParamValue = sorted(ParamValue)[-1]
        FileName = "run-%s_p=%s_%s_%s=%s:%s.png" %(args.run,args.period,args.estimator.replace('/','over'),args.ParamName,minParamValue,maxParamValue)
    else:
        FileName = "run-%s_p=%s_%s.png" %(args.run,args.period,args.estimator.replace('/','over'))

    return FileName

# -----------------------------------------------------------------------------
def getReduceName(args,plot,ParamValue):
    '''produce the names of a png file'''

    minParamValue = sorted(ParamValue)[0]
    maxParamValue = sorted(ParamValue)[-1]
    FileName = "run-%s_reduced%s_%s_%s=%s:%s.png" %(args.run,args.reduce,plot,args.reduce,minParamValue,maxParamValue)


    return FileName

# -----------------------------------------------------------------------------
def getZipName(run):
    '''produce the name of a zip file'''

    ZipName = "run-%s.zip" %(run) 

    return ZipName
