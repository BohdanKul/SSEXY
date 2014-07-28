#include <iostream>
#include <fstream>
#include "communicator.h"
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <time.h>

using namespace std;
//using boost::lexical_cast;

Communicator::Communicator(int _Nx, int _Ny, int _r, float _T, float _Beta, long _p, string rfName, int _Asize, bool _measTime)
{
    p = _p;
    
    //Define temperature in one of two ways
    if (_T == -1){
        _T = 1.0/(1.0*_Beta);
         dataName = boost::str(boost::format("%02d-%03d-%03d-%s%06.3f") %_r %_Nx %_Ny %"b" %_Beta);
    }
    else dataName = boost::str(boost::format("%02d-%03d-%03d-%s%06.3f") %_r %_Nx %_Ny %"t" %_T);
    
    if  (not(_Asize<0)){
        dataName += boost::str(boost::format("-%04d") %_Asize);
    }

    types  = vector<string> {"state","estimator"};
    outDir = "OUTPUT"; 
    
    //Generate or fetch the id
    if  (rfName=="") 
        GenerateId();
    else{
        string fileName = string(find( rfName.rbegin(), rfName.rend(), '/').base(), rfName.end());
        id = atol(fileName.substr(25,9).c_str());
        }

    string  fileName;
    for (vector<string>::iterator type=types.begin(); type!=types.end(); type++){
        fileName = boost::str(boost::format("%s/%s-%s-%09d.dat") %outDir %*type %dataName %id);
        //State files are read-only at first
        if ((rfName!="") and (*type=="state")){
            cout << fileName;    
            mFStreams[*type] = new fstream(fileName,ios_base::in);
        }
        //Other files are write by appending
        else
            mFStreams[*type] = new fstream(fileName,ios_base::out|ios_base::app);
        if (!*mFStreams[*type]){
           cerr << "Unable to process file: " << fileName << endl;
           exit(EXIT_FAILURE); 
        }
        //cout << boost::str(boost::format("%s-%|15t|-%s-%|40t|-%010.3f") %*type %dataName %10.1111) << endl; //%lexical_cast<string>(id));                  
    }
    
}

void Communicator::GenerateId()
{
    time_t seconds = long(time(NULL) - 39*365*24*60*60);
    id = long(seconds + p);
    string fName;
    fName = boost::str(boost::format("OUTPUT/estimator-%s-%09d.dat") % dataName %id);
    boost::filesystem::path logPath(fName);

    while(boost::filesystem::exists(logPath)) {
        id += 1;
        fName = boost::str(boost::format("OUTPUT/estimator-%s-%09d.dat") % dataName %id);
        logPath = boost::filesystem::path(fName);

    }

}

void Communicator::reset(string _fileName)
{
    mFStreams[_fileName]->close();
    string fileName = boost::str(boost::format("%s/%s-%s-%09d.dat") %outDir %_fileName %dataName %id);
    mFStreams[_fileName]->open(fileName, fstream::out|fstream::trunc);
}

fstream* Communicator::stream(string _fileName)
{
    return mFStreams[_fileName];

}


