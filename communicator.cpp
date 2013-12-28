#include <iostream>
#include <fstream>
#include "communicator.h"
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <time.h>

using namespace std;
//using boost::lexical_cast;

Communicator::Communicator(int _Nx, int _Ny, float _T, long _p)
{
    p = _p;
    types  = vector<string> {"spin","link","loop","operator","estimator","bond","vertex"};
    outDir = "OUTPUT"; 
    GenerateId();
    dataName = boost::str(boost::format("%03d-%03d-%06.3f") %_Nx %_Ny %_T);
    string  fileName;
    for (vector<string>::iterator type=types.begin(); type!=types.end(); type++){
        fileName = boost::str(boost::format("%s/%s-%s-%09d.dat") %outDir %*type %dataName %id);
        mFStreams[*type] = new fstream(fileName,ios_base::out);
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

fstream* Communicator::stream(string _fileName)
{
    return mFStreams[_fileName];


}

//int main()
//{
//    Communicator communicator(4,4,1.01);
//    cout << communicator.types[0] << endl;
//    *communicator.stream("spin") << "yo" << endl;
//    return 0;
//}

