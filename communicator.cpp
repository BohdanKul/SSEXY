#include <iostream>
#include <fstream>
#include "communicator.h"
#include <boost/format.hpp>

using namespace std;
//using boost::lexical_cast;

Communicator::Communicator(int Nx, int Ny, float T)
{
    types  = vector<string> {"spin","link","loop","operator","estimator","bond","vertex"};
    outDir = "OUTPUT"; 
    GenerateId();
    dataName = boost::str(boost::format("%03d-%03d-%06.3f-%09d") %Nx %Ny %T %id);
    string  fileName;
    for (vector<string>::iterator type=types.begin(); type!=types.end(); type++){
        fileName = boost::str(boost::format("%s/%s-%s.dat") %outDir %*type %dataName);
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
    id = 301021023;    
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

