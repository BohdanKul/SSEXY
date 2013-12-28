#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <fstream>
#include <vector>
#include <unordered_map>

using namespace std;

class Communicator
{
    public:
        Communicator(int Nx, int Ny, float T, long process);
        fstream* stream(string _fileName); 
        string dataName;
        string outDir;
        vector <string>   types;
        vector <fstream> files;   

    private:
        long id;
        long p;
        void GenerateId();
        unordered_map <string,fstream*> mFStreams;
}; 
#endif
