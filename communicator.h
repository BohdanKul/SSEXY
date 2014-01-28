#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <fstream>
#include <vector>
#include <unordered_map>

using namespace std;

class Communicator
{
    public:
        Communicator(int Nx, int Ny, int _r, float _T, float _Beta, long process, string rfName, int _maxSpin);
        fstream* stream(string _fileName); 
        void     reset(string _fileName);
        long     getId(){return id;};

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
