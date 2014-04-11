#ifndef LATTICE_H
#define LATTICE_H

#include <fstream>
#include <vector>

using namespace std;

class LATTICE
{
    public:
        LATTICE(string _Lname, const char* _lfile);
        vector<long>* getLattice(){return &lattice;}

        int           getSize(){return size;}        
        bool          isDefined(){return _isDefined;}
    private:
        int Generate(int _size);
        int Load(const char* _lfile);
    
        string       Lname;    //Lattice name
        vector<long> lattice;  //Lattice data structure 
        int          size;   //Lattice size 
        bool         _isDefined;
}; 
#endif
