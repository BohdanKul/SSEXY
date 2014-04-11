#include "lattice.h"
#include <string.h>
#include <iostream>
#include <sstream>

using namespace std;

LATTICE::LATTICE(string _Lname, const char* _cline)
{
    Lname    = _Lname;  //LATTICE name
    lattice  = {};      //LATTICE data structure 
    size     = -1;

    if  (strlen(_cline)==0){
        cout << "Taking " << _Lname << " region to be empty" << endl;
        _isDefined = false;
        return;
    }
        
    char *end;
    size = strtol(_cline, &end, 10);

    if  (size==0) {
        _isDefined = true;
        return;
    }

    if  (!*end){ 
        cout << "Automatic generation of region  " << Lname << " containing " << size << " spins" << endl;
        _isDefined = true;
        Generate(size);
    }
    else{
        cout << "Loading " << Lname << " region from: " << _cline << endl;        
        _isDefined = true;
        size = Load(_cline);
    }
}

int LATTICE::Generate(int _size)
{
    for (int i=0; i!=_size; i++)
        lattice.push_back(i); 
}

int LATTICE::Load(const char* _lfile)
{
    ifstream RAfile (_lfile);
    string line;
    string lline;
    int _size = 0;
    if  (RAfile.is_open()){
        while (getline(RAfile,line))
              lline=line;
        RAfile.close();
        istringstream sline(lline.erase(0,1));
        int n;
        while (sline >> n){
            lattice.push_back(n);
            _size += 1;
        }
    }
    else{
        cout << "Unable to process region" << Lname << " file" << endl;
        exit(0);
    }    
    return  _size;
}
