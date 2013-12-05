#ifndef REPLICA_H
#define REPLICA_H

#include <vector>
//#include "communicator.h"
#include "communicator.cpp"
#include "randombase.cpp"

using namespace std;


//Workhorse
class Replica: public RandomBase{
private:

    //Physical parameters
    unsigned short Nx,Ny;                //Lattice dimensions
    long   NBonds;               //Number of bonds
    long N;                      //Number of spin sites
    float T;                             //Temperature
    float Beta;                          //Inverse temperature
    vector <vector<long>> sites; //Sites connected by bonds  
    vector <long> spins;                  //Spins array
    vector <long> ap;                     //Propagated spins state

    //Algorithmic variables
    long ESteps;                 //Number of equilibration steps    
    long estep;
    long M;                     //Upper operator list size
    long n;                     //Current operator list size
    vector<long> sm;            //Operator vector list
    vector<long> first;         //First and last leg coordinates
    vector<long> last;          //for a particular spin site
    vector<long> links;         //Linked vertex list
    vector<long> vtx;           //Vertex type

    bool Debug;

    //File management
    Communicator communicator;
 
   
    //Methods
    float BondDiagonalEnergy(long b);
    void  LatticeGeometry();
    long  VertexType(long oper);

public:

    Replica(unsigned short _Nx, unsigned short _Ny, float _T, long seed);
    long DiagonalMove(float ratio);

    void Equilibrate();
    void ConstructLinks();

    vector<long>* getFirst() {return &first;}
    vector<long>* getLast()  {return &last;}
    vector<long>* getLink()  {return &links;}
    vector<long>* getVtx()   {return &vtx;}
    vector<long>* getSpin() {return &spins;}
    vector<long>* getOper() {return &sm;}
    long          getn()     {return n;}
};
#endif
