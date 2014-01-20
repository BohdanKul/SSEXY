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
    unsigned short Nx,Ny;        //Lattice spatial dimensions
    unsigned short ndim;         //Lattice dimension
    long   NBonds;               //Number of bonds
    long N;                      //Number of spin sites
    float T;                             //Temperature
    float Beta;                          //Inverse temperature
    vector <vector<long>> sites; //Sites connected by bonds  
    vector <long> spins;                  //Spins array
    vector <long> ap;                     //Propagated spins state

    //Algorithmic variables
    int  id;                    
    long ESteps;                //Number of equilibration steps    
    long estep;
    long M;                     //Operators list size
    long n;                     //Non-I operators list size
    vector<long> sm;            //Operator vector list
    vector<long> first;         //First and last leg coordinates
    vector<long> last;          //for a particular spin site
    vector<long> links;         //Linked vertex list
    vector<long> vtx;           //Vertex type
    vector<long> spinPart;      //Partion keeping the label of a loop
                                //that a particular spin belongs to
    long nLoop;                 //Number of distinct loops in a partition
    long nfLoop;                //Number of attempted loops during partitioning
    bool Debug;

 
   
    //Methods
    float BondDiagonalEnergy(long b);
    void  LatticeGeometry();
    long  VertexType(long oper);
    void  getSpinNeighs(long spin, vector<long>& sneighs);

public:

    Replica(unsigned short _Nx, unsigned short _Ny, float _T, long seed, int _id);
    long DiagonalMove();

    long  MeasureMagnetization();
    float MeasureSpinStiffness();

    long PartitionSpins();
    void RelabelLoop(long initSpin,long olabel,long nlabel);

    void ConstructLinks();
    void AdjustM();

    vector<vector<long>> getSites()  {return sites;}
    vector<long>* getFirst()  {return &first;}
    vector<long>* getLast()   {return &last;}
    vector<long>* getLink()   {return &links;}
    vector<long>* getVtx()    {return &vtx;}
    vector<long>* getSpin()   {return &spins;}
    vector<long>* getOper()   {return &sm;}
    vector<long>* getPart()   {return &spinPart;}
    long          getn()      {return n;}
    long          getM()      {return M;}
    long          getnLoops() {return nLoop;}
    long          getnfLoops(){return nfLoop;}
    void          setn(long _n) {n=_n;}
    void          setM(long _M) {M=_M;}
    void          setSpinAt(long i, long n){spins[i]=n;}
};
#endif
