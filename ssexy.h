#include <vector>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
//#include "communicator.h"
#include "communicator.cpp"

using namespace std;

//Define data-types of random objects 
typedef boost::mt19937 t_eng;                 //Mersenne twister
typedef boost::uniform_real<float> t_uReal;   //Uniform real distribution
typedef boost::uniform_int<int>    t_uInt;    //Uniform longeger distribution

//Workhorse
class SSEXY
{
private:
    long LegSpin[6][4];

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
    long Nloops;                 //Number of loops per MC move
    long M;                      //Upper operator list size
    long n;                      //Current operator list size
    vector <long> VisitLegs;              //
    vector <long> sm;                     //Operator vector list
    vector<long> first;          //First and last leg coordinates
    vector<long> last;           //for a particular spin site
    vector <long> NvisitedLegs;  //Number of visited legs per MC step

    bool Debug;

    //File management
    Communicator communicator;
 
    //Random generation objects
    t_eng   eng;                 //Mersenne twister
    t_uReal uReal;               //Uniform real distribution
    t_uInt  uInt;                //Uniform distribution 
    boost::variate_generator<t_eng&, t_uInt > uRandInt; //Generator of longegers
    boost::variate_generator<t_eng&, t_uReal> uRand;    //Generator of real number
    
    //Methods
    float BondDiagonalEnergy(long b);
    void  LatticeGeometry();
    long  VertexType(long oper);

public:
    SSEXY(unsigned short _Nx, unsigned short _Ny, float _T, long seed);
    void Equilibrate();
    void OffDiagonalMove();
    long DiagonalMove(float ratio);
    pair<long,long> SwitchLeg(long leg, long vtype);
};
