#include <vector>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;

//Define data-types of random objects 
typedef boost::mt19937 t_eng;                 //Mersenne twister
typedef boost::uniform_real<float> t_uReal;   //Uniform real distribution
typedef boost::uniform_int<int>    t_uInt;    //Uniform integer distribution

//Workhorse
class SSEXY
{
private:
    int LegSpin[6][4];

    unsigned short Nx,Ny;                //Lattice dimensions
    unsigned int   NBonds;               //Number of bonds
    unsigned int M;                      //Upper operator list size
    unsigned int n;                      //Current operator list size
    unsigned int N;                      //Number of spin sites
    unsigned int Nl;                     //Number of loops per MC step
    float T;                             //Temperature
    float Beta;                          //Inverse temperature
    vector <vector<unsigned int>> sites; //Sites connected by bonds  
    vector <int> sm;                     //Operator vector list
    vector <int> spins;                  //Spins array
    vector <int> ap;                     //Propagated spins state
    vector<unsigned int> first;    //First and last leg coordinates
    vector<unsigned int> last;    //for a particular spin site

    //Algorithmic variables
    unsigned int ESteps;         //Number of equilibration steps    
    
    //Random objects
    t_eng   eng;                 //Mersenne twister
    t_uReal uReal;               //Uniform real distribution
    t_uInt  uInt;                //Uniform distribution 
    boost::variate_generator<t_eng&, t_uInt > uRandInt; //Generator of integers
    boost::variate_generator<t_eng&, t_uReal> uRand;    //Generator of real number
    
    void  LatticeGeometry();
    float BondDiagonalEnergy(unsigned int b);
    unsigned int  VertexType(unsigned int b, int oper);

public:
    SSEXY(unsigned short _Nx, unsigned short _Ny, float _T, unsigned int seed);
    void Equilibrate();
    void OffDiagonalMove();
    int DiagonalMove(float ratio);
    pair<int,int> SwitchLeg(unsigned int leg, unsigned int vtype);
};
