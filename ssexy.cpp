#include <iostream>
#include "ssexy.h"
#include "stdlib.h"
#include <math.h>

using namespace std;
    

/**************************************************************
*  Constructor 
**************************************************************/
SSEXY::SSEXY(unsigned short _Nx, unsigned short _Ny, float _T, unsigned int seed):
//Initialize random objects with default values
eng(1),uReal(0,1),uInt(0,10240),uRandInt(eng,uInt),uRand(eng,uReal) 

{
    //Reinitialize randomness with the seed
    uRand.engine().seed(seed);
    uRand.distribution().reset();
    uRandInt.engine().seed(seed);
    uRandInt.distribution().reset();

    //Initialize variables
    Nx     = _Nx;
    Ny     = _Ny;
    N      = Nx*Ny;
    T    = _T;
    Beta = 1/T;

    //Initialize bonds
    NBonds = 2*N; 
    sites.resize(NBonds,vector<unsigned int>(2,0));
    LatticeGeometry();
    
    //Initialize spins
    spins.resize(N,0);
    for (vector<int>::iterator spin=spins.begin(); spin!=spins.end(); spin++) {
        *spin = pow(-1,1+uRandInt()%2);
        //cout << *spin << " ";
    }

    //Initialize operator list
    M = round(Beta*NBonds*1.5);
    n = 0;
    sm.resize(M,0);  

    //Initialize algorithmic variables
    ESteps = 15;              
}



/**************************************************************
* Equilibration. Parameters adjustment 
**************************************************************/
void SSEXY::Equilibrate(){

    for (int step=0; step<ESteps; step++){
         cout << "Operator list: \n";   
         for (vector<int>::iterator oper=sm.begin(); oper!=sm.end(); oper++) {
             cout << *oper << " ";
         }
         cout << endl;     
         cout << "Spins: \n";
         for (vector<int>::iterator spin=spins.begin(); spin!=spins.end(); spin++) {
             cout << *spin << " ";
         }
         cout << endl;     

         if (DiagonalMove(1.5) == 1){
            cout << "Resizing \n"; 
            M = int(1.2*M);
            sm.resize(M,0); 
        }
   }
}



/**************************************************************
* Measuring diogonal contribution to energy from a bond
**************************************************************/
float SSEXY::BondDiagonalEnergy(unsigned int b){
    return 1.0;
}


/***************************************************************
* Building bonds structure via definition of spin they act upon 
***************************************************************/
void SSEXY::LatticeGeometry()
/*
Define 2 arrays that remember where a b'th pair bond starts from (site1)
and where it goes to (site 2) for a  2-d periodic square lattice.
Below, the left bottom site initiates 2 bonds:
.
|
. __ .
*/
{
     unsigned int b=1;                             //Bond index
     for (int y=0; y<Ny; y++){
         for (int x=0; x<Nx; x++){
             sites[b+0][0] = y*Nx+x;
             sites[b+0][1] = y*Nx+(x+1)%Nx;
             sites[b+1][0] = y*Nx+x;
             sites[b+1][1] = ((y+1)%Ny)*Nx+x;
         //    cout << site2[b]<< " ";
         }
         //cout << endl;
         b += 2;
     }      
}


/**************************************************************
* Diogonal Monte Carlo move
**************************************************************/
int SSEXY::DiagonalMove(float ratio)
{
     unsigned int b=0;     //Bond index
     unsigned int nnew=0;  //Counter of the net added diagonal operators number
     float AccP = 0;       //Acceptance probability
     for (vector<int>::iterator oper=sm.begin(); oper!=sm.end(); oper++) {
         //Identity operator
         if (*oper==0){
            b = uRandInt()%NBonds;                                //Bond to be created 
            AccP = NBonds*Beta*BondDiagonalEnergy(b) / (M-n);     //Acceptance probability
            if (uRand()<AccP){
               *oper = 2*b;
               //Update the current operator list length
               n += 1;                                          
               //If it is too close to the max (M), halt the move 
               if (M/n<ratio){
                  cout << "Ratio "<< n << endl;
                  return 1;   
               }
            }   
         }
         //Diagonal operator
         else if (*oper%2==0){
            b = (unsigned int)(*oper/2);                          //Bond to be removed
            AccP = (M-n+1) / (NBonds*Beta*BondDiagonalEnergy(b)); //Acceptance probability
            if (uRand()<AccP){
               *oper = 0;
               n -= 1;  
            }
         }
         //Off-diagonal operator
         else{
            b = (unsigned int)((*oper-1)/2);                      //Bond being acted on
            //Flip the corresponding spins
            spins[sites[b][0]] = -spins[sites[b][0]]; 
            spins[sites[b][1]] = -spins[sites[b][1]]; 
         }
    }
    
return 0;
}


/**************************************************************
* Off-diagonal Monte Carlo move
**************************************************************/
void SSEXY::OffDiagonalMove()
{
    vector<unsigned int> first(N,-1), last(N,-1);
    vector<unsigned int> link(4*n,-1);
    int p=-1;
    unsigned int b=0;     //Bond index
    for (vector<int>::iterator oper=sm.begin(); oper!=sm.end(); oper++) {
        if (*oper != 0){
            b = (unsigned int)((*oper-(*oper)%2)/2);
            p++;
            for (int a=0; a<2; a++){
                if (last[sites[b][a]] != -1){
                    link[4*p+a]    = last[sites[b][a]];
                    link[last[sites[b][a]]] = 4*p+a;
                }
                else{
                    first[sites[b][a]] = 4*p+a;
                }
                last[sites[b][a]]  = 4*p+a+2;
            }
        }
    }
    for (int lspin=0; lspin<last.size(); lspin++){
        if (last[lspin] != -1){
            link[last[lspin]]  = first[lspin];
            link[first[lspin]] = last[lspin];
        }
    }
}

int main()
{
    SSEXY ssexy(5,1,2,1);
    ssexy.Equilibrate();
    return 0;
}
