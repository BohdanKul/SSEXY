#include <iostream>
#include "replica.h"
//#include "communicator.h"
//#include "communicator.cpp"
#include "stdlib.h"
#include <math.h>
#include <string.h>
#include <iomanip>
#include "helper.cpp"

using namespace std;
    

/**************************************************************
*  Constructor 
**************************************************************/
Replica::Replica(unsigned short _Nx, unsigned short _Ny, float _T, long seed):
//Initialize random objects with default values
communicator(_Nx,_Ny,_T), RandomBase(seed)

{
  


    first.size();
    //Initialize variables
    Nx     = _Nx;
    Ny     = _Ny;
    N      = Nx*Ny;
    T    = _T;
    Beta = 1.0/T;
    cout << "Beta = " << Beta << endl;

    //Initialize bonds
    if (Ny>1) NBonds = 2*N;
    else      NBonds = N; 
    sites.resize(NBonds+1,vector<long>(2,0));
    LatticeGeometry();
    

   //Initialize spins
    spins.resize(N,0);
    for (vector<long>::iterator spin=spins.begin(); spin!=spins.end(); spin++) {
        *spin = pow(-1,uRandInt()%2);
        //cout << *spin << " ";
    }

    first.resize(N,-1);
    last.resize(N,-1);

    //Initialize operator list
    M = round(((float) NBonds) *Beta*1.5);
    n = 0;
    sm.resize(M,0);  

    //Initialize algorithmic variables
    ESteps  = 10000000;              
    estep = 1000;
    Debug = false;

    //Initialize communicator class
    string eHeader = boost::str(boost::format("#%15s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s") %"n"%"dn"%"E"%"dE"%"Legs"%"dLegs"%"Magn"%"dMagn"%"M"%"dM"%"rho_s"%"drho_s");
    *communicator.stream("estimator") << eHeader <<endl;    
    
    for (auto bond = ++sites.begin(); bond != sites.end(); bond++)
        *communicator.stream("bond") << boost::str(boost::format("%6d%6d") %bond->at(0) %bond->at(1) ) ;
    *communicator.stream("bond") << endl;
}



/**************************************************************
* Equilibration. Parameters adjustment 
**************************************************************/
void Replica::Equilibrate(){

    long i;
    long  totalvLegs = 0;     //Average number of visited legs
    long Cumn = 0;            //Cummulative number of non-identity operators  
    long CumM = 0;            //Cummulative magnetization  
    long totalMagn = 0;       //Average magnetication  
    long CumMagn  = 0;        //Cummulative magnetization  
    long CumvLegs = 0;        //Cummulative number of visited legs  
    float E = 0;              //Average energy  

    //Spind stiffness variables
    long WNx  = 0;            //Number of off-diagonal shifts in x-direction 
    long WNy  = 0;            //Number of off-diagonal shifts in y-direction 
    long b      = 0;          //Bond index
    long shift  = 1;          //Takes only 2 values: +1 or -1 depending on an off-diaogonal operator
    long CumSS  = 0;          //Cummulative spin stiffness

    //Main loop
    for (long step=0; step<ESteps; step++){
        if (Debug){
            //Spins output
            for (vector<long>::iterator spin=spins.begin(); spin!=spins.end(); spin++) 
                *communicator.stream("spin") << boost::str(boost::format("%6d") %*spin) ;
            *communicator.stream("spin") << endl;
            //Operator output 
            for (vector<long>::iterator oper=sm.begin(); oper!=sm.end(); oper++) 
                if (not(*oper == 0))
                   *communicator.stream("operator") << boost::str(boost::format("%6d") %*oper) ;
            *communicator.stream("operator") << endl;
        } 
        //Diagonal move 
        if (DiagonalMove(1.5) == 1){
           M = long(1.2*M);
           sm.resize(M,0); 
        }
        Cumn += n;
        CumM += M;

        if (Debug){
        //Operator output 
            for (vector<long>::iterator oper=sm.begin(); oper!=sm.end(); oper++) 
                //if (not(*oper == 0))
                    *communicator.stream("operator") << boost::str(boost::format("%6d") %*oper) ;
            *communicator.stream("operator") << endl;
        }

        
        //Compute Magnetization   
        totalMagn = 0;
        for (vector<long>::iterator cspin = spins.begin(); cspin != spins.end(); cspin++)
            totalMagn += *cspin;
        CumMagn += totalMagn;

        Debug = false;
        //Compute Winding numbers
        WNx = 0;
        WNy = 0;
        ap = spins;                             //Reinitiate the propagated spins state                 
        for (vector<long>::iterator oper=sm.begin(); oper!=sm.end(); oper++){
            if  ((*oper)%2==1){
                //Determine whether the winding number needs to be increased or decreased
                if  (VertexType(*oper) == 5) shift = -1;
                else                         shift = +1;
                //Determine in what direction *oper operates and change corresponding WN
                b = (long)((*oper-(*oper)%2)/2);                      //Operator's bond
                if  (b%2) WNy += shift;
                else      WNx += shift;
                ap[sites[b][0]] = -ap[sites[b][0]];         //Update the propagated spins state
                ap[sites[b][1]] = -ap[sites[b][1]]; 
            }
        } 
        CumSS += long(WNx/Nx)*long(WNx/Nx) + long(WNy/Ny)*long(WNy/Ny);        

        //Estimators output         
        if (step%estep == 0){
           E = -((float) Cumn/((float) estep*N))/Beta + (float) 1;
           *communicator.stream("estimator") << boost::str(boost::format("%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E") %(Cumn/(1.0*estep)) %0.0 %E %0.0 %(CumvLegs/(1.0*estep)) %0.0 %(CumMagn/(1.0*estep)) %0.0 %(CumM/(1.0*estep)) %0.0 %(CumSS/(2.0*estep*Beta)) %0.0)<< endl;
           //cout << "n = " << setw(7) << (float) Cumn/((float) estep) << " E = " << setw(7) << E << " M = " << setw(9) << (float) CumM/estep << " Legs = " << setw(6) << (float) CumvLegs/estep << " Magnetization " << setw(6) << (float) CumMagn/estep << endl;
           Cumn = CumM = CumvLegs = CumMagn = CumSS = 0;
           }      

         
     }
}

/***************************************************************
* Building bonds structure via definition of spin they act upon 
***************************************************************/
void Replica::LatticeGeometry()
/*
Define a 2-d array to store the start (0th row) and the end (1st row) 
of the b'th pair bond  for a  2-d periodic square lattice.
Below, the left bottom site initiates 2 bonds:
.
|
. __ .

*/
{
    long b=1;                             //Bond index
                                                   //Required to start from a value > 0 
    if (Ny >1){                                    //A 2-d lattice
         for (long y=0; y<Ny; y++){
             for (long x=0; x<Nx; x++){
                 //Define 2 bonds per site
                 sites[b+0][0] = y*Nx+x;
                 sites[b+0][1] = y*Nx+(x+1)%Nx;
                 sites[b+1][0] = y*Nx+x;
                 sites[b+1][1] = ((y+1)%Ny)*Nx+x;
                 b += 2;
             //    cout << site2[b]<< " ";
             }
             //cout << endl;
         } 
    }
    else                                           //A 1-d chain
         for (long x=0; x<Nx; x++){
             //Define 1 bond per site
             b = x+1;
             sites[b][0] = x;
             sites[b][1] = (x+1)%Nx;
             }
}

float Replica::BondDiagonalEnergy(long b){
      return (float) 0.5;
}

/**************************************************************
* Diogonal Monte Carlo move
**************************************************************/
long Replica::DiagonalMove(float ratio)
{
     long b=0;    //Bond index
     float AccP = 0;       //Acceptance probability
     float uRan;
     ap = spins;
     for (vector<long>::iterator oper=sm.begin(); oper!=sm.end(); oper++) {
     //Identity operator
         if (*oper==0){
            b = uRandInt()%NBonds+1;                                                   //Bond to be created 
            AccP = ((float) NBonds)*Beta*BondDiagonalEnergy(b) / ((float) (M-n));     //Acceptance probability
            uRan =  uRand();
            if (uRan<AccP){
                //cout << "I -> D. Operator = " << *oper << endl;
               *oper = 2*b;
               n += 1;               //Update the current operator list length                                         
               if (((float) M/n)<ratio){       //If it is too close to the max (M), halt the move 
                   cout << "M/n = " << (float) M/n << " < " << ratio << endl;
                   return 1;   
               }
            }   
         }
     //Diagonal operator
         else if (*oper%2==0){
            b = (long)(*oper/2);                                   //Bond to be removed
            AccP = ((float) M-n+1)/(((float) NBonds)*Beta*BondDiagonalEnergy(b));  //Acceptance probability
            if (uRand()<AccP){
               //cout << "D -> I. Operator = " << *oper  <<endl;
               *oper = 0;
               n -= 1;  
            }
         }
     //Off-diagonal operator
         else{
//            b = (long)((*oper-1)/2);                      //Bond being acted on
//            ap[sites[b][0]] = -ap[sites[b][0]];           //Flip spins connected by the b'th bond
//            ap[sites[b][1]] = -ap[sites[b][1]]; 
         }
    }
    
return 0;
}



/**************************************************************
* Determine type of a vertex for a particular operator acting 
* on a given bond based on its type and current propagated state. 
**************************************************************/
long Replica::VertexType(long oper)
{
    long b = (long)((oper-oper%2)/2);    //Determine operator's bond based on operator's value
    long s0 = ap[sites[b][0]];  
    long s1 = ap[sites[b][1]];
    //cout << oper << endl;
    if  (oper%2 == 0){
        if  ((s0 ==-1) and (s1 ==-1))
            return 1;
        if  ((s0 == 1) and (s1 == 1))
            return 2;
        if  ((s0 == 1) and (s1 ==-1))
            return 3;
        if  ((s0 == -1) and (s1 == 1))
            return 4;
    }   
    else{
        if  ((s0 ==-1) and (s1 == 1))
            return 5;
        if  ((s0 == 1) and (s1 ==-1))
            return 6;
    }     
    cout << "Impossible configuration in VertexType" << endl;
    return -1; 
}



/**************************************************************
* Construct links in preparation for an off-diagonal move
**************************************************************/
void Replica::ConstructLinks()
{
    long p=-1;                              //Index of the current operator            
    long b=0;                               //Bond index
   
    //Reinitiate necessary data structures    
    ap = spins;                             //propagated spins state                 
    fill(first.begin(),first.end(),-1);     //spin's 1st  leg
    fill(last.begin(),last.end(),-1);       //spin's last leg
    links.assign(4*n,-1);                   //linked vertices
    vtx.assign(n,-1);                       //vertex type


    //Construct linked vertex list (link vector structure)
    //Coordiantes of a leg on a particular operator can be encoded 
    //by one number: 4*operator_index+leg_index (1 of 4). By construction, 
    //the value of the links element for a particular index (link[index]), 
    //encodes the leg coordinates to which the leg encoded by the 
    //index-value is linked.  
    for (vector<long>::iterator oper=sm.begin(); oper!=sm.end(); oper++) {
        if (*oper != 0){                        //Ignore unit operator
            b = (long)((*oper-(*oper)%2)/2);    //Determine operator's bond based on operator's value
            p++;                                //Increase the operator counter
            for (long a=0; a<2; a++){           //Go through both sites connected by the bond 
                if (last[sites[b][a]] != -1){   //Check whether the site has been already acted upon
                    //cout << sites[b][a] <<endl; 
                    links[4*p+a]    = last[sites[b][a]];    //Construct a link between 2 legs
                    links[last[sites[b][a]]] = 4*p+a;
                }
                else{
                    first[sites[b][a]] = 4*p+a;             //Record this first occurence
                }
                if (a==0) last[sites[b][a]]  = 4*p+3;       //Update what is the last leg acting on this site
                else      last[sites[b][a]]  = 4*p+2;
            }
           
            vtx[p] = VertexType(*oper);                   //Record what type of vertex is it
            
            if  (*oper%2 == 1){                             //If it is an off-diagonal operator
                ap[sites[b][0]] = -ap[sites[b][0]];         //Update the propagated spins state
                ap[sites[b][1]] = -ap[sites[b][1]]; 
            }
            
        }

    }
} 








