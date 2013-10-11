#include <iostream>
#include "ssexy.h"
#include "stdlib.h"
#include <math.h>
#include <string.h>

using namespace std;
    

/**************************************************************
*  Constructor 
**************************************************************/
SSEXY::SSEXY(unsigned short _Nx, unsigned short _Ny, float _T, unsigned int seed):
//Initialize random objects with default values
eng(1),uReal(0,1),uInt(0,65535),uRandInt(eng,uInt),uRand(eng,uReal) 

{
    int tmp[6][4] = {   {-1,-1,-1,-1},
                        { 1, 1, 1, 1},
                        { 1,-1,-1, 1},
                        {-1, 1, 1,-1},
                        {-1, 1,-1, 1},
                        { 1,-1, 1,-1}
                    };
    memcpy(LegSpin,tmp,sizeof tmp);
   
     //Initialize ranndomness with the seed
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

    first.resize(N,-1);
    last.resize(N,-1);

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
* Measuring diagonal contribution to energy from a bond
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
     ap = spins;
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
            ap[sites[b][0]] = -ap[sites[b][0]]; 
            ap[sites[b][1]] = -ap[sites[b][1]]; 
         }
    }
    
return 0;
}



/**************************************************************
* Determine type of a vertex for a particular operator acting 
* on a given bond based on its type and current propagated state. 
**************************************************************/
unsigned int SSEXY::VertexType(unsigned int b, int oper)
{
    int s0 = ap[sites[b][0]];  
    int s1 = ap[sites[b][1]];
    if  (oper%2 == 0){
        if  ((s0 ==-1) and (s1 ==-1))
            return 1;
        if  ((s0 ==-1) and (s1 == 1))
            return 2;
        if  ((s0 == 1) and (s1 ==-1))
            return 3;
        if  ((s0 == 1) and (s1 == 1))
            return 6;
    }   
    else{
        if  ((s0 ==-1) and (s1 == 1))
            return 4;
        if  ((s0 == 1) and (s1 ==-1))
            return 5;
    }      
}

int sgn(int val) {
    return ((0<=val) - (val<0));
}


/**************************************************************
* Off-diagonal Monte Carlo move
**************************************************************/
void SSEXY::OffDiagonalMove()
{
    int p=-1;                                       //Index of the current operator            
    unsigned int b=0;                               //Bond index
    vector<unsigned int> links(4*n,-1);             //Linked vertex list
    vector<unsigned int> vtx(n,-1);                 //Vertex type
    //Reinitiate the propagated spins state
    ap = spins;
 
    //Construct linked vertex list (link vector structure)
    //Coordiantes of a leg on a particular operator can be encoded 
    //by one number: 4*operator_index+leg_index (1 of 4). By construction, 
    //the value of the links element for a particular index (link[index]), 
    //encodes the leg coordinates to which the leg encoded by the 
    //index-value is linked.  
    for (vector<int>::iterator oper=sm.begin(); oper!=sm.end(); oper++) {
        //Ignore unit operator
        if (*oper != 0){
            //Determine operator's bond based on operator's value
            b = (unsigned int)((*oper-(*oper)%2)/2);
            //Increase the operator counter
            p++;
            //For each of 2 site bonds 
            for (int a=0; a<2; a++){
                //If it is not for the first time that an operator
                //acts on this site. We can construct a links.
                if (last[sites[b][a]] != -1){
                    links[4*p+a]    = last[sites[b][a]];
                    links[last[sites[b][a]]] = 4*p+a;
                }
                //Otherwise, just record this first occurence
                else{
                    first[sites[b][a]] = 4*p+a;
                }
                //Update what is the last leg acting on this site
                last[sites[b][a]]  = 4*p+a+2;
            }
           
            //Record what type of vertex is it
            vtx[p] = VertexType(b,*oper);

            //Update the propagated spins state if it is an off-diagonal operator
            if  (*oper%2 == 1){
                ap[sites[b][0]] = -ap[sites[b][0]]; 
                ap[sites[b][1]] = -ap[sites[b][1]]; 
            }
            
        }

    }

    //Construct links across the boundary in the p-expansion
    for (int lspin=0; lspin<last.size(); lspin++){
        if (last[lspin] != -1){
            links[last[lspin]]  = first[lspin];
            links[first[lspin]] = last[lspin];
        }
    }


//------------------------------------------------------------------------------------------------
    unsigned int j0;            //Loop entrance leg
    unsigned int j;             //Current leg
    pair<int,int> legtype;      //Struct with the leg, and operator type   
 
    //Construct Nl number of loops 
    for (int i=0; i!=Nl; i++) {
        j0 = uRandInt()%(4*n);    //Pick a random loop entrance leg among all possible legs
        j  = j;
        //Construct an operator-loop
        do  {
                p = (int) j/4;                      //Current operator index
                legtype = SwitchLeg(j%4,vtx[p]);    //Get the next leg and the new operator's type
                j       = legtype.first +4*p;       //Move to the next leg
                vtx[p]  = legtype.second;           //Update the type of the operator
                if  (j == j0) break;                 //If the loop is closed, we are done
                else j = links[j];                  //Else move to the next linked leg

        } while  (j !=j0);                          //Another way to close the loop
    }
    
//------------------------------------------------------------------------------------------------
    //Map back the changes to the operator list
    p = -1;
    int leg;
    for (vector<int>::iterator oper=sm.begin(); oper!=sm.end(); oper++) 
        if (*oper != 0){
            b = (unsigned int)((*oper-(*oper)%2)/2);     //Determine operator's bond
            p++;                                         //Increase the operator counter
            *oper = 2*b + (sgn(vtx[p]-4)+1)%2;           //The new operator bond is unchanged
        }                                                //but its type could have changed
        
    //Map back the changes to the spins state
    for (int i=0; i!=N; i++)
        if  ((first[i]==-1) and (uRand()<0.5))             //If the spin is not acted by an operator 
            spins[i] = -spins[i];                        //Flip it with probability 1/2
        else{                                            //Otherwise:               
            p   = (int) first[i]/4;                      //Get the first operator that acts on it   
            leg = first[i]%4;                            //And the leg that the spin is      
            spins[i] = LegSpin[vtx[p]-1][leg];           //Use the pregenerated leg/operator -> spins
                                                         //map to establish what is the new spin  
            }
     
}


//#################################################################
// There are 4 different types of moves for pair-interacting bonds.
// Below they are implemented as seperate functions. For a given
// enLeg, each function returns the exit leg after the move is done.
// Leg labelling convention : bottom left leg is 0,  the label # is 
// increased by one in the counter-clockwise direction.     
//#################################################################


int Bounce(unsigned int enLeg){
    return 0;
}

int ContinueStraight(unsigned int enLeg){
    return enLeg - sgn(enLeg-2)*(1-(sgn((enLeg%3)-1)-1));
}

unsigned int SwitchReverse(unsigned int enLeg){
    return enLeg - sgn(enLeg%2-1);
}

unsigned int SwitchContinue(unsigned int enLeg){
    return enLeg - 2*sgn(enLeg - 2);
}
//#################################################################
//#################################################################




int LegSpin(unsigned int leg, unsigned int vtype){
    
}


pair<int,int> SSEXY::SwitchLeg(unsigned int enLeg, unsigned int vtype){
    int exLeg   = -1;
    int newtype = -1;
    switch (vtype){
        case 1:
                    if  (uRand()<0.5){
                            if ((enLeg==1) or (enLeg==3)) newtype = 5;
                            else newtype = 6;
                            exLeg   = SwitchContinue(enLeg);
                    }                
                    else{
                            if ((enLeg==0) or (enLeg==3)) newtype = 3;
                            else newtype = 4;
                            exLeg   = ContinueStraight(enLeg);
                    }
                    break;
        case 2:
                    if  (uRand()<0.5){
                            if ((enLeg==1) or (enLeg==3)) newtype = 6;
                            else newtype = 5;
                            exLeg   = SwitchContinue(enLeg);
                    }                
                    else{
                            if ((enLeg==0) or (enLeg==3)) newtype = 4;
                            else newtype = 3;
                            exLeg   = ContinueStraight(enLeg);
                    }
                    break;
        case 3:
                    if  (uRand()<0.5){
                            if ((enLeg==0) or (enLeg==3)) newtype = 1;
                            else newtype = 2;
                            exLeg   = ContinueStraight(enLeg);
                    }                
                    else{
                            if ((enLeg==0) or (enLeg==1)) newtype = 5;
                            else newtype = 6;
                            exLeg   = SwitchReverse(enLeg);
                    }
                    break;
        case 4:
                    if  (uRand()<0.5){
                            if ((enLeg==0) or (enLeg==3)) newtype = 2;
                            else newtype = 1;
                            exLeg   = ContinueStraight(enLeg);
                    }                
                    else{
                            if ((enLeg==0) or (enLeg==1)) newtype = 6;
                            else newtype = 5;
                            exLeg   = SwitchReverse(enLeg);
                    }
                    break;
        case 5:
                    if  (uRand()<0.5){
                            if ((enLeg==0) or (enLeg==1)) newtype = 3;
                            else newtype = 4;
                            exLeg   = SwitchReverse(enLeg);
                    }                
                    else{
                            if ((enLeg==0) or (enLeg==2)) newtype = 2;
                            else newtype = 1;
                            exLeg   = SwitchContinue(enLeg);
                    }
                    break;
             
        case 6:
                    if  (uRand()<0.5){
                            if ((enLeg==0) or (enLeg==1)) newtype = 4;
                            else newtype = 3;
                            exLeg   = SwitchReverse(enLeg);
                    }                
                    else{
                            if ((enLeg==0) or (enLeg==2)) newtype = 1;
                            else newtype = 2;
                            exLeg   = SwitchContinue(enLeg);
                    }
                    break;
 
    }
return pair <int,int> (exLeg,newtype);
}

/*############################################
# Future Code
##############################################

unsigned int DeterministicSwitchLeg(unsigned int j, unsigned int vtype){
    switch (vtype){
        case 1:
            nettype = 6 - (enLeg%2);
            exLeg   = enLeg - 2*sgn(enLeg - 2);
            break;
        case 2:
            newtype = 5 + (enLeg%2);
            exLeg   = enLeg - 2*sgn(enLeg - 2);
            break;
        case 3:
            newtype = 5;
            exLeg = enLeg + sgn(enLeg%2-1);
        case 4:
            newtype = 6;
            exLeg = enLeg + sgn(enLeg%2-1); 
        case 5:
            
        case 6:

    }
}
*/



int main()
{
    SSEXY ssexy(5,1,2,1);
    ssexy.Equilibrate();
    return 0;
}
