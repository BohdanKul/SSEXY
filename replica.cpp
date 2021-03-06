#include <iostream>
#include "replica.h"
//#include "communicator.h"
//#include "communicator.cpp"
#include "stdlib.h"
#include <math.h>
#include <string.h>
#include <iomanip>
#include "helper.cpp"
#include <map>
#include <list>

using namespace std;
    

/**************************************************************
*  Constructor 
**************************************************************/
Replica::Replica(unsigned short _Nx, unsigned short _Ny, float _T, long seed, int _id):
//Initialize random objects with default values
 RandomBase(seed)

{
  
    //Data structure that associates possible
    //moves for every type of vertex. There are 
    //six types of vertices, and three types of
    //moves:
    //  0 - continue straight
    //  1 - switch and continue
    //  2 - switch and reverse
    long tmp[6][2] = {  {0,1},
                        {0,1},
                        {0,2},
                        {0,2},
                        {1,2},
                        {1,2}
                     };
    memcpy(PossibleVertexMoves,tmp,sizeof tmp);   


    first.size();
    //Initialize variables
    Nx     = _Nx;
    Ny     = _Ny;
    if (Ny>1) ndim = 2;
    else      ndim = 1;

    N      = Nx*Ny;
    T    = _T;
    Beta = 1.0/T;

    //Initialize bonds
    if (Ny>1) NBonds = 2*N;
    else      NBonds = N; 
    sites.resize(NBonds+1,vector<long>(2,0));
    LatticeGeometry();
    

   //Initialize spins
    spins.resize(N,0);
    spinPart.resize(2*N, -1);
    for (vector<long>::iterator spin=spins.begin(); spin!=spins.end(); spin++) {
        *spin = pow(-1,uRandInt()%2);
    }

    first.resize(N,-1);
    last.resize(N,-1);

    //Initialize operator list
    M = round(((float) NBonds) *Beta);
    n = 0;
    sm.resize(M,0);  


    //Initialize algorithmic variables
    id = _id;
    ESteps  = 10000000;              
    estep = 1000;

}


/**************************************************************
* Increase the length of the operator list 
**************************************************************/
void Replica::AdjustM(){
           M = M+7*long(sqrt(n));
           cout << id << ": Adjusting M to "<< M << endl;
           sm.resize(M,0); 
        }


/**************************************************************
* Measure Spin Stiffness 
**************************************************************/
float Replica::MeasureSpinStiffness(){
    //Spin stiffness variables
    long WNx  = 0;            //Number of off-diagonal shifts in x-direction 
    long WNy  = 0;            //Number of off-diagonal shifts in y-direction 
    long b    = 0;
    int shift = 0;

    //Compute Winding numbers
    WNx = 0;
    WNy = 0;
    ap = spins;               //Reinitiate propagated spins state                 
    for (vector<long>::iterator oper=sm.begin(); oper!=sm.end(); oper++){
        if  ((*oper)%2==1){
            //Determine whether the winding number 
            //needs to be increased or decreased
            if  (VertexType(*oper) == 5) shift = -1;
            else                         shift = +1;

            //Determine in what direction *oper 
            //operates and change corresponding WN
            b = (long)((*oper-(*oper)%2)/2);       //Operator's bond
            if  ((b%2) and (Ny >1)) WNy += shift;
            else                    WNx += shift;
            ap[sites[b][0]] = -ap[sites[b][0]];    //Update propagated spins state
            ap[sites[b][1]] = -ap[sites[b][1]];    //Necessarily step to  correctly track
                                                   //the type of an operator
        }
    } 
    return (long(WNx/Nx)*long(WNx/Nx) + long(WNy/Ny)*long(WNy/Ny))/(2.0*Beta);        
}

 
/**************************************************************
* Measure Spin Stiffness 
**************************************************************/
long Replica::MeasureMagnetization(){
     long Magn = 0;
     for (vector<long>::iterator cspin = spins.begin(); cspin != spins.end(); cspin++)
         Magn += *cspin;

     return Magn;
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
    if (ndim==2){                                    //A 2-d lattice
         for (long y=0; y<Ny; y++){
             for (long x=0; x<Nx; x++){
                 //Define 2 bonds per site
                 sites[b+0][0] = y*Nx+x;
                 sites[b+0][1] = y*Nx+(x+1)%Nx;
                 sites[b+1][0] = y*Nx+x;
                 sites[b+1][1] = ((y+1)%Ny)*Nx+x;
                 b += 2;
             }
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
* Diagonal Monte Carlo move
**************************************************************/
long Replica::DiagonalMove()
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
               *oper = 2*b;
               n += 1;               //Update the current operator list length                                         
               if (((float) M/n)<1.25){       //If it is too close to the max (M), halt the move 
                   return 1;   
               }
            }   
         }
     //Diagonal operator
         else if (*oper%2==0){
            b = (long)(*oper/2);                                   //Bond to be removed
            AccP = ((float) M-n+1)/(((float) NBonds)*Beta*BondDiagonalEnergy(b));  //Acceptance probability
            if (uRand()<AccP){
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
           // cout << "oper = " << *oper << endl;
            for (long a=0; a<2; a++){           //Go through both sites connected by the bond 
                if (last[sites[b][a]] != -1){   //Check whether the site has been already acted upon
                    //cout << sites[b][a] <<endl; 
                    links[4*p+a]    = last[sites[b][a]];    //Construct a link between 2 legs
                    links[last[sites[b][a]]] = 4*p+a;
                }
                else{
                    first[sites[b][a]] = 4*p+a;             //Record this first occurence
                }
            //    cout << "p = " << p << endl;
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
    tedge = ap;
}



long Replica::ContinueStraight(long enLeg){
    return enLeg - sgn(enLeg-2)*(1-(sgn((enLeg%3)-1)-1));
}

long Replica::SwitchReverse(long enLeg){
    return enLeg - sgn(enLeg%2-1);
}

long Replica::SwitchContinue(long enLeg){
    return enLeg - 2*sgn(enLeg - 2);
}
/**************************************************************
* Switch leg deterministically for the loop construction
**************************************************************/
long Replica::SwitchLegDeter(long enLeg, long vtype, long p){
     
    //"Memory" of deterministic path tracing. On first encounter
    //of a vertex, it remembers the move that has been made, so
    //that on the next pass, the same move can be done.
    
    bool DetDebug = false;
    long move = -1;  //Move to make
    //Choose a random move if the vertex hasnt been visited
    if  (VtxPreviousMove[p]==-1){
        //There are two move choices for a given vertex type
        move = PossibleVertexMoves[vtype-1][uRandInt()%2];
        //Remember it
        VtxPreviousMove[p] = move;
        if (DetDebug) cout << "New ";
    }
    //Otherwise, repeat the same move
    else{
         move = VtxPreviousMove[p];
         if (DetDebug) cout << "Memory "; 
        }
    //Do the move
    long exLeg = -1;
    switch (move){
            case 0: exLeg = ContinueStraight(enLeg); break;
            case 1: exLeg = SwitchContinue(enLeg);   break;
            case 2: exLeg = SwitchReverse(enLeg);    break;   
    }
    if (DetDebug) cout << "move=" << move <<" enLeg =" << enLeg << " exLeg=" << exLeg << " type=" << vtype << " operator#=" << p << endl;
    return exLeg;
}


//long Replica::SwitchLegDeter(long enLeg, long vtype, long p){
//    
//    //Go straight if it is type 1 vertex
//    if (vtype<3){
//        //cout << "type: " << setw(4) << vtype << " leg: " << setw(4) << enLeg << endl;
//        return ContinueStraight(enLeg);
//    }
//    //Switch and reverse if it type 2 vertex
//    else{
//        return SwitchReverse(enLeg);
//    }
//}
    
/**************************************************************
* Partition edge spins based on the loop they belong to. 
**************************************************************/
void Replica::GetDeterministicLinks(){


    bool DetDebug = false;

    //Reset the main datastructure    
    fill(spinPart.begin(),spinPart.end(),-1);

    //Reset the memory of performed vertex moves
    VtxPreviousMove.assign(n,-1);

    //A map from edge legs to spins they are associated with.
    //In order to distinguish the upper and lower edge spins,
    //spin index on the upper edge is shifted by the number of spins.
    map<long,int> LegToSpin;

    //Edge or not edge leg status vector 
    vector<bool> isELeg(4*n,false);

    //Fill the map and the list
    int spin = 0;
    if (DetDebug)
        cout << endl << "Before partitionning-------------------------------" << endl;
    for (auto leg=first.begin(); leg!=first.end(); leg++){
        if (DetDebug) cout << setw(4) << *leg;
        if  (*leg!=-1){
            LegToSpin[*leg] = spin;
            isELeg[*leg] = true;
        }
        spin += 1;
    }
    if (DetDebug) cout << endl;
    

    //Continue the filling with the 2nd edge legs
    for (auto leg=last.begin(); leg!=last.end(); leg++){
        if (DetDebug)  cout << setw(4) << *leg;
        if  (*leg!=-1){
            LegToSpin[*leg] = spin;
            isELeg[*leg] = true;
        }
        spin += 1;
    }
    if (DetDebug) cout << endl;

   
    // Marked or unmarked leg status vector
    vector<bool> isMLeg(4*n,false); 
    
    long p;          //Operator index
    long leg;        //Leg index
    long nLoop = 0;  //Number of constructed loops

    list<long> path;
    LoopPaths.clear();

    //Find the head and tail spins of each loop
    for (int ispin=0; ispin!=2*N; ispin++){
        
        //cout << "new spin: " << ispin << endl; 
        //Get the leg associated with the spin
        if (ispin<N) leg = first[ispin];
        else         leg = last[ispin-N];
    
        //If the leg is inactive, the loop is trivial
        if  (leg==-1){ 
            if  (ispin<N){
                path.clear();
                LoopPaths[ispin]   = path;
                LoopPaths[ispin+N] = path;

                nLoop += 1;
                if  (DetDebug) cout << "Loop: " << nLoop << endl;
                spinPart[ispin]   = N+ispin;
                spinPart[ispin+N] = ispin;
            }
            else continue;
        }
    
        //Otherwise, follow links until we hit an edge leg
        else{
            if (DetDebug) cout << endl << "spin: "<<setw(4)<<ispin<<" leg: "<<setw(4) << leg << endl;
            //If the leg hasnt been assigned to a loop yet
            if (isMLeg[leg] == false){
                // Set its status as marked
                isMLeg[leg] = true;
                
                // Construct a new loop 
                nLoop += 1; 
                if  (DetDebug) cout << "Loop: " << nLoop << endl;
                path.clear();
                path.push_back(leg);

                //First vertex move needs to be done out of loop
                //since the last move must also be a vertex move 
                p = (long) leg/4;   //index of the corresponding operator
                leg = p*4 + SwitchLegDeter(leg%4,vtx[p],p);  
                if (DetDebug) cout << "V switch: " << setw(4) << leg;

                while (isELeg[leg] == false){
                   //Move to the leg it is connected to
                   path.push_back(leg);
                   leg = links[leg];
                   path.push_back(leg);
                   p   = (long) leg/4; 
                   if (DetDebug) cout << " L switch: " << setw(4) << leg << endl;

                    //Switch to another leg on the same vertex
                    p = (long) leg/4;   //index of the corresponding operator
                    leg = p*4 + SwitchLegDeter(leg%4,vtx[p],p);  
                    if (DetDebug) cout << "V switch: " << setw(4) << leg;
               
                //Stop if  we have reached a leg at an edge
                }  
                isMLeg[leg] = true;
                path.push_back(leg);
            
                //Store loop's head and tail spins.
                spinPart[ispin] = LegToSpin[leg]; 
                spinPart[LegToSpin[leg]] = ispin; 
                LoopPaths[ispin]   = path;
                LoopPaths[LegToSpin[leg]] = path;
            }
        }
        if  (DetDebug){
            cout << "Spin partition: " << endl;
            long i=0;
            for (auto spin=spinPart.begin(); spin!=spinPart.end(); spin++){
                cout << setw(4) << *spin;
                i += 1;
                if  (i == N) cout << endl; 
            }
            cout << endl;
        }
    }            
} 
