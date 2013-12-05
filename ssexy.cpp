#include <iostream>
#include "ssexy.h"
#include "replica.cpp"
//#include "communicator.h"
//#include "communicator.cpp"
#include "stdlib.h"
#include <math.h>
#include <string.h>
#include <iomanip>
#include "helper.cpp"

using namespace std;



//**************************************************************************
SSEXY::SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, long seed):
RandomBase(seed)
{
    long tmp[6][4] = {  {-1,-1,-1,-1},
                        { 1, 1, 1, 1},
                        { 1,-1,-1, 1},
                        {-1, 1, 1,-1},
                        {-1, 1,-1, 1},
                        { 1,-1, 1,-1}
                    };
    memcpy(LegSpin,tmp,sizeof tmp);

    // Initialize replicas
    r = _r; 
    for(int i=0; i!=r; i++){
        Replicas.push_back(new Replica(_Nx,_Ny,_T,seed));
    }
    Debug = false;
    Nloops = 1;    

    // Initialize data structures
    firsts.resize(r,NULL);
    lasts .resize(r,NULL);
    links.resize(r,NULL);
    spins.resize(r,NULL);
    vtxs.resize(r,NULL);
    sms.resize(r,NULL);
    shifts.resize(r,0);
    ns.resize(r,0);
    NvisitedLegs.resize(Nloops,0);
    
    // Initialize region A
    Aregion = {};
} 
        

//**************************************************************************
int SSEXY::MCstep()
{

    //----------------------------------------------------------------------        
    // Perform diagonal update in each replica
    // while keeping track of relevant structures 
    //----------------------------------------------------------------------        
    nTotal = 0;
    for (int j=0; j!=r; j++){
        Replicas[j]->DiagonalMove(1.5);
        Replicas[j]->ConstructLinks();
        firsts[j] = Replicas[j]->getFirst();        
        lasts[j]  = Replicas[j]->getLast();        
        links[j]  = Replicas[j]->getLink();        
        spins[j]  = Replicas[j]->getSpin();      
        vtxs[j]   = Replicas[j]->getVtx();      
        sms[j]    = Replicas[j]->getOper();      
        ns[j]     = Replicas[j]->getn();
        if (j>0)
            shifts[j] = shifts[j-1] + ns[j-1];
        cout << "n = " << ns[j] << endl;
        nTotal   += ns[j];    
    } 
    
    //----------------------------------------------------------------------        
    // Shift the values of the leg coordinates by the 
    // total number of legs in previous replicas
    //----------------------------------------------------------------------        
    for (int j=1; j!=r; j++){
        for (auto link=links[j]->begin(); link!=links[j]->end(); link++){
            if  (*link != -1){
                *link += shifts[j]*4;
            }
        }
        for (int k=0; k!=firsts[0]->size(); k++){
            if  (firsts[j]->at(k) != -1){
                 firsts[j]->at(k) += shifts[j]*4;
            }
            if  (lasts[j]->at(k) != -1){
                lasts[j]->at(k)   += shifts[j]*4;
            }
        }
    }
    
    for (int j=0; j!=r; j++){
        cout << "=============================================" << endl;
        cout << "link1" << endl;
        cout << "=============================================" << endl;
        for (auto link=links[0]->begin(); link!=links[0]->end(); link++){
            cout << *link << " ";
        }
        cout << endl;
    }

    for (int j=0; j!=r; j++){
        cout << "=============================================" << endl;
        cout << "first1" << endl;
        cout << "=============================================" << endl;
        for (auto first=firsts[0]->begin(); first!=firsts[0]->end(); first++){
            cout << *first << " ";
        }
        cout << endl;
    }

    for (int j=0; j!=r; j++){
        cout << "=============================================" << endl;
        cout << "last" << j+1 << endl;
        cout << "=============================================" << endl;
        for (auto last=lasts[0]->begin(); last!=lasts[0]->end(); last++){
            cout << *last << " ";
        }
        cout << endl;
    }

   // Merge structures necessary for the loop construction
    LINK = MergeVectors(links);
    VTX  = MergeVectors(vtxs);
    //----------------------------------------------------------------------        
    // Connect replicas 
    //----------------------------------------------------------------------        
    int kfirst = -1;
    int klast  = -1;
    for (int j=0; j!=firsts[0]->size(); j++){
        //If spin belongs to the region of interest
        if  (find(Aregion.begin(),Aregion.end(),j)!=Aregion.end()){
            //Connect replicas together
            kfirst = -1;
            klast  = -1;
            for (int k=0; k!=r; k++){
                //If the spins is being acted upon in this replica
                if  (firsts[k]->at(j)!=-1){
                    //And it is not the first replica with the spin being acted upon    
                    if  (kfirst != -1){
                        //Connect 2 inner replicas by updating LINK structure
                        LINK[lasts[klast]->at(j)] = firsts[k]->at(j);
                        LINK[firsts[k]->at(j)]    = lasts[klast]->at(j);
                    }
                    else  kfirst = k;
                    klast = k;
                }
            }//Replica loop    
        
            //Finally, connect outer replicas
            if  (klast!=-1){
                LINK[lasts[klast]->at(j)]   = firsts[kfirst]->at(j);
                LINK[firsts[kfirst]->at(j)] = lasts[klast]->at(j);
            }
        }//If: the spin belong to A 
   
        //If the spin doesnt belong to the region of interest
        else{
            //Connect each replica to itself
            for (int k=0; k!=r; k++){
                //But only, if the spin is being acted upon
                if  (firsts[k]->at(j)!=-1){
                    LINK[lasts[k]->at(j)]  = firsts[k]->at(j);
                    LINK[firsts[k]->at(j)] = lasts[k]->at(j);
                }
            }//Replica loop
        }//If: spin doesnt belong to A 
    }//Spins loop

    cout << "=============================================" << endl;
    cout << "LINK" << endl;
    cout << "=============================================" << endl;
    for (auto link=LINK.begin(); link!=LINK.end(); link++){
        cout << *link << " ";
    }
    cout << endl;
    //----------------------------------------------------------------------        
    // Construct loops 
    //----------------------------------------------------------------------        
    long j0;                      //Loop entrance leg
    long j;                       //Current leg
    long p;
    pair<long,long> legtype;      //Struct with the leg, and operator type   
    
    //Construct Nl number of loops 
    if (nTotal>0){
        for (long i=0; i!=Nloops; i++) {
            j0 = uRandInt()%(4*nTotal);       //Pick a random loop entrance leg among all possible legs
            j  = j0;
            cout << j << " b " << endl;    
            //NvisitedLegs[i] = 0;         //Number of visited legs for i'th loop
            //Construct an operator-loop
            do  {
                    p = (long) j/4;                     //Current operator index
                    legtype = SwitchLeg(j%4,VTX[p]);    //Get the next leg and the new operator's type
                    j       = legtype.first +4*p;       //Move to the next leg
                    VTX[p]  = legtype.second;           //Update the type of the operator
                    NvisitedLegs[i] += 1;
                    cout << j << "  b" << endl;
                    if   (j == j0) break;               //If the loop is closed, we are done
                    else{ 
                        j = LINK[j];                   //Else move to the next linked leg
                        NvisitedLegs[i] += 1;
                        cout << j << " b" << endl;    
                    }
            } while  (j !=j0);                          //Another way to close the loop
        }    
    }
    //----------------------------------------------------------------------        
    //  Map back the changes to the operator list
    //----------------------------------------------------------------------        
    long leg;
    long b;
    for (int k=0; k!=r; k++){
        p = -1;
        for (auto oper=sms[k]->begin(); oper!=sms[k]->end(); oper++) 
            if (*oper != 0){
                b = (long)((*oper-(*oper)%2)/2);             //Determine operator's bond
                p++;                                         //Increase the operator counter
                
               *oper = 2*b + (long)((sgn(VTX[p+shifts[k]]-5)+1)/2);    //The new operator bond is unchanged
            
            }                                                //but its type could have changed
    }
    //----------------------------------------------------------------------        
    //  Map back the changes to the spins state
    //----------------------------------------------------------------------        

    for (int k=0; k!=r; k++){
        for (long j=0; j!=firsts[k]->size(); j++){
            if  (firsts[k]->at(j)==-1){                  //If the spin is not acted by an operator P
                if  (Replicas[k]->uRand()<0.5)             //Flip it randomly
                    spins[k]->at(j) = -spins[k]->at(j); 
            }
            else{
                p   = (long) firsts[k]->at(j)/4;          //Get the first operator acting on it
                leg = (long) firsts[k]->at(j)%4;          //Get the spin's leg on this operator  
                spins[k]->at(j) = LegSpin[VTX[p]-1][leg]; //Use the leg/operator-type -> spin-type map
            }
        }//Spins loop 
     }//Replica loop  
}

//#################################################################
// There are 4 different types of moves for pair-longeracting bonds.
// Below they are implemented as seperate functions. For a given
// enLeg, each function returns the exit leg after the move is done.
// Leg labelling convention : bottom left leg is 0,  the label # is 
// increased by one in the counter-clockwise direction.     
//#################################################################


long SSEXY::Bounce(long enLeg){
    return 0;
}

long SSEXY::ContinueStraight(long enLeg){
    return enLeg - sgn(enLeg-2)*(1-(sgn((enLeg%3)-1)-1));
}

long SSEXY::SwitchReverse(long enLeg){
    return enLeg - sgn(enLeg%2-1);
}

long SSEXY::SwitchContinue(long enLeg){
    return enLeg - 2*sgn(enLeg - 2);
}
//#################################################################
//#################################################################

pair<long,long> SSEXY::SwitchLeg(long enLeg, long vtype){
    long exLeg   = -1;
    long newtype = -1;
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
return pair <long,long> (exLeg,newtype);
}
//######################################################################
// Main
//######################################################################

 
int main()
{
    SSEXY ssexy(1,6,1,0.1,5);  
    for (int i=0; i!=1000; i++){
        ssexy.MCstep();
    }
    return 0;
} 
