#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include "ssexy.h"
#include "lattice.h"
#include "lattice.cpp"
#include "replica.cpp"
//#include "communicator.h"
//#include "communicator.cpp"
#include "stdlib.h"
#include <math.h>
#include <string.h>
#include <iomanip>
#include "helper.cpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;



//**************************************************************************
SSEXY:: SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, float _Beta, 
              long seed, bool _measSS, int _Asize, string frName, 
              vector<long>* _Anor, vector<long>* _Ared, vector<long>* _Aext): 
communicator(_Nx,_Ny,_r,_T,_Beta,seed,frName,_Asize), RandomBase(seed)
{
    long tmp[6][4] = {  {-1,-1,-1,-1},
                        { 1, 1, 1, 1},
                        { 1,-1,-1, 1},
                        {-1, 1, 1,-1},
                        {-1, 1,-1, 1},
                        { 1,-1, 1,-1}
                    };
    memcpy(LegSpin,tmp,sizeof tmp);
    //Store the simulation id
    SSEXID = communicator.getId();

    // Initialize replicas
    Nx = _Nx;
    Ny = _Ny;
    N  = Nx*Ny;
    r = _r;

    //Define temperature in one of two ways
    if (_T != -1){
        T = _T;
        Beta = 1.0/(1.0*T);
    }
    else{
        Beta = _Beta;
        T = 1.0/(1.0*_Beta);
    } 
    
    //Initiate the replicas
    for(int i=0; i!=r; i++){
        Replicas.push_back(new Replica(_Nx,_Ny,T,seed, i+1));
    }
    
    //Load replicas' datastructures if needed
    if (frName!="") LoadState();
    
    Debug         = false;
    DebugSRT      = true;
    DebugILRT     = true;
    Nloops        = 1;    
    nMeas         = 0;
    if  (DebugSRT) binSize = 100;
    else           binSize = 100;
    saveFreq      = 1; 
    nSaved        = 0;
    maxLoopSize   = 2400000;
    SpinStiffness = 0;
    nAred         = 0;
    nAext         = 0;
    nAredRT       = 0;
    nAextRT       = 0;
    measSS        = _measSS;
    measRatio     = not (_Aext->empty());

    if (measSS)    cout << "Measuring spin stiffness" << endl;
    if (measRatio) cout << "Measuring Z ratio" << endl;

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
    TNvisitedLegs.resize(Nloops,0);
    Tns.resize(r+1,0);              //+1 is for the totals
    Partitions.clear(); 

    //Initialize region A
    Aregion = _Anor; 
  
    //For ratio trick there are two regions A.
    //Initialize them and their set  difference.
    if  (measRatio){
        Ared = *_Ared;   //Reduced A region
        Aext = *_Aext;   //Extended A region
        Adif = {};       //Their difference: elements contained in Aext that arent in Ared 
        if  (Aext.size()>Ared.size()){
            for (auto spin=Aext.begin(); spin!=Aext.end(); spin++){
                if  (find(Ared.begin(),Ared.end(),*spin)==Ared.end())
                    Adif.push_back(*spin);
                }
            }
        else{
            for (auto spin=Ared.begin(); spin!=Ared.end(); spin++){
                if  (find(Aext.begin(),Aext.end(),*spin)==Aext.end())
                    Adif.push_back(*spin);
                }
            }
            
        if  (Adif.empty()){
            cout << "Error: region A and its extension are equal" << endl;
            exit(0);
        }
        
        cout << "A" << endl;
        for (auto spin=Aregion->begin(); spin!=Aregion->end(); spin++){
            cout << *spin << ' ';
        }
        cout << endl;
 
        cout << "Ared" << endl;
        for (auto spin=Ared.begin(); spin!=Ared.end(); spin++){
            cout << *spin << ' ';
        }
        cout << endl;
 
        cout << "Aext" << endl;
        for (auto spin=Aext.begin(); spin!=Aext.end(); spin++){
            cout << *spin << ' ';
        }
        cout << endl;
 
        cout << "Adiff" << endl;
        for (auto spin=Adif.begin(); spin!=Adif.end(); spin++){
            cout << *spin << ' ';
        }
        cout << endl;
    }

   //Write headers for a new estimator file
    string eHeader;
    if (frName==""){
        eHeader = boost::str(boost::format("#%15s%16s%16s")%"nT"%"ET"%"Legs");
        if (measSS)    eHeader += boost::str(boost::format("%16s")%"SS");
        if  (DebugSRT){
            if (measRatio) eHeader += boost::str(boost::format("%16s")%"nAred");
            if (measRatio) eHeader += boost::str(boost::format("%16s")%"nAext");
        }
        if  (DebugILRT){
            if (measRatio) eHeader += boost::str(boost::format("%16s")%"LRatio");
            if (measRatio) eHeader += boost::str(boost::format("%16s")%"nAredRT");
            if (measRatio) eHeader += boost::str(boost::format("%16s")%"nAextRT");
        }
        *communicator.stream("estimator") << eHeader;    


        if  (r>1)
            for (int j=0; j!=r; j++){
                eHeader = boost::str(boost::format("%15s%i%15s%i") %"n"%(j+1)%"E"%(j+1));
                *communicator.stream("estimator") << eHeader;    
            }
        *communicator.stream("estimator") << endl;    
    }
} 
        
int SSEXY::AdjustParameters()
{
    long Tn;    //Cummulative n
    long TLegs; //Cummulative Legs
    Tn    = 0;
    TLegs = 0;

    //Run simulation without recording results
    for (int i=0; i!=100; i++){
        //Perform a Monte-Carlo step
        MCstep();
        //Accumulate M's
        for (int j=0; j!=r; j++)
            Tn += Replicas[j]->getn();
        //Accumulate legs 
        for (int j=0; j!=Nloops; j++)
            TLegs += NvisitedLegs[j];
    } 
    
    //Increase Nloops if not enough legs
    // are being generated
    if  (TLegs > 2*Tn)
        return 0;
    else{
        Nloops += 1;
        NvisitedLegs.resize(Nloops,0);
        TNvisitedLegs.resize(Nloops,0);
        cout << SSEXID << ": Legs = " << TLegs/100.0 << " 2xnT = " << Tn*2/100.0 << endl;
        cout << SSEXID << ": Increase # of loops to " << Nloops << endl;
        return 1;
    }
}



/**************************************************************************
*Get the next spin as determined by boundary condition (connected or not)
***************************************************************************/
int SSEXY::BCnextSpin(int sindex,int& replica, bool connected){

    if (sindex<N) sindex +=N;
    else          sindex -=N;
    if  (connected)
        replica = !replica;

    return sindex;
}    

/**************************************************************************
*Measure the number of distinct loops existing between 2 partitions when
*their boundary condition are taken into account. 
***************************************************************************/
long SSEXY::LoopPartition(vector<long>& BC){

    Partitions.clear(); 
    Partitions.push_back(*(Replicas[0]->getPart()));
    Partitions.push_back(*(Replicas[1]->getPart()));
   
    int i;
    if  (Debug){
        cout << endl << "Before connection=================================" << endl;
        cout << "Boundary conditions: " << endl;
        for (auto ibc=BC.begin(); ibc!=BC.end(); ibc++)
            cout << *ibc << " ";
        cout << endl;
        for (int ir=0; ir!=2; ir++){
            i=0;
            for (auto spin=Partitions[ir].begin(); spin!=Partitions[ir].end(); spin++){
                cout << setw(4) << *spin;
                i += 1;
                if  (i == N) cout << endl; 
            }
        cout << endl << endl;
        }
    }
    int  replica;   
    int  spin;
    int  nspin;
    bool connected;
    long nLoop = 0;
    int  ospin;
    int  oreplica;
    //Repeat for all edge spins in both replicas
    DeterPaths.clear();

    for (auto oreplica=0; oreplica!=2; oreplica++){
        for (auto ispin=0; ispin!=2*N; ispin++){
            
            //If the spin hasnt been visited
            if  (not(Partitions[oreplica][ispin]<0)){
                //Follow the BC loop until it comes back to the initial spin
                nLoop += 1;
                ospin   = ispin;
                spin    = ospin;
                replica = oreplica;
                DeterPaths[nLoop];
                if (Debug) cout << "(r,s) = (" << replica << "," << spin << ")" << endl;
                do{
                    //Switch to the other end of the loop the spin belongs to
                    nspin = Partitions[replica][spin];
                    if (Debug) cout << "L: (r,s) = (" << replica << "," << nspin << ")" << endl;

                    //Mark the visited spins by the negative 
                    //of the loop index that they belong to
                    Partitions[replica][spin]  = -nLoop;
                    Partitions[replica][nspin] = -nLoop;
                    
                    //Switch to the spin connected by BC
                    connected = not(find(BC.begin(),BC.end(),nspin%N)==BC.end());
                    spin = BCnextSpin(nspin, replica, connected);       
         
                    //Add vertices path segment to the loop 
                    //Necessary for a deterministic loop update
                    list<long> tpath = Replicas[replica]->getLoopPaths()->at(spin); 
                    if  (replica == 1)
                        for (auto aleg=tpath.begin(); aleg!=tpath.end(); aleg++)
                            DeterPaths[nLoop].push_back(*aleg+4*ns[0]);
                    else
                        DeterPaths[nLoop].insert(DeterPaths[nLoop].end(),tpath.begin(),tpath.end());
                    

                    if (Debug) cout << "B: (r,s) = (" << replica << "," << spin << ")" << endl;


                } while ((spin!=ospin) or (replica!=oreplica));
                if  (Debug){
                     for (int ir=0; ir!=2; ir++){
                         i=0;
                         for (auto spin=Partitions[ir].begin(); spin!=Partitions[ir].end(); spin++){
                             cout << setw(4) << *spin;
                             i += 1;
                             if  (i == N) cout << endl; 
                         }
                     cout << endl;
                     }
                     cout << "Deterministic path " << nLoop << " : " ;
                     for (auto iloop=DeterPaths[nLoop].begin(); iloop!=DeterPaths[nLoop].end(); iloop++)
                         cout << *iloop << " ";
                     cout << endl << endl;
                 }
            }
        }
    }
  
    //cout << "#Loops: " << nLoop << endl;
    return nLoop; 
}        




/**************************************************************************
* Deterministic loops off-diagonal update.
* It builds loops based on deterministic vertex moves and then flips each
* one of them with 50% probability. The end product is an updated VTX vector.
***************************************************************************/
long SSEXY::DeterministicOffDiagonalMove(){
    long nLoops  = LoopPartition(*Aregion);
    long enleg;
    long exleg;
    long p;
    pair<long,long> legvtx;
    bool DetDebug = false;   
 
    if  (DetDebug){
        cout << "Offdiagonal deterministic update" << endl;
        cout << "Deterministic loops: " << nLoops << endl;
    }
        
    for (auto iLoop=1; iLoop!=(nLoops+1); iLoop++){
        if  (uRand() < 0.5){ 
            if  (DetDebug)
                cout << "Loop: " << iLoop << endl;
            long i = 0;
            for (auto leg=DeterPaths[iLoop].begin(); leg!=DeterPaths[iLoop].end(); ++leg){
                //cout << "i: " << ++i << endl;
                if (*leg!=exleg){
                    enleg = *leg;
                    exleg = *(++leg);
                    p     = (long) (*leg)/4;
                    if  (DetDebug)
                        cout << enleg << " " << exleg << " " ;
                    
                    //Deterministic out the new vertex type by hacking SwitchLeg method
                    legvtx = SwitchLeg(enleg%4,VTX[p],-1);
                    if  ((legvtx.first-(exleg%4)) != 0)
                        legvtx = SwitchLeg(enleg%4,VTX[p],2);
                    
                    //Flip the vertex type
                    VTX[p] = legvtx.second;
                }
            }
            if  (DetDebug)
                cout << endl;
        }
    }
}





/**************************************************************************
* Random loops off-diagonal update.
* It generates a worm that moves randomly through the legs-space until 
* it closes on itself by reaching its own tail. Every leg on worm's path
* gets flipped. The end product is an updated VTX vector.
***************************************************************************/
long SSEXY::RandomOffDiagonalUpdate(){

    // Connect replicas 
    int kfirst = -1;
    int klast  = -1;
    for (int j=0; j!=firsts[0]->size(); j++){
        //If spin belongs to the region of interest
        if  (find(Aregion->begin(),Aregion->end(),j)!=Aregion->end()){
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
            NvisitedLegs[i] = 0;         //Number of visited legs for i'th loop
            //Construct an operator-loop
            do  {
                    p = (long) j/4;                     //Current operator index
                    legtype = SwitchLeg(j%4,VTX[p],0.5);//Get the next leg and the new operator's type
                    j       = legtype.first +4*p;       //Move to the next leg
                    VTX[p]  = legtype.second;           //Update the type of the operator
                    NvisitedLegs[i] += 1;
                    if   (j == j0) break;               //If the loop is closed, we are done
                    else{ 
                        j = LINK[j];                   //Else move to the next linked leg
                        NvisitedLegs[i] += 1;
                    }
                    
                    //Limit the maximum loop size
                    if  (NvisitedLegs[i]>maxLoopSize){
                        cout << SSEXID << ": Extremely large loop" << endl;
                        return 1;
                    }          
                         
            } while  (j !=j0);                          //Another way to close the loop
        }    
    }

}


/**************************************************************************
* Measure Z[Aregion]/Z[Aextended], i.e. the ratio of modified part-
*tion function connected at Aregion and the one connected at Aextended. 
***************************************************************************/
float SSEXY::ALRTrick(){
    //Measure the number of loops formed in each
    //modified geometry.
    //When we call MeasureNLoop for the 2nd time,
    //all necessary structures are already computed.
    //That is why we need a special boolean flag in 
    //order to reduce the computational effort.
    long AnLoops  = LoopPartition(Ared);
    long EAnLoops = LoopPartition(Aext);
//    cout << "Ratio: " << (1.0*EAnLoops)/(1.0*AnLoops) << endl;
    return pow(2,AnLoops-EAnLoops);
}


/**************************************************************************
* Attemp to switch the size of region A based on the boundary conditions
***************************************************************************/
float SSEXY::ILRTrick(){

    //Partition edge spins according to the loops they belong to
    LoopPartition(Ared);

    bool lDebug = false;
    if  (lDebug){
        cout << endl << "Loops at deltaA=================================" << endl;
        for (int r=0; r!=2; r++){
            for (int bs=0; bs!=2; bs++){
                for (auto spin=Adif.begin(); spin!=Adif.end(); spin++)
                    cout << setw(4) << Partitions[r][*spin+N*bs];
                cout << endl;
            } 
        cout << endl << endl;
        }
    }


    long downLoop;
    long upLoop;
    vector<long> twoLoops (2,0);
    map<long, set<long>> ConnectedLoops;
    set<long> AllLoops;
    //Merge loops that share common spins from Adif
    for (auto ADspin=Adif.begin(); ADspin!=Adif.end(); ADspin++){
        //If region A is increased, BCs exist between replicas 
        if  (Ared.size()<Aext.size()){    
            twoLoops[0] = Partitions[0][*ADspin]; 
            twoLoops[1] = Partitions[1][*ADspin+N]; 
        }
        //If region A is decreased, BCs exist within replicas
        else{
            twoLoops[0] = Partitions[0][*ADspin]; 
            twoLoops[1] = Partitions[0][*ADspin+N]; 
        }

        //Create a map of loop connections. 
        //Key = loop index, value = set of connected loops
        for (int i=0; i!=2; i++){
            downLoop = twoLoops[i];
            upLoop   = twoLoops[(i+1)%2];
            if  (not ConnectedLoops.count(downLoop)){
                set<long> lvalue {upLoop};  
                ConnectedLoops[downLoop] =lvalue; 
                AllLoops.insert(downLoop);
            }
            else
                ConnectedLoops[downLoop].insert(upLoop);
            }
    }

    if  (lDebug){
        cout << "Connected loops map: " << endl;
        for (auto loop=ConnectedLoops.begin(); loop!=ConnectedLoops.end(); loop++){
            cout << loop->first << ": ";
            for (auto loop2=loop->second.begin(); loop2!=loop->second.end(); loop2++){
                cout << *loop2 << " ";
            }
            cout << endl;
        }
    }

    long L0 = AllLoops.size();
    long Ls = 0;

    set<long> MarkedLoops; 
    set<long> path;

    if  (lDebug) cout << "Paths: " << endl;

    for (auto loop=AllLoops.begin(); loop!=AllLoops.end(); loop++)
        if  (not MarkedLoops.count(*loop)){
            GetConnectedSubraph(ConnectedLoops, path, *loop);
            for (auto mloop=path.begin(); mloop!=path.end(); mloop++){
                MarkedLoops.insert(*mloop);
                if  (lDebug) cout << *mloop << " ";
            }
            if  (lDebug) cout << endl;

            path.clear();  
            Ls += 1;
        }          

    if  (lDebug) cout << "Ls = " << Ls << " L0 = " << L0 << endl;
     
    nAredRT += pow(2,L0);
    nAextRT += pow(2,Ls);

    return pow(2,Ls-L0);
            
}


/**************************************************************************
* Attemp to switch the size of region A based on the boundary conditions
***************************************************************************/
vector<long>* SSEXY::SwitchAregion()
{
    bool Switch = true;    
    int  S0;
    int  S1;
    //if  (Aregion == &Ared){
    if (Ared.size() < Aext.size()){   
        for (auto sindex=Adif.begin(); sindex!=Adif.end(); sindex++){
            S0 = Replicas[1]->getSpin()->at(*sindex);
            S1 = Replicas[0]->getTedge()->at(*sindex);
            if  (S0!=S1){
                Switch = false;
                break;
            }
        }
        if  (Switch)
            //Aregion = &Aext;
            nAext +=1;
        nAred +=1;
        }
    else{
        for (auto sindex=Adif.begin(); sindex!=Adif.end(); sindex++){
            S0 = Replicas[0]->getSpin()->at(*sindex);
            S1 = Replicas[0]->getTedge()->at(*sindex);
            if  (S0!=S1){
                Switch = false;
                break;
            }
        } 
        if  (Switch)
            //Aregion = &Ared;
            nAred += 1;
        nAext += 1;
    }        
    return Aregion; 
}
        

/**************************************************************************
* Perform measuements
**************************************************************************/

int SSEXY::Measure()
{
    //Accumulate the new measurement
    nMeas += 1;
    for (int j=0; j!=r; j++)
        Tns[j] += ns[j];
    Tns[r] += nTotal;

    //Accumulate the legs 
    for (int j=0; j!=Nloops; j++)
        TNvisitedLegs[j] += NvisitedLegs[j];

    //Accumulate the spin stiffness estimator
    if (measSS)
       SpinStiffness += Replicas[0]->MeasureSpinStiffness(); 
    

    //Accumulate the partition function ratio estimator
    if  (measRatio)
        if  (DebugSRT){
            //if  (Aregion == &Ared) nAred += 1;
            //else                   nAext += 1; 
        }
        if  (DebugILRT)
            LRatio += ILRTrick();
    
    // If we've collected enough of measurements
    float E;
    long  TNLegs = 0;
    if  (nMeas == binSize){
        //Record total Energy
        E = -((float) Tns[r]/((float) binSize*N))/Beta + r*1.0;   //r*1.0 term represents the added energy offset per bond
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E%16.8E") %(Tns[r]/(1.0*binSize)) %E);
        
        //Record loop data
        for (int j=0; j!=Nloops; j++)
            TNLegs += TNvisitedLegs[j];
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(TNLegs*1.0/(1.0*binSize)));

        //Record spin stiffness if needed
        if (measSS){
           *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(SpinStiffness/(1.0*binSize)));
           SpinStiffness = 0;
        }

        //Record partition function ratio if needed
        if  (measRatio){
            if  (DebugSRT){
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(1.0*nAred/(1.0*binSize)));
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(1.0*nAext/(1.0*binSize)));
                nAred = 0;
                nAext = 0;
            }
            if  (DebugILRT){
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(1.0*LRatio/(1.0*binSize)));
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(1.0*nAredRT/(1.0*binSize)));
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E") %(1.0*nAextRT/(1.0*binSize)));
                LRatio = 0;
                nAredRT = 0;
                nAextRT = 0;
            }
        }

       //Record energy of each replica
        if  (r>1)         
        
            for (int j=0; j!=r; j++){
                E = -((float) Tns[j]/((float) binSize*N))/Beta + 1.0;
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E%16.8E") %(Tns[j]/(1.0*binSize)) %E);
            }

        //Carry to a new line
        *communicator.stream("estimator") << endl;    

        //Set accumulating variables to 0
        for (int j=0; j!=(r+1); j++)
            Tns[j] = 0;
        nMeas  = 0;
        
        for (int j=0; j!=Nloops; j++)
            TNvisitedLegs[j]=0;
        
        //Increase the counter of recordered bins
        nSaved += 1;
        if  (nSaved == saveFreq){
            SaveState();
            nSaved = 0;
        
        //Notify about a new measurement taken
        cout << SSEXID << ": Measurement taken" << endl;
        }
    }

    return 0;
}

int SSEXY::SaveState(){
    //Erase previous state info
    communicator.reset("state");

    //Save the state of the random generator
    uRandInt.distribution().reset(); 
    uRand.distribution().reset(); 
    *communicator.stream("state") << eng << endl;
    
    //Save each replica's state
    for (int j=0; j!=r; j++){
        //Save random-engine internal state
        Replicas[j]->uRandInt.distribution().reset(); 
        Replicas[j]->uRand.distribution().reset(); 
        *communicator.stream("state") << Replicas[j]->eng;
        *communicator.stream("state") << endl;

        //Save operators
        for (auto oper = Replicas[j]->getOper()->begin(); oper!=Replicas[j]->getOper()->end(); oper++){
            *communicator.stream("state") << *oper << " ";
        }
        *communicator.stream("state") << endl;

        //Save spins
        for (auto spin= Replicas[j]->getSpin()->begin(); spin!=Replicas[j]->getSpin()->end(); spin++){
            *communicator.stream("state") << *spin << " ";
        }
        *communicator.stream("state") << endl;
    }   
}

int SSEXY::LoadState(){
    int n;
    int i;
    long N;    

    //First load ssexy generator state

    //Buffer variables
    string       sBuf0;
    stringstream ssBuf0;    
    
    //Load random-engine internal state
    uRandInt.distribution().reset(); 
    uRand.distribution().reset(); 
    getline(*communicator.stream("state"),sBuf0);     //1st line of the state file 
    ssBuf0 << sBuf0;                                   
    ssBuf0 >> eng;

    //Save each replica's state
    for (int j=0; j!=r; j++){
        //Buffer variables
        string       sBuf;
        stringstream ssBuf;    
        //cout << "Replica #" << j+1 << endl;
        //Load random-engine internal state
        Replicas[j]->uRandInt.distribution().reset(); 
        Replicas[j]->uRand.distribution().reset(); 
        getline(*communicator.stream("state"),sBuf);     //1st line of the state file 
        ssBuf << sBuf;                                   
        ssBuf >> Replicas[j]->eng;

        //Load operators
        istringstream issBuf1;
        getline(*communicator.stream("state"),sBuf);     //2nd line of the state file
        issBuf1.str(sBuf);
        Replicas[j]->getOper()->resize(0,0);
        N = 0;                                           //Non-identity operator counter
        while (issBuf1 >> n){
              if (n!=0) N += 1;
              Replicas[j]->getOper()->push_back(n);
        }
        Replicas[j]->setn(N);
        Replicas[j]->setM(Replicas[j]->getOper()->size());
        ////cout << "n = " << N << " M = " << Replicas[j]->getOper()->size() << endl;
        //cout << "Operators line: " << sBuf << endl;
        //cout << "Operators:" << endl;
        //for (auto oper=Replicas[j]->getOper()->begin(); oper!=Replicas[j]->getOper()->end(); oper++)
        //    cout << *oper << " ";
        //cout << endl;

        //Load spins
        istringstream issBuf2;
        getline(*communicator.stream("state"),sBuf);     //3rd line of the state file
        issBuf2.str(sBuf);
        i=0;
        //cout << "----- spins------" << endl;
        while (issBuf2 >> n){
            Replicas[j]->getSpin()->at(i)=n;
            i += 1;
        }
        //cout << "Spins line: " << sBuf << endl;
        //cout << "Spins:" << endl;
    }
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
        if  (Replicas[j]->DiagonalMove()==1)
            Replicas[j]->AdjustM();
        Replicas[j]->ConstructLinks();
        //if  (measRatio and DebugILRT)
        //    Replicas[j]->GetDeterministicLinks();
        Replicas[j]->GetDeterministicLinks();
        firsts[j] = Replicas[j]->getFirst();        
        lasts[j]  = Replicas[j]->getLast();        
        links[j]  = Replicas[j]->getLink();        
        spins[j]  = Replicas[j]->getSpin();      
        vtxs[j]   = Replicas[j]->getVtx();      
        sms[j]    = Replicas[j]->getOper();      
        ns[j]     = Replicas[j]->getn();
        if (j>0)
            shifts[j] = shifts[j-1] + ns[j-1];
        nTotal   += ns[j];    
    }

    if  (measRatio and DebugSRT){
        Aregion = SwitchAregion();
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
    
   // Merge structures necessary for the loop construction
    LINK = MergeVectors(links);
    VTX  = MergeVectors(vtxs);

    //RandomOffDiagonalUpdate();
    //DeterministicOffDiagonalMove();
    RandomOffDiagonalUpdate();

    //----------------------------------------------------------------------        
    //  Map back the changes to the operator list
    //----------------------------------------------------------------------        
    long leg;
    long b;
    long p;
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
    bool previous;
    bool flipRand;
    for (long j=0; j!=firsts[0]->size(); j++){
        if  (find(Aregion->begin(),Aregion->end(),j)!=Aregion->end()){
            previous = false;
            flipRand = true;
            for (int k=0; k!=r; k++){
                //If the spin is inactive in the current replica
                if  (firsts[k]->at(j)==-1){                        //If the spin is not acted upon
                    //But if it was active in the previous replica, the
                    //current replica spin must get changed according 
                    //to its last leg type in the preivous replica
                    if  (previous == true){
                        p   = (long) lasts[k-1]->at(j)/4;          //Get the last operator acting on it
                        leg = (long) lasts[k-1]->at(j)%4;          //Get the spin's leg on this operator  
                        spins[k]->at(j) = LegSpin[VTX[p]-1][leg];  //Use the leg/operator-type -> spin-type map
                    }
                }
                //If the spin is active in the current replica
                else{
                    //If the spin was active in the previous replica,
                    //just update the current replica spin type.
                    if  (previous == true){
                        p   = (long) firsts[k]->at(j)/4;          //Get the first operator acting on it
                        leg = (long) firsts[k]->at(j)%4;          //Get the spin's leg on this operator  
                        spins[k]->at(j) = LegSpin[VTX[p]-1][leg]; //Use the leg/operator-type -> spin-type map
                    }
                    //Otherwise, synchronize the spin type in all previous
                    //replicas having this spin inactive with the current 
                    //replica spin type.
                    else{
                        //Determine the new spin type
                        p   = (long) firsts[k]->at(j)/4;          //Get the first operator acting on it
                        leg = (long) firsts[k]->at(j)%4;          //Get the spin's leg on this operator  
                        spins[k]->at(j) = LegSpin[VTX[p]-1][leg]; //Use the leg/operator-type -> spin-type map
                        //Synchronize it with previous replicas
                        for (int t=k; t!=0; t--)
                            if   (firsts[t-1]->at(j) != -1) break;
                            else spins[t-1]->at(j) = spins[k]->at(j);
                    } 
                    flipRand = false;
                    previous = true;
                }//else if the spin is active     
            }//replica loop

            //if the spin is inactive in all replicas
            //flip it randomly 
            if  (flipRand == true)
                if  (Replicas[0]->uRand()<0.5)
                    for (int k=0; k!=r; k++)
                        spins[k]->at(j) = -spins[k]->at(j); 
        }//if the spin doesnt belong to region A
        else    
            for (int k=0; k!=r; k++){
                if  (firsts[k]->at(j)==-1){                    //If the spin is not acted upon
                    if  (Replicas[k]->uRand()<0.5)             //Flip it randomly
                        spins[k]->at(j) = -spins[k]->at(j); 
                }
                else{
                    p   = (long) firsts[k]->at(j)/4;          //Get the first operator acting on it
                    leg = (long) firsts[k]->at(j)%4;          //Get the spin's leg on this operator  
                    spins[k]->at(j) = LegSpin[VTX[p]-1][leg]; //Use the leg/operator-type -> spin-type map
                }
            }
     }//Spins loop
     return 0; 
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

pair<long,long> SSEXY::SwitchLeg(long enLeg, long vtype, float prob){
    long exLeg   = -1;
    long newtype = -1;
    switch (vtype){
        case 1:
                    if  (uRand()<prob){
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
                    if  (uRand()<prob){
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
                    if  (uRand()<prob){
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
                    if  (uRand()<prob){
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
                    if  (uRand()<prob){
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
                    if  (uRand()<prob){
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

 
int main(int argc, char *argv[])
{
    po::options_description cmdLineOptions("Command line options"); 
    po::options_description simulationOptions("Simulation options"); 
    po::options_description measurementOptions("Measurements options"); 
    po::options_description physicalOptions("Physical parameters options"); 
    po::variables_map params;
    simulationOptions.add_options()
            ("help,h", "produce help message")
            ("state,s",      po::value<string>()->default_value(""),"path to the state file")
            ("process_id,p", po::value<int>()->default_value(0),"process id")
            ("replica,r",    po::value<int>()->default_value(1),"number of replicas")
            ("region_A,a",   po::value<string>()->default_value(""),"path to the file defining region A.\n"
                             "If set to an integer value, "
                             "it defines the number of consecutif spins in region A. ")
            ;
    measurementOptions.add_options()
            ("temperature,T",po::value<double>()->default_value(-1), "temperature")
            ("beta,b",       po::value<double>()->default_value(-1),"inverse temperature")
            ("width,x",      po::value<int>(),"lattice width")
            ("height,y",     po::value<int>()->default_value(1),"lattice height")
            ;
    physicalOptions.add_options()
            ("measn,m",      po::value<long>(),"number of measurements to take")
            ("super,w",      "turn on the spin stifness measurement. \n(r must be set to 1)")
            ("rtrick,t",     po::value<string>()->default_value(""),"path to the file defining extended region A.")
            ("region_Ared",  po::value<string>()->default_value(""),"path to the file defining reduced region A. ")
            ;
    cmdLineOptions.add(simulationOptions).add(measurementOptions).add(physicalOptions);
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), params);
    po::notify(params);

    if (params.count("help")) {
       cout << cmdLineOptions << "\n";
       return 1;
    } 
   
    if  ((params.count("super")) && (params["replica"].as<int>() != 1)){
        cerr << "Error: cannot measure spin stiffness for a multiple replicas simulation" << endl;
        return 1; 
    }
    
    if  ((params["temperature"].as<double>() != -1) && (params["beta"].as<double>() != -1)){
        cerr << "Error: simultanious definition of temperature via T and beta parameters" << endl;
        return 1; 
    }
    
    if  (!(params.count("measn"))){
        cerr << "Error: define the number of measurements to take" << endl;
        return 1; 
    }
    
    if  (!(params.count("width"))){
        cerr << "Error: define lattice width" << endl;
        return 1; 
    }
    
    if  ((params["temperature"].as<double>() == -1) && (params["beta"].as<double>() == -1)){
        cerr << "Error: define temperature via beta or T paramater" << endl;
        return 1; 
    }

    // Attempt to define region A and its extensions--------------------------------
    if  ((params["region_Ared"].as<string>()!="") and (params["rtrick"].as<string>()=="")){
        cout << "Error: reduced A region can be defined only when --rtrick flag is set" << endl;
        return 1;
    }
    
    
    string region_Ared;
    region_Ared = params["region_Ared"].as<string>();
    if  (region_Ared=="")
         region_Ared = params["region_A"].as<string>();
    //if  (params.count("region_A") and params.count("rtrick")){
        //char* endAred
        //char* endAext;
        //int  sizeAred = strtol(params["region_A"].as<string>().c_str(), &endAred, 10);
        //int  sizeAext = strtol(params["rtrick"].as<string>().c_str(),   &endAext, 10);
        //if  ((!*endAred) xor (!*endAext)){
        //    cout << "Error: when rtrick and region A are both specified, they must be specified via the same method (number or file)" << endl;
        //    exit(0);
        //}
    //} 
    LATTICE Anor("A",          params["region_A"].as<string>().c_str());
    LATTICE Aext("A extended", params["rtrick"].as<string>().c_str());
    LATTICE Ared("A reduced",  region_Ared.c_str());
  
    //------------------------------------------------------------------------------


    SSEXY ssexy(params["replica"].as<int>(), params["width"].as<int>(),
                params["height"].as<int>(),  params["temperature"].as<double>(), 
                params["beta"].as<double>(), params["process_id"].as<int>(), 
                params.count("super"),       Anor.getSize(),  
                params["state"].as<string>(), 
                Anor.getLattice(), Ared.getLattice(), Aext.getLattice());  

    cout << endl << "Equilibration stage" << endl << endl;
    
    int NoAdjust;
    NoAdjust=0;
    for (int i=0; i!=0; i++){
        if  (ssexy.AdjustParameters() == 0) NoAdjust += 1;
        else NoAdjust = 0;
        
        if (NoAdjust == 10) break;
    }
        
    cout << endl << "Measurement stage" << endl << endl;
    for (long i=0; i!=100*params["measn"].as<long>(); i++){
        ssexy.MCstep();
        ssexy.Measure();
    }
    return 0;
} 
