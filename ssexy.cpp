#include <iostream>
#include <sstream>
#include "ssexy.h"
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
SSEXY::SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, float _Beta, long seed, vector<long>* _Aregion):
communicator(_Nx,_Ny,_T,seed), RandomBase(seed)
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
    for(int i=0; i!=r; i++){
        Replicas.push_back(new Replica(_Nx,_Ny,_T,seed));
    }
    Debug   = false;
    Nloops  = 1;    
    nMeas   = 0;
    binSize = 100;

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
    Tns.resize(r+1,0);              //+1 is for the totals
    
    //Debugging data structues

    sites = Replicas[0]->getSites();

    // Initialize region A
    Aregion = *_Aregion; 
    
    //Write headers
    string eHeader = boost::str(boost::format("#%15s%16s%16s%16s%16s%16s")%"nT"%"dnT"%"ET"%"dET"%"Legs"%"dLegs");
    *communicator.stream("estimator") << eHeader;    
    if  (r>1)
        for (int j=0; j!=r; j++){
            string eHeader = boost::str(boost::format("%15s%i%15s%i%15s%i%15s%i") %"n"%(j+1)%"dn"%(j+1)%"E"%(j+1)%"dE"%(j+1));
            *communicator.stream("estimator") << eHeader;    
        }
    *communicator.stream("estimator") << endl;    
} 
        

//**************************************************************************
int SSEXY::Measure()
{
    float E;
    // If we've collected enough of measurements
    if  (nMeas == binSize){
        //Record total Energy
        E = -((float) Tns[r]/((float) binSize*N))/Beta + r*1.0;   //r*1.0 term represents the added energy offset per bond
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E%16.8E%16.8E%16.8E") %(Tns[r]/(1.0*binSize)) %0.0 %E %0.0);
        
        //Record loop data
        *communicator.stream("estimator") << boost::str(boost::format("%16.8E%16.8E") %(0.0/(1.0*binSize)) %0.0);

        //Record energy of each replica
        if  (r>1)         
        
            for (int j=0; j!=r; j++){
                E = -((float) Tns[j]/((float) binSize*N))/Beta + 1.0;
                *communicator.stream("estimator") << boost::str(boost::format("%16.8E%16.8E%16.8E%16.8E") %(Tns[j]/(1.0*binSize)) %0.0 %E %0.0);
            }
        *communicator.stream("estimator") << endl;    

        //Set accumulating variables to 0
        for (int j=0; j!=(r+1); j++)
            Tns[j] = 0;
        nMeas = 0;
    }

    //Otherwise accumulate the new measurement
    else{
        for (int j=0; j!=r; j++)
            Tns[j] += ns[j];
        Tns[r] += nTotal;
        nMeas += 1;
    } 
    return 0;
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
        firsts[j] = Replicas[j]->getFirst();        
        lasts[j]  = Replicas[j]->getLast();        
        links[j]  = Replicas[j]->getLink();        
        spins[j]  = Replicas[j]->getSpin();      
        vtxs[j]   = Replicas[j]->getVtx();      
        sms[j]    = Replicas[j]->getOper();      
        ns[j]     = Replicas[j]->getn();
        if (j>0)
            shifts[j] = shifts[j-1] + ns[j-1];
        //cout << "j = " << j << " n = " << links[j]->size()/4.0 << endl;
        nTotal   += ns[j];    
    } 
    
//    ap = *spins[0];                             //propagated spins state                 
//    cout << endl << "==========================" <<endl;
//    cout << "Initial State " << endl;
//    for (auto nspin = ap.begin(); nspin!=ap.end(); nspin++)
//               cout << *nspin << " ";
//    cout << endl; 
//    long bond = 0;
//    for (vector<long>::iterator oper=sms[0]->begin(); oper!=sms[0]->end(); oper++) {
//
//        if (*oper != 0)                        //Ignore unit operator
//            bond = (long)((*oper-(*oper)%2)/2);    //Determine operator's bond based on operator's value
//        
//        if  (*oper%2 == 1){                             //If it is an off-diagonal operator
//            cout << "Operator: " << *oper << endl;
//            ap[sites[bond][0]] = -ap[sites[bond][0]];         //Update the propagated spins state
//            ap[sites[bond][1]] = -ap[sites[bond][1]]; 
//        }
//    }
//    
//    cout << endl << "Propagated State " << endl;
//    for (auto nspin = ap.begin(); nspin!=ap.end(); nspin++)
//               cout << *nspin << " ";
//
//    cout << endl; 
//    bool alert = false; 
//    for (auto aspin=Aregion.begin(); aspin!=Aregion.end(); aspin++) 
//        if  (spins[1]->at(*aspin) != ap[*aspin])
//            alert = true;
//            //cout << "A region spin defect" << endl;  
//
//    if (alert){
//       cout << "-----------------------------Alert-----------------------------" << endl;
//       for (auto nspin = spins[1]->begin(); nspin!=spins[1]->end(); nspin++)
//                cout << *nspin << " ";
//       cout << endl;
//       for (auto nspin = ap.begin(); nspin!=ap.end(); nspin++)
//                cout << *nspin << " ";
//       cout << endl;
//       }
//    else
//        cout << endl << endl << "No Alert" << endl;
    
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
    
//    for (int j=0; j!=r; j++){
//        cout << "=============================================" << endl;
//        cout << "link"<< j << endl;
//        cout << "=============================================" << endl;
//        for (auto link=links[j]->begin(); link!=links[j]->end(); link++){
//            cout << *link << " ";
//        }
//        cout << endl;
//    }
//
//    for (int j=0; j!=r; j++){
//        cout << "=============================================" << endl;
//        cout << "first"<< j << endl;
//        cout << "=============================================" << endl;
//        for (auto first=firsts[j]->begin(); first!=firsts[j]->end(); first++){
//            cout << *first << " ";
//        }
//        cout << endl;
//    }
//
//    for (int j=0; j!=r; j++){
//        cout << "=============================================" << endl;
//        cout << "last" << j+1 << endl;
//        cout << "=============================================" << endl;
//        for (auto last=lasts[j]->begin(); last!=lasts[j]->end(); last++){
//            cout << *last << " ";
//        }
//        cout << endl;
//    }

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

//    cout << "=============================================" << endl;
//    cout << "LINK" << endl;
//    cout << "=============================================" << endl;
//    for (auto link=LINK.begin(); link!=LINK.end(); link++){
//        cout << *link << " ";
//    }
//    cout << endl;
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
//            cout << j << " b " << endl;    
            //NvisitedLegs[i] = 0;         //Number of visited legs for i'th loop
            //Construct an operator-loop
            do  {
                    p = (long) j/4;                     //Current operator index
                    legtype = SwitchLeg(j%4,VTX[p]);    //Get the next leg and the new operator's type
                    j       = legtype.first +4*p;       //Move to the next leg
                    VTX[p]  = legtype.second;           //Update the type of the operator
                    NvisitedLegs[i] += 1;
//                  cout << j << "  b" << endl;
                    if   (j == j0) break;               //If the loop is closed, we are done
                    else{ 
                        j = LINK[j];                   //Else move to the next linked leg
                        NvisitedLegs[i] += 1;
//                        cout << j << " b" << endl;    
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
    bool previous;
    bool flipRand;
    for (long j=0; j!=firsts[0]->size(); j++){
        if  (find(Aregion.begin(),Aregion.end(),j)!=Aregion.end()){
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
        }//if the spin belongs to region A
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

 
int main(int argc, char *argv[])
{
    po::options_description cmdLineOptions("Command line options"); 
    po::variables_map params;
    cmdLineOptions.add_options()
            ("help,h", "produce help message")
            ("temperature,T",po::value<double>()->default_value(-1), "temperature")
            ("beta,b",       po::value<double>()->default_value(-1),"inverse temperature")
            ("width,x",      po::value<int>(),"lattice width")
            ("height,y",     po::value<int>()->default_value(1),"lattice height")
            ("process_id,p", po::value<int>()->default_value(0),"process id")
            ("replica,r",    po::value<int>()->default_value(1),"number of replicas")
            ("state,s",      po::value<string>()->default_value(""),"path to the state file")
            ("region_A,a",   po::value<string>()->default_value(""),"path to the file defining region A");
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), params);
    po::notify(params);

    if (params.count("help")) {
       cout << cmdLineOptions << "\n";
       return 1;
    } 
    
    if  ((params["temperature"].as<double>() != -1) && (params["beta"].as<double>() != -1)){
        cerr << "Error: simultanious definition of temperature via T and beta parameters" << endl;
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

    vector<long> Aregion ={};
    if  (params["region_A"].as<string>()=="" ){
        cout << "Taking A region to be empty" << endl;
    }
    else{
        ifstream RAfile (params["region_A"].as<string>());
        string line;
        string lline;
        if  (RAfile.is_open()){
            while (getline(RAfile,line))
                  lline=line;
            RAfile.close();
            istringstream sline(lline.erase(0,1));
            int n;
            while (sline >> n)
                Aregion.push_back(n);
        }
        else{
            cout << "Unable to process region A file" << endl;
            return 1;
            }    
        }


    SSEXY ssexy(params["replica"].as<int>(), params["width"].as<int>(),
                params["height"].as<int>(),  params["temperature"].as<double>(), 
                params["beta"].as<double>(),5, &Aregion);  
    for (int i=0; i!=20055000; i++){
        ssexy.MCstep();
        ssexy.Measure();
    }
    return 0;
} 
