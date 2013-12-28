#ifndef SSEXY_H
#define SSEXY_H
#include "randombase.h"
#include "replica.h"
#include <vector>
using namespace std;


//Workhorse
class SSEXY: public RandomBase{


    private:
        int r; // order of Renyi entropy
        int Nloops;
        int binSize;
        int Nx;
        int Ny;
        int N;
        float Beta;
        float T;
        vector<vector<long>> sites;
        vector<long> ap;
        vector<Replica*> Replicas; 
        vector<long>  Aregion;
        vector<vector<long>*> firsts ;
        vector<vector<long>*> lasts  ;
        vector<vector<long>*> links  ;
        vector<vector<long>*> spins  ;
        vector<vector<long>*> vtxs   ;
        vector<vector<long>*> sms    ;
        vector<long> shifts;
        vector<long> ns;
        vector<long> Tns;
        
        vector<long> NvisitedLegs;  //Number of visited legs per MC step
        vector<long> LINK; 
        vector<long> VTX; 

        long LegSpin[6][4];
        long nTotal;
        long nMeas;
        bool Debug;
        //File management
        Communicator communicator;

        long Bounce(long enLeg);
        long ContinueStraight(long enLeg);
        long SwitchReverse(long enLeg);
        long SwitchContinue(long enLeg);
        pair<long,long> SwitchLeg(long leg, long vtype);

    public:
       SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, float _Beta, long seed, vector<long>* _Aregion); 
       int  MCstep(); 
       int Measure(); 
};
#endif
