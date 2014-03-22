#ifndef SSEXY_H
#define SSEXY_H
#include "randombase.h"
#include "replica.h"
#include <vector>
using namespace std;


//Workhorse
class SSEXY: public RandomBase{


    private:
        //Algorithmic parameters
        int  Nloops;
        int  binSize;
        long maxLoopSize;
        bool measSS;
        bool measRatio;

        //Physical parameters
        int r; // order of Renyi entropy
        int Nx;
        int Ny;
        int N;
        float Beta;
        float T;

        //Measurements
        float SpinStiffness;
        float ZRatio;
        long  nAred;
        long  nAext;
        
        vector<Replica*> Replicas; 
        vector<vector<long>> sites;
        vector<vector<long>*> firsts ;
        vector<vector<long>*> lasts  ;
        vector<vector<long>*> links  ;
        vector<vector<long>*> spins  ;
        vector<vector<long>*> vtxs   ;
        vector<vector<long>*> sms    ;
        
        vector<long> ap;
        vector<long>* Aregion;
        vector<long> Aext;
        vector<long> Ared;
        vector<long> Adif;
        vector<long> shifts;
        vector<long> ns;
        vector<long> Tns;

        vector<long> NvisitedLegs;   //Number of visited legs per MC step
        vector<long> TNvisitedLegs;  //Accumulated number of visited legs per MC step
        vector<long> LINK; 
        vector<long> VTX; 

        long LegSpin[6][4];
        long nTotal;
        long nMeas;
        long saveFreq;
        long nSaved; 
        bool Debug;
        bool DebugSRT;

        //File management
        Communicator communicator;
        long SSEXID;

        long Bounce(long enLeg);
        long ContinueStraight(long enLeg);
        long SwitchReverse(long enLeg);
        long SwitchContinue(long enLeg);
        pair<long,long> SwitchLeg(long leg, long vtype);
        int BCnextSpin(int sindex, int& replica,bool connected);
        vector<long>* SwitchAregion();
        
        float MeasureZRatio();
        long  MeasureNLoop(vector<long>& BC);

    public:
       SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, float _Beta, long seed, bool _measSS, int _Asize, string rfName, vector<long>* _Anor, vector<long>* _Ared, vector<long>* _Aext); 
       int   AdjustParameters();
       int   MCstep(); 
       int   Measure(); 
       int   SaveState();
       int   LoadState();
};
#endif
