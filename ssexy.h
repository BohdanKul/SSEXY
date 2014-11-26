#ifndef SSEXY_H
#define SSEXY_H
#include "randombase.h"
#include "replica.h"
#include "lattice.h"
#include <vector>
#include <set>
#include <map>
#include <list>
//#include "timer.h"
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
        bool measTime;
        bool detVerbose;

        //Physical parameters
        int r; // order of Renyi entropy
        int Nx;
        int Ny;
        int N;
        int dim;
        float Beta;
        float T;

        //Measurements
        float SpinStiffness;
        float LRatio;
        double ALRatio;
        long  nAred;
        long  nAext;
        long  nAredRT;
        long  nAextRT;
        
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
        long         AnLoops;
        long         EnLoops;
        vector<long> shifts;
        vector<long> ns;
        vector<long> Tns;

        vector<long> NvisitedLegs;   //Number of visited legs per MC step
        vector<long> TNvisitedLegs;  //Accumulated number of visited legs per MC step
        vector<long> LINK; 
        vector<long> VTX; 
        vector<vector<long>> Partitions;
        map<long,vector<long>> DeterPaths;

        long LegSpin[6][4];
        long nTotal;
        long nMeas;
        long saveFreq;
        long nSaved; 
        bool Debug;
        bool SRTon;
        bool ILRTon;
        bool ALRTon;
        bool RandOffUpdate;

        //File management
        Communicator communicator;
        //Timer        timer;
        long SSEXID;

        long Bounce(long enLeg);
        long ContinueStraight(long enLeg);
        long SwitchReverse(long enLeg);
        long SwitchContinue(long enLeg);
        pair<long,long> SwitchLeg(long leg, long vtype, float prob);
        int BCnextSpin(int sindex, int& replica,bool connected);
        vector<long>* SwitchAregion();
        
        double ALRTrick();
        float ILRTrick();
        long  LoopPartition(vector<long>& BC);
        long  DeterministicOffDiagonalMove();
        long  RandomOffDiagonalUpdate();
        long  HardSwitchAregion();
        long  FlipLoop(vector<long>& BC, vector<long>& visited, long ispin, int ireplica){; 
        //long  GetConnectedSubraph(map<long,set<long>>& graph, set<long>& path,long cpos);

    public:
       SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, float _Beta, long seed, bool _measSS, bool _measTime, bool _detVerbose, int _Asize, string rfName, LATTICE* _Anor, LATTICE* _Ared, LATTICE* _Aext); 
       ~SSEXY();
       int   AdjustParameters();
       int   MCstep(); 
       int   Measure(); 
       int   SaveState();
       int   LoadState();
};
#endif
