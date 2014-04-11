#ifndef SSEXY_H
#define SSEXY_H
#include "randombase.h"
#include "replica.h"
#include "lattice.h"
#include <vector>
#include <set>
#include <map>
#include <list>
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
        float LRatio;
        float ALRatio;
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
        long SSEXID;

        long Bounce(long enLeg);
        long ContinueStraight(long enLeg);
        long SwitchReverse(long enLeg);
        long SwitchContinue(long enLeg);
        pair<long,long> SwitchLeg(long leg, long vtype, float prob);
        int BCnextSpin(int sindex, int& replica,bool connected);
        vector<long>* SwitchAregion();
        
        float ALRTrick();
        float ILRTrick();
        long  LoopPartition(vector<long>& BC);
        long  DeterministicOffDiagonalMove();
        long  RandomOffDiagonalUpdate();
        //long  GetConnectedSubraph(map<long,set<long>>& graph, set<long>& path,long cpos);

    public:
       SSEXY(int _r, unsigned short _Nx, unsigned short _Ny, float _T, float _Beta, long seed, bool _measSS, int _Asize, string rfName, LATTICE* _Anor, LATTICE* _Ared, LATTICE* _Aext); 
       int   AdjustParameters();
       int   MCstep(); 
       int   Measure(); 
       int   SaveState();
       int   LoadState();
};
#endif
