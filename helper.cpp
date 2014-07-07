#ifndef HELPER_CPP
#define HELPER_CPP
 
#include "stdlib.h"
#include <vector>
#include <set>
#include <map>
using namespace std;

//**************************************************************************
// Merge vectors (stored by reference in vectors)
vector<long> MergeVectors(vector<vector<long>*> vectors){
    vector<long> merged;
    long totalSize = 0;
    //Calculate the total size of resulting merged vector
    for (long i=0; i!=vectors.size(); i++){
        totalSize += vectors[i]->size();
    }

    //Reserve enough memory for it
    merged.reserve(totalSize);

    //Merge in sequence
    for (long i=0; i!=vectors.size(); i++){
        merged.insert(merged.end(),vectors[i]->begin(),vectors[i]->end());
    }
    return merged;
}

//**************************************************************************

long GetConnectedSubraph(map<long,set<long>>& graph, set<long>& path,long cpos){

    path.insert(cpos);
    for (auto npos=graph[cpos].begin(); npos!=graph[cpos].end(); npos++)
        if  (not path.count(*npos))
            GetConnectedSubraph(graph,path,*npos);
}

long sgn(long val) {
    return ((0<=val) - (val<0));
}
#endif
