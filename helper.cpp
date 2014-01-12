#ifndef HELPER_CPP
#define HELPER_CPP
 
#include "stdlib.h"
#include <vector>
using namespace std;

//**************************************************************************
// Merge vectors (stored by reference in vectors)
vector<long> MergeVectors(vector<vector<long>*> vectors){
    vector<long> merged;
    int totalSize = 0;
    //Calculate the total size of resulting merged vector
    for (int i=0; i!=vectors.size(); i++){
        totalSize += vectors[i]->size();
    }

    //Reserve enough memory for it
    merged.reserve(totalSize);

    //Merge in sequence
    for (int i=0; i!=vectors.size(); i++){
        merged.insert(merged.end(),vectors[i]->begin(),vectors[i]->end());
    }
    return merged;
}
//**************************************************************************

long sgn(long val) {
    return ((0<=val) - (val<0));
}
#endif
