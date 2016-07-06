#ifndef CLUSTERHELPER_H_
#define CLUSTERHELPER_H_
// C++ includes.
#include <iostream>
#include <vector>
// Local includes.
#include "g4cemc/RawCluster.h"
#include "g4cemc/RawClusterContainer.h"

using std::vector;

class ClusterHelper {
    public:
        static void NewEvent();
        static void NewCluster(RawCluster*, RawClusterContainer*);
        static vector<float> energy;
        static vector<float> eta;
        static vector<float> phi;
};
vector<float> ClusterHelper::energy;
vector<float> ClusterHelper::eta;
vector<float> ClusterHelper::phi;


void ClusterHelper::NewEvent() {
    energy.clear();
    eta.clear();
    phi.clear();
}

void ClusterHelper::NewCluster(RawCluster* cluster, RawClusterContainer* _clusters) {
    cluster = new RawClusterv1();
    _clusters->AddCluster(cluster);
    energy.push_back(0.0);
    eta.push_back(0.0);
    phi.push_back(0.0);
}

#endif
