#ifndef CLUSTERHELPER_H_
#define CLUSTERHELPER_H_
// C++ includes.
#include <iostream>
#include <vector>
// Local includes.
#include "RTHelper.h"
#include "g4cemc/RawCluster.h"
#include "g4cemc/RawClusterContainer.h"
#include "g4cemc/RawTowerContainer.h"
#include "g4cemc/RawClusterContainer.h"

using std::vector;
typedef RawTowerContainer RTContainer;
typedef RawTowerGeomContainer RTGeomContainer;

class ClusterHelper {
    public:
        static void NewEvent();
        static void NewCluster(RawCluster*, RawClusterContainer*);
        static void SetNClusters(int n);
        static vector<float> energy;
        static vector<float> eta;
        static vector<float> phi;
        static int nClusters;
};
int ClusterHelper::nClusters = 0;
vector<float> ClusterHelper::energy;
vector<float> ClusterHelper::eta;
vector<float> ClusterHelper::phi;

void ClusterHelper::SetNClusters(int n) {
    nClusters = n;
}

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
