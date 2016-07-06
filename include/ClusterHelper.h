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

using std::vector;
typedef std::multimap<int, RTHelper> TowerMap;
typedef TowerMap::iterator ClusterItr;
typedef RawTowerContainer RTContainer;

class ClusterHelper {
    public:
        static void NewEvent();
        static void NewCluster(RawCluster*, RawClusterContainer*);
        static vector<float> GetClustersEnergy(TowerMap, RTContainer*);
        static vector<float> GetClustersEta(TowerMap, RTContainer*);
        static vector<float> GetClustersPhi(TowerMap, RTContainer*);
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

vector<float> ClusterHelper::GetClustersEnergy(TowerMap clusteredTowers, RTContainer* towers) {
    ClusterItr ctitr = clusteredTowers.begin();
    ClusterItr lastct = clusteredTowers.end();
    for (; ctitr != lastct; ctitr++) {
        energy[ctitr->first] += RTHelper::GetRawTower(ctitr->second, towers)->get_energy();
    }
    return energy;
}

vector<float> ClusterHelper::GetClustersEta(TowerMap clusteredTowers, RTContainer* towers) {
    ClusterItr ctitr = clusteredTowers.begin();
    ClusterItr lastct = clusteredTowers.end();
    for (; ctitr != lastct; ctitr++) {
        RawTower *rawTower   = RTHelper::GetRawTower(ctitr->second, towers);
        eta[ctitr->first]   += rawTower->get_energy() * ctitr->second.getEtaCenter();
    }
    for (int i = 0; i < nClusters; i++) {
        if (energy[i] > 0)   eta[i] /= energy[i];
        else                        eta[i] = 0.0;
    }
    return eta;
}

vector<float> ClusterHelper::GetClustersPhi(TowerMap clusteredTowers, RTContainer* towers) {
    ClusterItr ctitr = clusteredTowers.begin();
    ClusterItr lastct = clusteredTowers.end();
    for (; ctitr != lastct; ctitr++) {
        RawTower *rawTower = RTHelper::GetRawTower(ctitr->second, towers);
        phi[ctitr->first] += rawTower->get_energy() * ctitr->second.getPhiCenter();
    }
    for (int i = 0; i < nClusters; i++) {
        if (energy[i] > 0)   phi[i] /= energy[i];
        else                        phi[i] = 0.0;
        if (phi[i] > M_PI) phi[i] -= 2. * M_PI;
    }
    return phi;
}

#endif
