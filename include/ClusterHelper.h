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
typedef std::multimap<int, RTHelper> TowerMap;
typedef TowerMap::iterator ClusterItr;
typedef RawTowerContainer RTContainer;
typedef RawTowerGeomContainer RTGeomContainer;

class ClusterHelper {
    public:
        static void NewEvent();
        static void NewCluster(RawCluster*, RawClusterContainer*);
        static vector<float> GetClustersEnergy(TowerMap, RTContainer*);
        static vector<float> GetClustersEta(TowerMap, RTContainer*);
        static vector<float> GetClustersPhi(TowerMap, RTContainer*, RawClusterContainer*, RTGeomContainer*);
        static void SetNClusters(int n);
        static bool _CorrectPhi(RawCluster*, RTContainer*, RTGeomContainer*);
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

vector<float> ClusterHelper::GetClustersPhi(TowerMap clusteredTowers, RTContainer* towers, RawClusterContainer* clusters, RTGeomContainer* towerGeom) {
    ClusterItr ctitr = clusteredTowers.begin();
    ClusterItr lastct = clusteredTowers.end();
    // First, get all constitutent tower phi's as an energy-weighted sum.
    for (; ctitr != lastct; ctitr++) {
        RawTower *rawTower = RTHelper::GetRawTower(ctitr->second, towers);
        phi[ctitr->first] += rawTower->get_energy() * ctitr->second.getPhiCenter();
    }
    // Then divide by the total cluster energy.
    for (int i = 0; i < nClusters; i++) {
        if (energy[i] > 0)   phi[i] /= energy[i];
        else                        phi[i] = 0.0;
        if (phi[i] > M_PI) phi[i] -= 2. * M_PI;
    }
    // Finally, correct the mean Phi calculation for clusters at Phi discontinuity
    // Assumes that Phi goes from -pi to +pi
    for (int iCluster = 0; iCluster < nClusters; iCluster++) {
        RawCluster *cluster = clusters->getCluster(iCluster);
        float oldphi = cluster->get_phi();
        bool corr = _CorrectPhi(cluster, towers, towerGeom);
        if (corr) {
            cout << PHWHERE << " Cluster Phi corrected: " 
                      << oldphi << " " << cluster->get_phi() << endl;
        }
    }
    return phi;
}

bool ClusterHelper::_CorrectPhi(RawCluster* cluster, RTContainer* towers, RTGeomContainer *towerGeom) {
    double sum = cluster->get_energy();
    double phimin = 999.;
    double phimax = -999.;
    RawCluster::TowerConstRange begin_end = cluster->get_towers();
    RawCluster::TowerConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter) { 
        RawTower* tmpt = towers->getTower(iter->first);
        double phi = towerGeom->get_phicenter(tmpt->get_binphi());
        if(phi > M_PI) phi = phi - 2.*M_PI; 
        if (phi < phimin) {
            phimin = phi;
        }
        if (phi > phimax) {
            phimax = phi;
        }
    }

    // cluster is not at phi discontinuity
    if ((phimax - phimin) < 3.) return false; 

    float mean = 0.;
    for (iter =begin_end.first; iter != begin_end.second; ++iter) { 
        RawTower* tmpt = towers->getTower(iter->first);
        double e = tmpt->get_energy();
        double phi = towerGeom->get_phicenter(tmpt->get_binphi());
        if(phi > M_PI) phi = phi - 2.*M_PI; 
        if (phi < 0.) {
            phi = phi + 2.*M_PI;  // shift phi range for correct mean calculation
        }
        mean += e * phi;
    }
    mean = mean / sum;
    if (mean > M_PI) {
        mean = mean - 2.*M_PI;  // shift back
    }
    cluster->set_phi(mean);
    return true; // mean phi was corrected
}

#endif
