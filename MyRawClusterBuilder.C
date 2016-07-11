#include "MyRawClusterBuilder.h"
#include "RTHelper.h"
#include "ClusterHelper.h"
//#include "PHMakeGroups.h"
#include "IslandAlgorithm.h"
#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

/* ------------------------------------------------------ *
 * MyRawClusterBuilder::MyRawClusterBuilder()                   *
 * Initialize values of first four attributes and           *
 * call parent constructor (SubsysReco).                    *
 * -------------------------------------------------------- */
MyRawClusterBuilder::MyRawClusterBuilder(const string& name)
    : SubsysReco(name),
      _clusters(NULL),
      _min_tower_e(0.0),
      chkenergyconservation(0),
      detector("NONE") {}

/* ----------------------------------------------------- *
 * MyRawClusterBuilder::InitRun()                          *
 * ----------------------------------------------------- */
int MyRawClusterBuilder::InitRun(PHCompositeNode *topNode) {
    try {
        _CreateNodes(topNode);
    } catch (std::exception &e) {
        cout << PHWHERE << ": " << e.what() << endl;
        throw;
    }

    _file = new TFile("rootFiles/rcb.root","RECREATE"); 
    _tree = new TTree("tree","Tree of RCB information.");

    _tree->Branch("energy",  &_f_energy);
    _tree->Branch("eta", &_f_eta);
    _tree->Branch("phi", &_f_phi);
    return Fun4AllReturnCodes::EVENT_OK;
}

/* ------------------------------------------------------------ *
 * MyRawClusterBuilder::process_event(...)                        *
 * ------------------------------------------------------------ */
int MyRawClusterBuilder::process_event(PHCompositeNode *topNode) {
    namespace IAlgorithm = IslandAlgorithm;

    // Clear any previously used helper objects. 
    ClusterHelper::NewEvent();
    string nodeName;

    // Grab the container of RawTowers (kinematic info).
    nodeName = "TOWER_CALIB_" + detector;
    _towers  = findNode::getClass<RTContainer>(topNode, nodeName.c_str());
    if (!_towers) return _NodeError(nodeName, Fun4AllReturnCodes::DISCARDEVENT);

    // Grab the container of RawTowerGeoms (tower geometry info).
    nodeName    = "TOWERGEOM_" + detector;
    _towerGeom  = findNode::getClass<RTGeomContainer>(topNode, nodeName.c_str());
    if (!_towerGeom) return _NodeError(nodeName, Fun4AllReturnCodes::ABORTEVENT);
    
    // Store the number of bins in phi as a static value, as it should be. 
    RTHelper::setMaxPhiBin(_towerGeom->get_phibins());
    RTHelper::setMaxEtaBin(_towerGeom->get_etabins());

    // Make the list of _towers above minimum energy threshold.
    set_threshold_energy(0.1);
    std::list<RTHelper> allTowers = _GetAllTowers();
    std::list<RTHelper> seedTowers = IAlgorithm::GetSeedTowers(_towers, _towerGeom, _min_tower_e);
    cout << "seedTowers.size() = " << seedTowers.size() << endl;

    // Cluster the towers. 
    TowerMap clusteredTowers = IAlgorithm::GetClusteredTowers(seedTowers, _towers, _towerGeom);
    cout << "clusteredTowers.size() = " << clusteredTowers.size() << endl;
    
    // Fill _clusters (now empty) with the clusteredTowers and calculate their values.
    foreach (TowerPair& ctitr, clusteredTowers) {
        // Store this cluster's id and the associated RawTower.
        int clusterID            = ctitr.first;
        RawTower *clusteredTower = RTHelper::GetRawTower(ctitr.second, _towers);
        // If this tower belongs to a cluster we haven't seen yet, then
        // add the new cluster to _clusters and push_back 0.0 for eta, phi, and energy.
        RawCluster *rawCluster = _clusters->getCluster(clusterID); 
        if (!rawCluster) ClusterHelper::NewCluster(rawCluster, _clusters);
        // Finally, add the tower to this cluster.
        rawCluster->addTower(clusteredTower->get_id(), clusteredTower->get_energy());
        if (verbosity) _PrintCluster(ctitr);
    }

    // Calculate/store energy, eta, phi of clusters given clusteredTowers information.
    _energy = _GetClustersEnergy(clusteredTowers);
    _eta    = _GetClustersEta(clusteredTowers);
    _phi    = _GetClustersPhi(clusteredTowers);
    for (unsigned i = 0; i < _clusters->size(); i++) {
        _AssignClusterValues(i);
    }

    if (chkenergyconservation) _CheckEnergyConservation();
    return Fun4AllReturnCodes::EVENT_OK;

}

int MyRawClusterBuilder::End(PHCompositeNode *topNode) {
    return Fun4AllReturnCodes::EVENT_OK;
}

// -----------------------------------------------------------------------------
// -------------------------- PRIVATE HELPER METHODS. --------------------------
// -----------------------------------------------------------------------------

void MyRawClusterBuilder::_AssignClusterValues(int iCluster) {
    RawCluster* cluster = _clusters->getCluster(iCluster);
    cluster->set_energy(_energy[iCluster]);
    cluster->set_eta(_eta[iCluster]);
    cluster->set_phi(_phi[iCluster]);
    if (verbosity) {
        cout << " (eta,phi,e) = (" << cluster->get_eta() << ", "
             << cluster->get_phi() << ","
             << cluster->get_energy() << ")"
             << endl;
    }
}

// Serves to make ugly code in process_event less ugly.
int MyRawClusterBuilder::_NodeError(string nodeName, int retCode) {
    cout << PHWHERE << ": Could not find node " 
         << nodeName.data() << endl;
    return retCode;
}

void MyRawClusterBuilder::_PrintCluster(TowerPair ctitr) {
    cout << "MyRawClusterBuilder id: " << (ctitr.first) 
        << " Tower: " << " (iEta,iPhi) = (" << ctitr.second.getBinEta() 
        << "," << ctitr.second.getBinPhi() << ") " << " (eta,phi,e) = (" 
        << ctitr.second.getEtaCenter() << ","
        << ctitr.second.getPhiCenter() << ","
        //<< clusteredTower->get_energy() << ")"
        << endl;
}

// Return list of all RawTower pairs in _towers->getTowers() converted to RTHelpers.
std::list<RTHelper> MyRawClusterBuilder::_GetAllTowers() {
    std::list<RTHelper> allTowers;
    foreach (RawTowerPair& towerPair, _towers->getTowers()) {
        // TODO : change order of arguments. 
        _InsertTower(allTowers, towerPair);
    }
    return allTowers;
}

// Given iterator to a seed tower, place relevant info into std::vector of seed towers.
void MyRawClusterBuilder::_InsertTower(std::list<RTHelper>&  towerList, RawTowerPair towerPair)  {
    RTHelper rtHelper(towerPair.second);
    rtHelper.setCenter(_towerGeom);
    towerList.push_back(rtHelper);
}


// 1.
std::vector<float> MyRawClusterBuilder::_GetClustersEnergy(TowerMap clusteredTowers) {
    std::vector<float> energy;
    foreach (TowerPair& ctitr, clusteredTowers) {
        energy[ctitr.first] += RTHelper::GetRawTower(ctitr.second, _towers)->get_energy();
    }
    return energy;
}

// 2.
std::vector<float> MyRawClusterBuilder::_GetClustersEta(TowerMap clusteredTowers) {
    std::vector<float> eta;
    foreach (TowerPair& ctitr, clusteredTowers) {
        RawTower *rawTower  = RTHelper::GetRawTower(ctitr.second, _towers);
        eta[ctitr.first]   += rawTower->get_energy() * ctitr.second.getEtaCenter();
    }

    for (unsigned int i = 0; i < _clusters->size(); i++) {
        if (_energy[i] > 0) eta[i] /= _energy[i];
        else                eta[i] = 0.0;
    }
    return eta;
}

// 3.
std::vector<float> MyRawClusterBuilder::_GetClustersPhi(TowerMap clusteredTowers) {
    std::vector<float> phi;
    // First, get all constitutent tower phi's as an energy-weighted sum.
    foreach (TowerPair& ctitr, clusteredTowers) {
        RawTower *rawTower = RTHelper::GetRawTower(ctitr.second, _towers);
        phi[ctitr.first] += rawTower->get_energy() * ctitr.second.getPhiCenter();
    }

    // Then divide by the total cluster energy.
    for (unsigned int i = 0; i < _clusters->size(); i++) {
        if (_energy[i] > 0)  phi[i] /= _energy[i];
        else                phi[i] = 0.0;
        if (phi[i] > M_PI)  phi[i] -= 2. * M_PI;
    }

    // Finally, correct the mean Phi calculation for clusters at Phi discontinuity.
    for (unsigned int iCluster = 0; iCluster < _clusters->size(); iCluster++) {
        RawCluster *cluster = _clusters->getCluster(iCluster);
        float oldPhi        = cluster->get_phi();
        bool corr           = _CorrectPhi(cluster);
        if (corr) {
            cout << PHWHERE << " Cluster Phi corrected: " 
                 << oldPhi  << " " << cluster->get_phi() << endl;
        }
    }
    return phi;
}

/* ----------------------------------------------------- *
 * MyRawClusterBuilder::CreateNodes()                          *
 * Initialize values of first four attributes and        *
 * call parent constructor (SubsysReco).                 *
 * ----------------------------------------------------- */
void MyRawClusterBuilder::_CreateNodes(PHCompositeNode *topNode) {
    PHNodeIterator iter(topNode);
    // Grab the cEMC node
    PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode) {
        std::cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
        throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
    }
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode",detector ));
    if(!DetNode){
        DetNode = new PHCompositeNode(detector);
        dstNode->addNode(DetNode);
    }

    _clusters       = new RawClusterContainer();
    string nodeName = "CLUSTER_" + detector;
    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, nodeName.c_str(), "PHObject");
    DetNode->addNode(clusterNode);
}


bool MyRawClusterBuilder::_CorrectPhi(RawCluster* cluster) {
    double sum      = cluster->get_energy();
    double phimin   = 999.;
    double phimax   = -999.;
    RawCluster::TowerConstRange begin_end = cluster->get_towers();
    RawCluster::TowerConstIterator iter;
    for (iter = begin_end.first; iter != begin_end.second; ++iter) { 
        RawTower* tmpt = _towers->getTower(iter->first);
        double phi = _towerGeom->get_phicenter(tmpt->get_binphi());
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
        RawTower* tmpt = _towers->getTower(iter->first);
        double e = tmpt->get_energy();
        double phi = _towerGeom->get_phicenter(tmpt->get_binphi());
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

void MyRawClusterBuilder::_CheckEnergyConservation() {
    double ecluster = _clusters->getTotalEdep();
    double etower   = _towers->getTotalEdep();
    if (ecluster > 0 && (fabs(etower - ecluster) / ecluster) > 1e-9) {
        cout << "energy conservation violation: ETower: " << etower
             << " ECluster: " << ecluster 
             << " diff: " << etower - ecluster << endl;
    } else if (etower != 0) {
        cout << "energy conservation violation: ETower: " << etower
             << " ECluster: " << ecluster << endl;
    }
}

