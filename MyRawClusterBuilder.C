#include "MyRawClusterBuilder.h"
#include "RTHelper.h"
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
int MyRawClusterBuilder::Init/*Run*/(PHCompositeNode *topNode) {
    try {
        _CreateNodes(topNode);
    } catch (std::exception &e) {
        cout << PHWHERE << ": " << e.what() << endl;
        throw;
    }

    const string fileName = PATH + "rootFiles/rcb.root";
    _file = new TFile(fileName.c_str(),"RECREATE"); 
    /*_tree = new TTree("tree","Tree of RCB information.");
    _tree->Branch("energy",  &_f_energy);
    _tree->Branch("eta", &_f_eta);
    _tree->Branch("phi", &_f_phi);*/

    _tree = new TNtuple("tree", "cluster values", "energy:eta:phi");
    return Fun4AllReturnCodes::EVENT_OK;
}

/* ------------------------------------------------------------ *
 * MyRawClusterBuilder::process_event(...)                        *
 * ------------------------------------------------------------ */
int MyRawClusterBuilder::process_event(PHCompositeNode *topNode) {

    namespace IAlgorithm = IslandAlgorithm;
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
    std::list<RTHelper> seedTowers = IAlgorithm::GetSeedTowers(_towers, _towerGeom, _min_tower_e);
    cout << "seedTowers.size() = " << seedTowers.size() << endl;

    // Cluster the towers. 
    TowerMap clusteredTowers = IAlgorithm::GetClusteredTowers(seedTowers, _towers, _towerGeom);
    cout << "clusteredTowers.size() = " << clusteredTowers.size() << endl;
    
    // Fill _clusters (now empty) with the clusteredTowers and calculate their values.
    foreach (TowerPair& towerPair, clusteredTowers) {
        // Either get the cluster with this ID or create a new one.
        RawCluster *rawCluster = _clusters->getCluster(towerPair.first); 
        if (!rawCluster) _CreateNewCluster(&rawCluster);
        // Add the tower to this cluster.
        rawCluster->addTower(towerPair.second.getID(), towerPair.second.getEnergy());
        if (verbosity) _PrintCluster(towerPair);
    }

    cout << "_clusters->size() is " << _clusters->size() << endl;

    // Calculate/store energy, eta, phi of clusters given clusteredTowers information.
    _FillClustersEnergy(clusteredTowers);
    _FillClustersEta(clusteredTowers);
    _FillClustersPhi(clusteredTowers);
    for (unsigned i = 0; i < _clusters->size(); i++) {
        _AssignClusterValues(i);
        _tree->Fill(_energy[i], _eta[i], _phi[i]);
    }

    if (chkenergyconservation) _CheckEnergyConservation();
    return Fun4AllReturnCodes::EVENT_OK;

}

int MyRawClusterBuilder::End(PHCompositeNode *topNode) {
    cout << "WRITING TO FILE " << _file->GetName() << endl;
    _file->Write();
    _file->Close();
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
void MyRawClusterBuilder::_FillClustersEnergy(TowerMap clusteredTowers) {
    foreach (TowerPair& towerPair, clusteredTowers) {
        _energy[towerPair.first] += towerPair.second.getEnergy();
    }
}

// 2.
void MyRawClusterBuilder::_FillClustersEta(TowerMap clusteredTowers) {
    foreach (TowerPair& towerPair, clusteredTowers) {
        _eta[towerPair.first]   += towerPair.second.getEnergy() * towerPair.second.getEtaCenter();
    }
    for (unsigned i = 0; i < _clusters->size(); i++) {
        _eta[i] = (_energy[i] > 0) ? _eta[i] / _energy[i] : 0.0;
    }
}

// 3.
void MyRawClusterBuilder::_FillClustersPhi(TowerMap clusteredTowers) {
    // First, get all constitutent tower phi's as an energy-weighted sum.
    foreach (TowerPair& towerPair, clusteredTowers) {
        _phi[towerPair.first] += towerPair.second.getEnergy() * towerPair.second.getPhiCenter();
    }
    // Then divide by the total cluster energy.
    for (unsigned int i = 0; i < _clusters->size(); i++) {
        _phi[i] = (_energy[i] > 0) ? _phi[i] / _energy[i] : 0.0;
        if (_phi[i] > M_PI)  _phi[i] -= 2. * M_PI;
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

void MyRawClusterBuilder::_CreateNewCluster(RawCluster** rawCluster) {
    *rawCluster = new RawClusterv1();
    _clusters->AddCluster(*rawCluster);
    _energy.push_back(0.0);
    _eta.push_back(0.0);
    _phi.push_back(0.0);
}

// Serves to make ugly code in process_event less ugly.
int MyRawClusterBuilder::_NodeError(string nodeName, int retCode) {
    cout << PHWHERE << ": Could not find node " 
         << nodeName.data() << endl;
    return retCode;
}

void MyRawClusterBuilder::_PrintCluster(TowerPair towerPair) {
    cout << "MyRawClusterBuilder id: " << (towerPair.first) 
        << " Tower: " << " (iEta,iPhi) = (" << towerPair.second.getBinEta() 
        << "," << towerPair.second.getBinPhi() << ") " << " (eta,phi,e) = (" 
        << towerPair.second.getEtaCenter() << ","
        << towerPair.second.getPhiCenter() << ","
        //<< clusteredTower->get_energy() << ")"
        << endl;
}
