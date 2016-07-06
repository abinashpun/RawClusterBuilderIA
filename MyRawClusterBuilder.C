#include "MyRawClusterBuilder.h"
#include "include/RTHelper.h"
#include "include/ClusterHelper.h"
#include "PHMakeGroups.h"

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
    return Fun4AllReturnCodes::EVENT_OK;
}

/* ------------------------------------------------------------ *
 * MyRawClusterBuilder::process_event(...)                        *
 * ------------------------------------------------------------ */
int MyRawClusterBuilder::process_event(PHCompositeNode *topNode) {
    // Clear any previously used helper objects. 
    ClusterHelper::NewEvent();
    string nodeName;

    // Grab the container of RawTowers (kinematic info).
    nodeName                    = "TOWER_CALIB_" + detector;
    RawTowerContainer* towers   = findNode::getClass<RawTowerContainer>(topNode, nodeName.c_str());
    if (!towers) return _NodeError(nodeName, Fun4AllReturnCodes::DISCARDEVENT);

    // Grab the container of RawTowerGeoms (tower geometry info).
    nodeName                    = "TOWERGEOM_" + detector;
    RTGeomContainer *towerGeom  = findNode::getClass<RTGeomContainer>(topNode,nodeName.c_str());
    if (!towerGeom) return _NodeError(nodeName, Fun4AllReturnCodes::ABORTEVENT);
    
    // Store the number of bins in phi as a static value, as it should be. 
    RTHelper::set_maxphibin(towerGeom->get_phibins());

    // make the list of towers above threshold
    vector<RTHelper> seedTowers;
    RTCItr itr;
    RTCRange itrPair = towers->getTowers();
    for (itr = itrPair.first; itr != itrPair.second; itr++) {
        // Note: itr has the form of pair(int id, RawTower*).
        if (itr->second->get_energy() > _min_tower_e) {
            _InsertSeed(seedTowers, itr, towerGeom);
        }
    }
    cout << "seedTowers.size() = " << seedTowers.size() << endl;

    // Cluster the towers. 
    std::multimap<int, RTHelper> clusteredTowers;
    PHMakeGroups(seedTowers, clusteredTowers);

    ClusterItr ctitr  = clusteredTowers.begin();
    ClusterItr lastct = clusteredTowers.end();
    // Fill _clusters (now empty) with the clusteredTowers and calculate their values.
    for (; ctitr != lastct; ++ctitr) {
        // Store this cluster's id and the associated RawTower.
        int clusterID            = ctitr->first;
        RawTower *clusteredTower = RTHelper::GetRawTower(ctitr->second, towers);

        // If this tower belongs to a cluster we haven't seen yet, then
        // add the new cluster to _clusters and push_back 0.0 for eta, phi, and energy.
        RawCluster *rawCluster = _clusters->getCluster(clusterID); 
        if (!rawCluster) ClusterHelper::NewCluster(rawCluster, _clusters);

        // Finally, add the tower to this cluster.
        rawCluster->addTower(clusteredTower->get_id(), clusteredTower->get_energy());
        if (verbosity) _PrintCluster(ctitr);
    }

    // Get cluster id-indexed properties.
    ClusterHelper::SetNClusters(_clusters->size());
    vector<float> clustersEnergy = ClusterHelper::GetClustersEnergy(clusteredTowers, towers);
    vector<float> clustersEta    = ClusterHelper::GetClustersEta(clusteredTowers, towers);
    vector<float> clustersPhi    = ClusterHelper::GetClustersPhi(clusteredTowers, towers, _clusters, towerGeom);

    unsigned nClusters = _clusters->size();
    for (unsigned iCluster = 0; iCluster < nClusters; iCluster++) {
        // Set energy, eta, phi of the cluster.
        float e   = clustersEnergy[iCluster];
        float eta = clustersEta[iCluster];
        float phi = clustersPhi[iCluster];
        _AssignClusterValues(iCluster, e, eta, phi);
    }

    if (chkenergyconservation) _CheckEnergyConservation(towers);     
    
    return Fun4AllReturnCodes::EVENT_OK;
}

int MyRawClusterBuilder::End(PHCompositeNode *topNode) {
    return Fun4AllReturnCodes::EVENT_OK;
}

// ----------------------------------------------------------------------------
// Private helper methods.
// ----------------------------------------------------------------------------

void MyRawClusterBuilder::_CheckEnergyConservation(RTContainer* towers) {
    double ecluster = _clusters->getTotalEdep();
    double etower   = towers->getTotalEdep();
    if (ecluster > 0 && (fabs(etower - ecluster) / ecluster) > 1e-9) {
        cout << "energy conservation violation: ETower: " << etower
             << " ECluster: " << ecluster 
             << " diff: " << etower - ecluster << endl;
    } else if (etower != 0) {
        cout << "energy conservation violation: ETower: " << etower
             << " ECluster: " << ecluster << endl;
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
    PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(
            iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode) {
        std::cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
        throw std::runtime_error(
                "Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
    }
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(
            dstiter.findFirst("PHCompositeNode",detector ));
    if(!DetNode){
        DetNode = new PHCompositeNode(detector);
        dstNode->addNode(DetNode);
    }
    _clusters = new RawClusterContainer();
    string ClusterNodeName = "CLUSTER_" + detector;
    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(
            _clusters, ClusterNodeName.c_str(), "PHObject");
    DetNode->addNode(clusterNode);
}

// Serves to make ugly code in process_event less ugly.
int MyRawClusterBuilder::_NodeError(string nodeName, int retCode) {
    cout << PHWHERE << ": Could not find node " 
         << nodeName.data() << endl;
    return retCode;
}

void MyRawClusterBuilder::_AssignClusterValues(int iCluster, float e, float eta, float phi) {
    RawCluster* cluster = _clusters->getCluster(iCluster);
    cluster->set_energy(e);
    cluster->set_eta(eta);
    cluster->set_phi(phi);
    if (verbosity) {
        cout << " (eta,phi,e) = (" << cluster->get_eta() << ", "
             << cluster->get_phi() << ","
             << cluster->get_energy() << ")"
             << endl;
    }
}

// Given iterator to a seed tower, place relevant info into vector of seed towers.
void MyRawClusterBuilder::_InsertSeed(vector<RTHelper>&  vec, 
                                     RTCItr             seedItr, 
                                     RTGeomContainer*    towerGeom) {
        // Store the RTC pair elements separately.
        KeyType seedID      = seedItr->first;
        RawTower* seedTower = seedItr->second;
        // Wrap RawTower info inside an RTHelper. 
        RTHelper rtHelper(seedTower);
        rtHelper.set_id(seedID);
        rtHelper.setCenter(towerGeom);
        vec.push_back(rtHelper);
}

void MyRawClusterBuilder::_PrintCluster(ClusterItr ctitr) {
    cout << "MyRawClusterBuilder id: " << (ctitr->first) 
        << " Tower: " << " (iEta,iPhi) = (" << ctitr->second.get_bineta() 
        << "," << ctitr->second.get_binphi() << ") " << " (eta,phi,e) = (" 
        << ctitr->second.getEtaCenter() << ","
        << ctitr->second.getPhiCenter() << ","
        //<< clusteredTower->get_energy() << ")"
        << endl;
}
