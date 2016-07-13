#ifndef __MYRAWCLUSTERBUILDER_H__
#define __MYRAWCLUSTERBUILDER_H__

// C++ includes.
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <list>
// Fun4All/PHENIX includes.
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
// Local includes.
#include "g4cemc/RawTower.h"
#include "g4cemc/RawTowerGeomContainer.h"
#include "g4cemc/RawTowerContainer.h"
#include "include/RawClusterv1.h"
#include "g4cemc/RawClusterContainer.h"
// ROOT Includes.
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

using std::cout;
using std::endl;
using std::string;

const float GEN_ENERGY = 20.0; // TODO: change when know how to get gen info. 
const string PATH = "/phenix/u/bmckinz/sphenix_src/";

// Forward class declarations.
class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class RTHelper;
class ClusterHelper;

typedef RawTowerContainer                        RTContainer;
typedef RawTowerGeomContainer                    RTGeomContainer;
typedef std::multimap<int, RTHelper>             TowerMap;
typedef std::pair<const int, RTHelper>           TowerPair;
typedef std::pair<const unsigned int, RawTower*> RawTowerPair;


class MyRawClusterBuilder : public SubsysReco {
    public:
        MyRawClusterBuilder(const std::string& name = "MyRawClusterBuilder"); 
        virtual ~MyRawClusterBuilder() {}
        //int InitRun(PHCompositeNode *topNode);
        int Init(PHCompositeNode *topNode);
        int process_event(PHCompositeNode *topNode);
        int End(PHCompositeNode *topNode);
        void Detector(const string &d)              { detector = d; }
        void set_threshold_energy(const float e)    { _min_tower_e = e; }
        void checkenergy(const int i = 1)           { chkenergyconservation = i; }
    private:
        RawClusterContainer*_clusters;
        RTContainer*        _towers;
        RTGeomContainer*    _towerGeom;
        float               _min_tower_e;
        int                 chkenergyconservation;
        string              detector;
        std::vector<float>  _energy; 
        std::vector<float>  _eta; 
        std::vector<float>  _phi;

        // ROOT I/O Objects. 
        TFile* _file;
        TNtuple* ntp_cluster;
        TNtuple* ntp_tower;

        int  _NodeError(string nodeName, int retCode);
        void _AssignClusterValues(int iCluster);
        void _CreateNodes(PHCompositeNode *topNode);
        void _PrintCluster(TowerPair);
        void _CheckEnergyConservation();
        void _FillClustersEnergy(TowerMap);
        void _FillClustersEta(TowerMap);
        void _FillClustersPhi(TowerMap);
        std::list<RTHelper> _GetAllTowers();
        void _InsertTower(std::list<RTHelper>&, RawTowerPair);
        bool _CorrectPhi(RawCluster*);
        void _CreateNewCluster(RawCluster**);
        void _ShowTreeEntries();
        void _FillTowerTree(std::list<RTHelper>);
        void _FillClusterTree();
};

#endif
