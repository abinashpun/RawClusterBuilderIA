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
// ROOT Includes.
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TROOT.h"

// Forward class declarations.
class RawTower;
class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class TowerHelper;

typedef RawTowerContainer       RTContainer;
typedef RawTowerGeomContainer   RTGeomContainer;

typedef std::multimap<int, TowerHelper>             TowerMap;
typedef std::pair<const int, TowerHelper>           TowerPair;
typedef std::pair<const unsigned int, RawTower*> RawTowerPair;


class MyRawClusterBuilder : public SubsysReco {
    public:
        MyRawClusterBuilder(const std::string& name = "MyRawClusterBuilder"); 
        virtual ~MyRawClusterBuilder() {}
        int Init(PHCompositeNode *topNode);
        int process_event(PHCompositeNode *topNode);
        int End(PHCompositeNode *topNode);
        void Detector(const std::string &d)              { _detector = d; }
        void set_threshold_energy(const float e)    { _min_tower_e = e; }
        void checkenergy(const int i = 1)           { _checkEnergyConserv = i; }
        void SetGenPT(float pt)                     { _genPT = pt; }
        void SetParticleType(std::string s)              { _particleType = s; }
        void SetEvent(int i)                        { _iEvent = i; }
        void ClusterSimple(bool b)                  { _clusterSimple = b; }
    private:
        // Variables initialized in constructor list.
        RawClusterContainer*_clusters;
        float       _min_tower_e;
        int         _checkEnergyConserv;
        bool        _clusterSimple;
        std::string      _fileName;
        std::string      _detector;

        int         _iEvent;
        // Other sPHENIX private variables.
        RTContainer*        _towers;
        RTGeomContainer*    _towerGeom;

        std::string  _particleType;

        // ROOT I/O Objects. 
        TFile*   _file;
        TNtuple* ntp_tower;
        TTree* _tCluster;
        std::vector<int> towerIDs;
        int     _clusterID;
        float   _genPT;
        float   _f_energy;
        float   _f_ET;
        float   _f_eta;
        float   _f_phi;
        int     _nClusters;
        int     _nTowers;

        std::vector<float>  _energyVec; 
        std::vector<float>  _ETVec; 
        std::vector<float>  _etaVec; 
        std::vector<float>  _phiVec;

        // Private helper methods. 
        void _AssignClusterValues(int iCluster);
        void _PrintCluster(TowerPair);
        void _CheckEnergyConservation();
        void _FillClustersEnergy(TowerMap);
        void _FillClustersEta(TowerMap);
        void _FillClustersPhi(TowerMap);
        std::list<TowerHelper> _GetAllTowers();
        void _InsertTower(std::list<TowerHelper>&, RawTowerPair);
        bool _CorrectPhi(RawCluster*);
        void _CreateNewCluster(RawCluster**);
        void _ShowTreeEntries();
        void _FillTowerTree(std::list<TowerHelper>);
        void _FillClusterTree();
        int  _NodeError(std::string nodeName, int retCode);
        void _CreateNodes(PHCompositeNode *topNode);
};

#endif
