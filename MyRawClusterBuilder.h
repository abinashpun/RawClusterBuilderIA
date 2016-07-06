#ifndef __MYRAWCLUSTERBUILDER_H__
#define __MYRAWCLUSTERBUILDER_H__

// C++ includes.
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
// Fun4All/PHENIX includes.
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
// Local includes.
#include "g4cemc/RawTower.h"
#include "g4cemc/RawTowerGeomContainer.h"
#include "g4cemc/RawTowerContainer.h"
#include "RawClusterv1.h"
#include "g4cemc/RawClusterContainer.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

// Forward class declarations.
class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class RTHelper;
class ClusterHelper;

typedef RawTowerDefs::keytype               KeyType;
typedef RawTowerContainer::ConstRange       ConstRange;
typedef RawTowerContainer::ConstIterator    ConstIterator;

typedef RawTowerContainer                   RTContainer;
typedef RawTowerGeomContainer               RTGeomContainer;
typedef RawTowerContainer::ConstIterator    RTCItr;
typedef RawTowerContainer::ConstRange       RTCRange;
typedef std::multimap<int, RTHelper>::iterator ClusterItr;

class MyRawClusterBuilder : public SubsysReco {
    public:
        MyRawClusterBuilder(const std::string& name = "MyRawClusterBuilder"); 
        virtual ~MyRawClusterBuilder() {}
        int InitRun(PHCompositeNode *topNode);
        int process_event(PHCompositeNode *topNode);
        int End(PHCompositeNode *topNode);
        void Detector(const std::string &d)       {detector = d;}
        void set_threshold_energy(const float e)  {_min_tower_e = e;}
        void checkenergy(const int i = 1)         {chkenergyconservation = i;}
    private:
        RawClusterContainer*  _clusters;
        float                 _min_tower_e;
        int                   chkenergyconservation;
        string                detector;

        void AssignClusterValues(RawCluster* cluster, int iCluster);
        void CreateNodes(PHCompositeNode *topNode);
        bool CorrectPhi(RawCluster* cluster, RTContainer* towers, RTGeomContainer* towerGeom);
        int  NodeError(string nodeName, int retCode);
        void InsertSeed(vector<RTHelper>&, RTCItr, RTGeomContainer*);
        void PrintCluster(ClusterItr);
};


#endif /* RAWCLUSTERBUILDER_H__ */
