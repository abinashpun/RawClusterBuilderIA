#include <vector>
#include <map>
#include "RTHelper.h"
#include <g4cemc/RawTower.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTowerGeomContainer.h>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#ifndef __MYRAWCLUSTERBUILDER_H__
typedef RawTowerContainer       RTContainer;
typedef RawTowerGeomContainer   RTGeomContainer;
#endif

typedef std::pair<const unsigned int, RawTower*> RawTowerPair;

namespace IslandAlgorithm {
    using std::vector;
    using std::cout;
    using std::endl;

    class GridFinder;
    bool                comp(RTHelper, RTHelper);
    vector<RTHelper>    GetSeedTowers(RTContainer*, RTGeomContainer*, float);
    void                _RemoveAdjacentSeeds(vector<RTHelper> seedTowers);

    // compare function for std::sort.
    bool comp(RTHelper tower1, RTHelper tower2) { 
        return tower1.get_energy() > tower2.get_energy(); 
    }

    // 1. The island algorithm starts by a search for seeds. Seeds are defined as 
    // crystals with an energy above a certain threshold on transverse energy. 
    vector<RTHelper> GetSeedTowers(RawTowerContainer* _towers, 
                                   RawTowerGeomContainer* _towerGeom, 
                                   float _threshold=0.) {


        vector<RTHelper> seedTowers;
        foreach (RawTowerPair& towerMap, _towers->getTowers()) {
            if (towerMap.second->get_energy() > _threshold) {
                RTHelper rtHelper(towerMap.second);
                rtHelper.set_id(towerMap.first);
                rtHelper.setCenter(_towerGeom);
                seedTowers.push_back(rtHelper);
            }
        }

        seedTowers = _RemoveAdjacentSeeds(seedTowers);

        // Order seeds in decreasing energy.
        std::sort(seedTowers.begin(), seedTowers.end(), comp);
        return seedTowers;
    }

    vector<RTHelper> _RemoveAdjacentSeeds(vector<RTHelper> seedTowers) {
        vector<RTHelper> ret;
        GridFinder gridFinder(seedTowers);
        foreach (RTHelper& tower, seedTowers) {
            bool hasHighestEnergy = true;
            foreach (RTHelper& adjTower, gridFinder.AdjacentTo(tower)) {
                if (adjTower.get_energy() > tower.get_energy()) {
                    hasHighestEnergy = false;
                }
            }
            if (hasHighestEnergy) ret.push_back(tower);
        }
        return ret;
    }

    // helps quickly find adjacent towers.
    class GridFinder {
        public:
            GridFinder(vector<RTHelper>& v);
            vector<RTHelper> AdjacentTo(float binEta, float binPhi);
            vector<RTHelper> AdjacentTo(RTHelper tower);
            RTHelper TowerAt(float binEta, float binPhi);
            vector<RTHelper> getTowers() const { return _towerVector; }
            int IndexOf(float binEta, float binPhi);
        private:
            std::pair<bool,int>** boolMap;
            vector<RTHelper> _towerVector;
            int nBinsPhi, nBinsEta;

    };

    GridFinder::GridFinder(vector<RTHelper>& v) { 
        // Store private variables.
        _towerVector = v; 
        nBinsPhi = RTHelper::get_maxphibin();
        nBinsEta = RTHelper::get_maxetabin();
        std::pair arrBool[nBinsEta][nBinsPhi];
        for (int i=0; i < nBinsEta; i++) {
            for (int j=0; j < nBinsPhi; j++) {
                arrBool[i][j] = std::make_pair(false, -1);
            }
        }
        foreach (RTHelper& tower, _towerVector) {
            arrBool[tower.get_bineta()-1, tower.get_binphi()-1] = std::make_pair(true, -2);
        }

        int towerIndex = 0;
        for (int i=0;i<nBinsEta;i++) {
            for (int j=0;j<nBinsEta;j++) {
                if (arrBool[i][j].first == true) {
                    arrBool[i][j] = std::make_pair(true, towerIndex++);
                }
            }
        }
        boolGrid = arrBool;
    };

    vector<RTHelper> GridFinder::AdjacentTo(float binEta, float binPhi) {
        vector<RTHelper> adjacentTowers;
        int phiRight = (binPhi != nBins) ? binPhi + 1 : 1;
        int phiLeft  = (binPhi != 1)     ? binPhi - 1 : nBins;

        if (boolGrid[binEta-1][phiLeft-1]) adjacentTowers.push_back(TowerAt(binEta, phiLeft));
        if (boolGrid[binEta-1][phiRight-1])adjacentTowers.push_back(TowerAt(binEta, phiRight));
        if (binEta != 1) {
            if (boolGrid[binEta-2][phiLeft-1]) adjacentTowers.push_back(TowerAt(binEta - 1, phiLeft));
            if (boolGrid[binEta-2][binPhi -1]) adjacentTowers.push_back(TowerAt(binEta - 1, binPhi));
            if (boolGrid[binEta-2][phiRight-1]) adjacentTowers.push_back(TowerAt(binEta - 1, phiRight));
        } 
        if (binEta != nBinsEta) {
            if (boolGrid[binEta][phiLeft-1]) adjacentTowers.push_back(TowerAt(binEta + 1, phiLeft));
            if (boolGrid[binEta][binPhi-1]) adjacentTowers.push_back(TowerAt(binEta + 1, binPhi));
            if (boolGrid[binEta][phiRight-1]) adjacentTowers.push_back(TowerAt(binEta + 1, phiRight));
        }
        return adjacentTowers;
    }

    vector<RTHelper> GridFinder::AdjacentTo(RTHelper tower) {
        return AdjacentTo(tower.get_bineta(), tower.get_binphi());
    }

    RTHelper GridFinder::TowerAt(float binEta, float binPhi) {
        return _towerVector[IndexOf(binEta, binPhi)];
    }

    int GridFinder::IndexOf(float binEta, float binPhi) {
        // Note: I'm assuming bins start at one.
        return (binEta - 1) * nBinsPhi + binPhi;
    }

}
