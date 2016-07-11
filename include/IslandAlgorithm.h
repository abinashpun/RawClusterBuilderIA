#include <map>
#include "RTHelper.h"
#include "PHMakeGroups.h"
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
typedef RawTowerContainer            RTContainer;
typedef RawTowerGeomContainer        RTGeomContainer;
typedef std::multimap<int, RTHelper> TowerMap;
typedef std::pair<const int, RTHelper> TowerPair;
typedef std::pair<const unsigned int, RawTower*> RawTowerPair;
#endif


namespace IslandAlgorithm {
    using std::cout;
    using std::endl;

    // Required forward declarations. Definitions at end of file.
    bool lessEnergy(RTHelper tower1, RTHelper tower2);
    bool moreEnergy(RTHelper tower1, RTHelper tower2);
    void PrintSeeds(std::list<RTHelper>& seeds);

    /* ------------------------------------------------------------------------------------------ *
       1.   The island algorithm starts by a search for seeds. Seeds are defined as                
            crystals with an energy above a certain threshold on transverse energy.                
     * ------------------------------------------------------------------------------------------ */
    std::list<RTHelper> GetSeedTowers(RTContainer* _towers, RTGeomContainer* _towerGeom, float _threshold=0.) {

        // Collect all towers above threshold.
        std::list<RTHelper> seedTowers;
        foreach (RawTowerPair& rawTowerPair, _towers->getTowers()) {
            if (rawTowerPair.second->get_energy() > _threshold) {
                RTHelper rtHelper(rawTowerPair.second);
                rtHelper.setID(rawTowerPair.first);
                rtHelper.setCenter(_towerGeom);
                seedTowers.push_back(rtHelper);
            }
        }

        // Find towers with higher-energy adjacent towers.(aka bad seeds)
        std::set<RTHelper> badSeeds;
        seedTowers.sort(lessEnergy);
        // Todo: this could _definitely_ be optimized. 
        foreach (RTHelper& tower1, seedTowers) {
            foreach (RTHelper& tower2, seedTowers) {
                if (tower1.isAdjacent(tower2) && tower1.getEnergy() < tower2.getEnergy()) {
                    badSeeds.insert(tower1);
                }
            }
        }

        // Remove those seeds that are adjacent to higher energy ones.
        foreach (const RTHelper badSeed, badSeeds) {
            seedTowers.remove(badSeed);
        }

        // Order from hi-to-lo energy.
        seedTowers.sort(moreEnergy);

        PrintSeeds(seedTowers);
        return seedTowers;
    }

    int moveEta(std::string direction, int& currBinEta);
    int movePhi(std::string direction, int& currBinPhi);
    void _SearchEta(std::string direction,       RTHelper           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom);
    void _SearchPhi(std::string direction,       RTHelper           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom);
    /* ------------------------------------------------------------------------------------------ *
       2.   Starting from the most energetic seed, the algorithm collects crystals belonging to 
            a certain cluster. Moves both directions in phi, collecting all towers until it sees
            a rise in energy, or a hole. The algorithm then steps in eta and makes another phi
            search. The eta-steps are stopped when a rise in energy, or a hole, is encountered.
            When in direction in eta is completed, the algorithm goes back to the seed position 
            and works in the other eta direction. All towers are makred as belonging to that one 
            cluster and can't be subsequently used to seed another. 
     * ------------------------------------------------------------------------------------------ */
    TowerMap GetClusteredTowers(std::list<RTHelper> seedTowers, 
                                RTContainer*        _towers, 
                                RTGeomContainer*    _towerGeom) {

        int clusterID = 0;
        TowerMap clusteredTowers;
        foreach (RTHelper& seed, seedTowers) {
            // Begin by inserting the seed tower, which defines a cluster.
            clusteredTowers.insert(std::make_pair(clusterID, seed));
            _SearchPhi("north", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            _SearchPhi("south", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            _SearchEta("west", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            _SearchEta("east", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            clusterID++;
        }
        return clusteredTowers;
    }

    /* -------------------------------------------------------------------------
       _SearchPhi
       ------------------------------------------------------------------------- */
    void _SearchPhi(std::string direction,       RTHelper           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom) {

        // Get current tower info.
        float currEnergy = currentTower.getEnergy();
        int currBinPhi   = currentTower.getBinPhi();
        int currBinEta   = currentTower.getBinEta();

        // Get next tower info to decide if we should add it.
        RawTower* nextTower = _towers->getTower(currBinEta, movePhi(direction, currBinPhi));
        // Keep doing this until energy increase or hole.
        if (nextTower && currEnergy > nextTower->get_energy()) {
            RTHelper rtHelper(nextTower);
            rtHelper.setCenter(_towerGeom);
            clusteredTowers.insert(std::make_pair(clusterID, rtHelper));
            _SearchPhi(direction, rtHelper, clusterID, _towers, clusteredTowers, _towerGeom);
        }
    }

    /* -------------------------------------------------------------------------
       _SearchEta
       ------------------------------------------------------------------------- */
    void _SearchEta(std::string direction,       RTHelper           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom) {
        // Get next tower info to decide if we should add it.
        int currBinEta = currentTower.getBinEta();
        int currBinPhi = currentTower.getBinPhi();
        float currEnergy = currentTower.getEnergy();
        if (moveEta(direction, currBinEta) != -1) {
            RawTower* nextTower = _towers->getTower(currBinEta, currBinPhi);
            if (nextTower && currEnergy > nextTower->get_energy()) {
                RTHelper rtHelper(nextTower);
                rtHelper.setCenter(_towerGeom);
                clusteredTowers.insert(std::make_pair(clusterID, rtHelper));
                _SearchPhi("north", rtHelper, clusterID, _towers, clusteredTowers, _towerGeom);
                _SearchPhi("south", rtHelper, clusterID, _towers, clusteredTowers, _towerGeom);
                _SearchEta(direction, rtHelper, clusterID, _towers, clusteredTowers, _towerGeom);
            }
        }
    }

    int movePhi(std::string direction, int& currBinPhi) { 
        if (direction == "north") currBinPhi = (currBinPhi != RTHelper::getMaxPhiBin()) ? currBinPhi++ : 1;
        else if (direction == "south") currBinPhi = (currBinPhi != 1) ? currBinPhi-- :
            RTHelper::getMaxPhiBin();
        return currBinPhi;
    }

    int moveEta(std::string direction, int& currBinEta) { 
        if (direction == "west") currBinEta = (currBinEta != 1) ? currBinEta-- : -1;
        else if (direction == "east") currBinEta = (currBinEta != RTHelper::getMaxEtaBin()) ? currBinEta++ : -1;
        return currBinEta;
    }

    // A simple comparator that orders RTHelpers in order of INCREASING energy.
    bool lessEnergy(RTHelper tower1, RTHelper tower2) { 
        return tower1.getEnergy() < tower2.getEnergy(); 
    }

    // A simple comparator that orders RTHelpers in order of DECREASING energy.
    bool moreEnergy(RTHelper tower1, RTHelper tower2) {
        return !lessEnergy(tower1, tower2);
    }

    // Essentially a 'ToString' method for a list of seed towers.
    void PrintSeeds(std::list<RTHelper>& seeds) {
        foreach (RTHelper& seed, seeds) {
            cout << "seed (energy, eta, phi) = (" 
                 << seed.getEnergy() << ", "
                 << seed.getEtaCenter() << ", "
                 << seed.getPhiCenter() << ")" 
                 << endl;
        }
    }

}
