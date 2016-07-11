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
        foreach (RawTowerPair& towerMap, _towers->getTowers()) {
            if (towerMap.second->get_energy() > _threshold) {
                RTHelper rtHelper(towerMap.second);
                rtHelper.setID(towerMap.first);
                rtHelper.setCenter(_towerGeom);
                seedTowers.push_back(rtHelper);
            }
        }

        // Find towers with higher-energy adjacent towers.
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

    int movePhi(std::string direction, int& currBinPhi);
    void _SearchEast(RTHelper seed, RTContainer* _towers, TowerMap& clusteredTowers, RTGeomContainer* _towerGeom);
    void _Search(std::string direction, RTHelper seed, RTContainer* _towers, TowerMap& clusteredTowers, RTGeomContainer* _towerGeom);
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
                                RTContainer* _towers, 
                                RTGeomContainer* _towerGeom) {

        //int clusterID = 0;
        TowerMap clusteredTowers;
        /*
        foreach (RTHelper& seed, seedTowers) {
            // Begin by inserting the seed tower, which defines a cluster.
            clusteredTowers.insert(std::make_pair(clusterID, seed));
            _Search("north", seed, _towers, clusteredTowers, _towerGeom);
            _Search("south", seed, _towers, clusteredTowers, _towerGeom);
            _SearchEast(seed, _towers, clusteredTowers, _towerGeom);
            _SearchWest(seed, _towers, clusteredTowers, _towerGeom);
            clusterID++;
        }
        */
        return clusteredTowers;
    }


    void _SearchEast(RTHelper seed, RTContainer* _towers, TowerMap& clusteredTowers, RTGeomContainer* _towerGeom) {
        /*
        // Get next tower info to decide if we should add it.
        int currBinEta = seed.getBinEta();
        if (currBinEta != RTHelper::getMaxEtaBin()) {
            int nextBinEta = currBinEta++;
            RawTower* nextTower = _towers->getTower(nextBinEta, currBinPhi);
            if (nextTower) {
                RTHelper rtHelper(nextTower);
                rtHelper.setCenter(_towerGeom);
                _Search("north", rtHelper, _towers, clusteredTowers, _towerGeom);
                _Search("south", rtHelper, _towers, clusteredTowers, _towerGeom);
                _SearchEast(rtHelper, _towers, clusteredTowers, _towerGeom);
            } else {
                return;
            }
        }
        */
    }

    void _Search(std::string direction, RTHelper seed, RTContainer* _towers, TowerMap& clusteredTowers, RTGeomContainer* _towerGeom) {
        /*
         // Get current tower info.
        float currentEnergy = seed.getEnergy();
        int currBinPhi = seed.getBinPhi();
        int currBinEta = seed.getBinEta();
        // Get next tower info to decide if we should add it.
        int nextBinPhi = (currBinPhi != RTHelper::getMaxPhiBin()) ? currBinPhi++ : 1;
        RawTower* nextTower = _towers->getTower(currBinEta, nextBinPhi);
        if (!nextTower) return;
        float nextEnergy = nextTower->getEnergy();
        // Keep doing this until energy increase or hole.
        while (currentEnergy > nextEnergy) {
            // Add the nextTower to cluster.
            RTHelper rtHelper(nextTower);
            rtHelper.setCenter(_towerGeom);
            clusteredTowers.insert(std::make_pair(clusterID, rtHelper));
            // Setup to check the tower after that.
            nextTower = _towers->getTower(currBinEta, movePhi(direction, nextBinPhi));
            if (!nextTower) return;
            currentEnergy = nextEnergy;
            nextEnergy = nextTower->getEnergy();
        }
        */
    }

    int movePhi(std::string direction, int& currBinPhi) { 
        if (direction == "north") currBinPhi = (currBinPhi != RTHelper::getMaxPhiBin()) ? currBinPhi++ : 1;
        else if (direction == "south") currBinPhi = (currBinPhi != 1) ? currBinPhi-- :
            RTHelper::getMaxPhiBin();
        return currBinPhi;

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
