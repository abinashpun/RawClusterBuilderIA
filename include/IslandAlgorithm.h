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
                rtHelper.set_id(towerMap.first);
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
                if (tower1.is_adjacent(tower2) && tower1.get_energy() < tower2.get_energy()) {
                    badSeeds.push_back(tower1);
                }
            }
        }

        // Remove those seeds that are adjacent to higher energy ones.
        foreach (RTHelper& badSeed, badSeeds) {
            seedTowers.remove(badSeed);
        }

        // Order from hi-to-lo energy.
        seedTowers.sort(moreEnergy);

        PrintSeeds(seedTowers);
        return seedTowers;
    }


    /* ------------------------------------------------------------------------------------------ *
       2.   Starting from the most energetic seed, the algorithm collects crystals belonging to 
            a certain cluster. Moves both directions in phi, collecting all towers until it sees
            a rise in energy, or a hole. The algorithm then steps in eta and makes another phi
            search. The eta-steps are stopped when a rise in energy, or a hole, is encountered.
            When in direction in eta is completed, the algorithm goes back to the seed position 
            and works in the other eta direction. All towers are makred as belonging to that one 
            cluster and can't be subsequently used to seed another. 
     * ------------------------------------------------------------------------------------------ */
    TowerMap GetClusteredTowers(std::list<RTHelper> seedTowers, std::list<RTHelper> allTowers) {

        int ClusterID = 0;
        TowerMap clusteredTowers;

        foreach (RTHelper& seed, seedTowers) {

            // Begin by inserting the seed tower, which defines a cluster.
            clusteredTowers.insert(std::make_pair(clusterId, seed));
            float currentEnergy = seed.get_energy();
            float currBinPhi = seed.getBinPhi();
            float currBinEta = seed.getEtaPhi();
            while (currentEnergy > _towers->getTower(currBinEta, ++currBinPhi)->get_energy()) {
            }
            clusterID++;
        }
    }


    // A simple comparator that orders RTHelpers in order of INCREASING energy.
    bool lessEnergy(RTHelper tower1, RTHelper tower2) { 
        return tower1.get_energy() < tower2.get_energy(); 
    }

    // A simple comparator that orders RTHelpers in order of DECREASING energy.
    bool moreEnergy(RTHelper tower1, RTHelper tower2) {
        return !lessEnergy(tower1, tower2);
    }

    // Essentially a 'ToString' method for a list of seed towers.
    void PrintSeeds(std::list<RTHelper>& seeds) {
        foreach (RTHelper& seed, seeds) {
            cout << "seed (energy, eta, phi) = (" 
                 << seed.get_energy() << ", "
                 << seed.getEtaCenter() << ", "
                 << seed.getPhiCenter() << ")" 
                 << endl;
        }
    }

}
