#include <vector>
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
    using std::list;
    using std::vector;
    using std::cout;
    using std::endl;

    // compare function for std::sort.
    bool comp(RTHelper tower1, RTHelper tower2) { 
        return tower1.get_energy() < tower2.get_energy(); 
    }

    void PrintSeeds(list<RTHelper>& seeds) {
        foreach (RTHelper& seed, seeds) {
            cout << "seed (energy, eta, phi) = (" 
                 << seed.get_energy() << ", "
                 << seed.getEtaCenter() << ", "
                 << seed.getPhiCenter() << ")" 
                 << endl;
        }
    }

    // 1. The island algorithm starts by a search for seeds. Seeds are defined as 
    // crystals with an energy above a certain threshold on transverse energy. 
    list<RTHelper> GetSeedTowers(RTContainer* _towers, RTGeomContainer* _towerGeom, float _threshold=0.) {
        // Collect all towers above threshold.
        list<RTHelper> seedTowers;
        foreach (RawTowerPair& towerMap, _towers->getTowers()) {
            if (towerMap.second->get_energy() > _threshold) {
                RTHelper rtHelper(towerMap.second);
                rtHelper.set_id(towerMap.first);
                rtHelper.setCenter(_towerGeom);
                seedTowers.push_back(rtHelper);
            }
        }
        // Find towers with higher-energy adjacent towers.
        list<RTHelper> toRemove;
        seedTowers.sort(comp);
        foreach (RTHelper& tower1, seedTowers) {
            foreach (RTHelper& tower2, seedTowers) {
                if (tower1.is_adjacent(tower2) && tower1.get_energy() < tower2.get_energy()) {
                    toRemove.push_back(tower1);
                }
            }
        }
        // Remove those seeds that are adjacent to higher energy ones.
        foreach (RTHelper& badSeed, toRemove) {
            seedTowers.remove(badSeed);
        }
        PrintSeeds(seedTowers);
        // todo: re-sort so that high energy first.
        return seedTowers;
    }


    void ClusterTowers(list<RTHelper> seedTowers, TowerMap clusteredTowers, RTContainer* _towers) {
        int ClusterID = 0;
        foreach (RTHelper& seed, seedTowers) {
            //PhiCluster(seed, clusteredTowers, _towers);
            /*
            binEta = seed.get_bineta();
            binPhi = seed.get_binphi();
            RawTower* tower;
            while ( 
            ClusterID++;
            */
            while (haven't found a hole or rise in energy):
                etaRight++;
                searchNorth();
                searchSouth();
            while (haven't found a hole or rise in energy):
                etaLeft++;
                searchNorth();
                searchSouth();
        }
    }

    /*
    void PhiCluster(RTHelper& seed, TowerMap clusteredTowers, RTContainer* _towers) {
    }

        groups.insert(std::make_pair(component[i], hits[i]));
        */

}
