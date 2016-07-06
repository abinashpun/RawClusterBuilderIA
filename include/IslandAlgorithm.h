#include <vector>
#include <map>
#include "RTHelper.h"
#include <g4cemc/RawTowerContainer.h>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

namespace IslandAlgorithm {
    // 1. The island algorithm starts by a search for seeds. Seeds are defined as 
    // crystals with an energy above a certain threshold on transverse energy. 
    vector<RTHelper> GetSeedTowers(RawTowerContainer* _towers, float _threshold=0.) {

        vector<RTHelper> seedTowers;
        BOOST_FOREACH(RawTowerContainer::Map& towerMap, _towers) {
            if (towerMap.second->get_energy() > threshold) {
            }
        }
    
    }


}
