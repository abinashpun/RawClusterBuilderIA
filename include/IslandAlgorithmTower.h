// Raw Tower Helper Class.
#ifndef __ISLANDALGORITHMTOWER_H__
#define __ISLANDALGORITHMTOWER_H__

// C++ includes.
#include <iostream>
// Local includes.
#include "g4cemc/RawTower.h"
#include "g4cemc/RawTowerGeomContainer.h"
#include "g4cemc/RawTowerContainer.h"


class IslandAlgorithmTower {
    public:
        IslandAlgorithmTower(RawTower *);
        virtual ~IslandAlgorithmTower() {}
        bool isAdjacent(IslandAlgorithmTower &);
        int getID() const                       { return id; }
        int getBinEta() const                   { return bineta; }
        int getBinPhi() const                   { return binphi; }
        float getET() const                     { return energy / std::cosh(etaCenter);}
        float getEnergy() const                 { return energy; }
        float getEtaCenter() const              { return etaCenter; }
        float getPhiCenter() const              { return phiCenter; }
        void setID(const int i)                 { id = i; }
        void setEnergy(float e)                 { energy = e; }
        void setEtaCenter(float eta)            { etaCenter = eta; }
        void setPhiCenter(float phi)            { phiCenter = phi; }
        void setCenter(RawTowerGeomContainer*);
        static int getMaxPhiBin()               { return maxPhiBin; }
        static int getMaxEtaBin()               { return maxEtaBin; }
        static void setMaxPhiBin(const int i)   { maxPhiBin = i; }
        static void setMaxEtaBin(const int i)   { maxEtaBin = i; }
        static RawTower* GetRawTower(IslandAlgorithmTower, RawTowerContainer*);
    protected:
        RawTowerDefs::keytype id;
        int bineta; 
        int binphi;
        float energy;
        float etaCenter; 
        float phiCenter;
        static int maxPhiBin;
        static int maxEtaBin;
        static void ExitOnIDMismatch(int id1, int id2);
};
int IslandAlgorithmTower::maxPhiBin = -10;
int IslandAlgorithmTower::maxEtaBin = -10;

void IslandAlgorithmTower::setCenter(RawTowerGeomContainer* towerGeom) {
    etaCenter = towerGeom->get_etacenter(bineta);
    phiCenter = towerGeom->get_phicenter(binphi);
}

/* ----------------------------------------
   IslandAlgorithmTower(RawTower*) constructor.
   ---------------------------------------- */
IslandAlgorithmTower::IslandAlgorithmTower(RawTower *rt) : id(-1) {
    bineta = rt->get_bineta();
    binphi = rt->get_binphi();
    energy = rt->get_energy();
    id = rt->get_id();
}


// note: true for diagonally adjacent
bool IslandAlgorithmTower::isAdjacent(IslandAlgorithmTower &tower) {
    if (bineta - 1 <= tower.getBinEta() && tower.getBinEta()<=bineta+1) {
        if(binphi-1<=tower.getBinPhi() && tower.getBinPhi()<=binphi+1) {
            return true;
        } else if(((tower.getBinPhi() == maxPhiBin-1) && (binphi == 0)) || 
                ((tower.getBinPhi() == 0) && (binphi == maxPhiBin-1))) {
            return true;
        }
    }
    return false;
}

// Comparison first on bineta if not equal, else on binphi.
bool operator<(const IslandAlgorithmTower& a, const IslandAlgorithmTower& b) {
    if (a.getBinEta() != b.getBinEta()) {
        return a.getBinEta() < b.getBinEta();
    }
    return a.getBinPhi() < b.getBinPhi();
}

bool operator==(const IslandAlgorithmTower& a, const IslandAlgorithmTower& b) {
    return a.getID() == b.getID();
}

RawTower* IslandAlgorithmTower::GetRawTower(IslandAlgorithmTower towerHelper, RawTowerContainer* towers) {
    int iPhi = towerHelper.getBinPhi();
    int iEta = towerHelper.getBinEta();
    // Ensure ids match. TODO: not really sure what this means?
    int eid = (int) RawTowerDefs::encode_towerid(towers->getCalorimeterID(), iEta, iPhi);
    if (towerHelper.getID() != eid) ExitOnIDMismatch(towerHelper.getID(), eid);
    return towers->getTower(iEta, iPhi);
}

void IslandAlgorithmTower::ExitOnIDMismatch(int id1, int id2) {
    std::cout <<__PRETTY_FUNCTION__
        << " - Fatal Error - id mismatch. internal: " 
        << id1
        << ", towercontainer: " 
        << id2
        << std::endl;
    exit(1);
}

#endif
