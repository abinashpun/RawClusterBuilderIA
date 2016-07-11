// Raw Tower Helper Class.
#ifndef RTHELPER_H_
#define RTHELPER_H_
// C++ includes.
#include <iostream>
// Local includes.
#include "g4cemc/RawTower.h"
#include "g4cemc/RawTowerGeomContainer.h"
#include "g4cemc/RawTowerContainer.h"

class RTHelper {
    public:
        RTHelper(RawTower *);
        virtual ~RTHelper() {}
        bool isAdjacent(RTHelper &);
        void setID(const int i)                 { id = i; }
        void setEnergy(float e)                 { energy = e; }
        void setCenter(RawTowerGeomContainer*);
        int getID() const                       { return id; }
        int getBinEta() const                   { return bineta; }
        int getBinPhi() const                   { return binphi; }
        float getEnergy() const                 { return energy; }
        void setEtaCenter(float eta)            { etaCenter = eta; }
        void setPhiCenter(float phi)            { phiCenter = phi; }
        float getEtaCenter() const              { return etaCenter; }
        float getPhiCenter() const              { return phiCenter; }
        static void setMaxPhiBin(const int i)   { maxPhiBin = i; }
        static void setMaxEtaBin(const int i)   { maxEtaBin = i; }
        static int getMaxPhiBin()               { return maxPhiBin; }
        static int getMaxEtaBin()               { return maxEtaBin; }
        static RawTower* GetRawTower(RTHelper, RawTowerContainer*);
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
int RTHelper::maxPhiBin = -10;
int RTHelper::maxEtaBin = -10;

void RTHelper::setCenter(RawTowerGeomContainer* towerGeom) {
    etaCenter = towerGeom->get_etacenter(bineta);
    phiCenter = towerGeom->get_phicenter(binphi);
}

/* ----------------------------------------
   RTHelper(RawTower*) constructor.
   ---------------------------------------- */
RTHelper::RTHelper(RawTower *rt) : id(-1) {
    bineta = rt->get_bineta();
    binphi = rt->get_binphi();
    energy = rt->get_energy();
    id = rt->get_id();
}

// note: true for diagonally adjacent
bool RTHelper::isAdjacent(RTHelper &tower) {
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
bool operator<(const RTHelper& a, const RTHelper& b) {
    if (a.getBinEta() != b.getBinEta()) {
        return a.getBinEta() < b.getBinEta();
    }
    return a.getBinPhi() < b.getBinPhi();
}

bool operator==(const RTHelper& a, const RTHelper& b) {
    return a.getID() == b.getID();
}

RawTower* RTHelper::GetRawTower(RTHelper towerHelper, RawTowerContainer* towers) {
    int iPhi = towerHelper.getBinPhi();
    int iEta = towerHelper.getBinEta();
    // Ensure ids match. TODO: not really sure what this means?
    int eid = (int) RawTowerDefs::encode_towerid(towers->getCalorimeterID(), iEta, iPhi);
    if (towerHelper.getID() != eid) ExitOnIDMismatch(towerHelper.getID(), eid);
    return towers->getTower(iEta, iPhi);
}

void RTHelper::ExitOnIDMismatch(int id1, int id2) {
    cout <<__PRETTY_FUNCTION__
         << " - Fatal Error - id mismatch. internal: " 
         << id1
         << ", towercontainer: " 
         << id2
         << endl;
    exit(1);
}

#endif
