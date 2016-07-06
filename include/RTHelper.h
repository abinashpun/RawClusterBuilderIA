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
        bool is_adjacent(RTHelper &);
        void set_id(const int i)                { id = i; }
        void set_etaCenter(float eta)           { etaCenter = eta; }
        void set_phiCenter(float phi)           { phiCenter = phi; }
        void setCenter(RawTowerGeomContainer*);
        int get_id() const                      { return id; }
        int get_bineta() const                  { return bineta; }
        int get_binphi() const                  { return binphi; }
        float getEtaCenter() const              { return etaCenter; }
        float getPhiCenter() const              { return phiCenter; }
        static void set_maxphibin(const int i)  { maxphibin = i; }
        static int get_maxphibin()              { return maxphibin; }
        static RawTower* GetRawTower(RTHelper, RawTowerContainer*);
    protected:
        static void ExitOnIDMismatch(int id1, int id2);
        static int maxphibin;
        RawTowerDefs::keytype id;
        float   etaCenter, phiCenter;
        int     bineta, binphi;
};
int RTHelper::maxphibin = -10;

void RTHelper::setCenter(RawTowerGeomContainer* towerGeom) {
    etaCenter = towerGeom->get_etacenter(bineta);
    phiCenter = towerGeom->get_phicenter(binphi);
}

RawTower* RTHelper::GetRawTower(RTHelper towerHelper, RawTowerContainer* towers) {
    int iPhi = towerHelper.get_binphi();
    int iEta = towerHelper.get_bineta();

    // Ensure ids match. TODO: not really sure what this means?
    int eid = (int) RawTowerDefs::encode_towerid(towers->getCalorimeterID(), iEta, iPhi);
    if (towerHelper.get_id() != eid) ExitOnIDMismatch(towerHelper.get_id(), eid);

    return towers->getTower(iEta, iPhi);
}

RTHelper::RTHelper(RawTower *rt) : id(-1) {
    bineta = rt->get_bineta();
    binphi = rt->get_binphi();
}

// note: true for diagonally adjacent
bool RTHelper::is_adjacent(RTHelper &tower) {
    if (bineta - 1 <= tower.get_bineta() && tower.get_bineta()<=bineta+1) {
        if(binphi-1<=tower.get_binphi() && tower.get_binphi()<=binphi+1) {
            return true;
        } else if(((tower.get_binphi() == maxphibin-1) && (binphi == 0)) || 
                  ((tower.get_binphi() == 0) && (binphi == maxphibin-1))) {
            return true;
        }
    }
    return false;
}

// Comparison first on bineta if not equal, else on binphi.
bool operator<(const RTHelper& a, const RTHelper& b) {
    if (a.get_bineta() != b.get_bineta()) {
        return a.get_bineta() < b.get_bineta();
    }
    return a.get_binphi() < b.get_binphi();
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
