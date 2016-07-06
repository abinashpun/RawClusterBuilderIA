#ifndef __TRUTHJETTRIGGER_H__
#define __TRUTHJETTRIGGER_H__


// --- need to check all these includes...
#include <fun4all/SubsysReco.h>
#include <vector>

#include "TMath.h"

#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"

class PHCompositeNode;

class MyRawClusterBuilder: public SubsysReco
{

 public:

  MyRawClusterBuilder();

  int Init(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

 private:

};

#endif // __TRUTHJETTRIGGER_H__
