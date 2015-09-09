//-------------------------------------------------------------

#ifndef JetAlgo_h
#define JetAlgo_h

// C++ includes
#include <string>
#include <math.h>
#include <iostream>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TObject.h>
#include <TH1F.h>
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TAxis.h"

// includes for jet's algorithm 
#include "RecoTools.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

using namespace std;

class JetAlgo : public RecoTools
{
public:

  //! constructor
  JetAlgo();

  //! destructor
  virtual ~JetAlgo();

  int setInputTree(string, string, vector<vector<string>>, string);

  void setVerbose(bool);
  void setAlgorithm(string);
  void setRadius(float, double, double);
  void setJetOfJet(float, float, double, double);
  void setWideJet(float, float, double, double);
  
  float jRadius(vector<fastjet::PseudoJet>);
  float jetOfJet(vector<fastjet::PseudoJet>);
  float wideJet(vector<fastjet::PseudoJet>);

  //! loop over events (run algorithm)
  int loop(string, unsigned short, float, double, double); 

  //! analysis
  void integraDraw(TH1F*, TCanvas*, double, double, double, double, double, float);
  void integra(TH1F*, float, float, double*, double*, double*);
  void integra(vector<double>*, float, double*, double*); 
  void meanWeight(TH1F*, double*);
  int analysis(string, unsigned short, float, double, double, float);

private:

  TFile *file;
  TTree *tree;
  int *nParticles;
 
  string mode;
  
  int AllParticles(vector<fastjet::PseudoJet>*, unsigned long);
  int AllGenParticles(vector<fastjet::PseudoJet>*, unsigned long); 

  unsigned short nType;
  unsigned short *nParam;
  TLeaf ***leafs;

  string algorithm;
  bool rFlag;
  float rRadius; 
  double rPT; double rEta;

  bool jjFlag;
  float jjRadius; float jjSteps; float jjLastRadius;
  double jjPT; double jjEta;

  bool wjFlag;
  float wjRadius; float wjLastRadius;
  float wjPT; float wjEta;
 
  bool verbose;
};

#endif
