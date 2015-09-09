//--------------------------------------------------------------
// Description:  Collection of useful reconstruction algorithms
// Authors:      Maurizio Pierini, CERN
//--------------------------------------------------------------

/// The RecoTools class provides generic reconstruction algorithms
#ifndef RecoTools_h
#define RecoTools_h

// std includes
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TLorentzVector.h>

// FASTJET includes     
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"

using namespace std;

class RecoTools {

public:

  /// Class Constructor
  RecoTools();
  /// Class Destructor
  virtual ~RecoTools();

protected:

  /// Angular matching between a list of jets and a reference jet
  /// The closests jet to the reference jet is returned.
  fastjet::PseudoJet AngularMatching(vector<fastjet::PseudoJet> j, fastjet::PseudoJet V);

  /// Angular matching between a list of jets and a reference jet
  /// The closests jet to the reference jet is returned.
  TLorentzVector AngularMatching(vector<TLorentzVector> j, TLorentzVector V);

  /// Cluster sequence from a generic list of constituents
  fastjet::ClusterSequence JetMaker(vector<fastjet::PseudoJet> p, fastjet::JetDefinition jet_def);

  /// Convert a vector of TLorentzVectors to a vector of fastjet pseudojets
  vector<fastjet::PseudoJet> ConvertToPseudoJet(vector<TLorentzVector> v);

  /// Convert a TLorentzVectors to a fastjet pseudojets
  fastjet::PseudoJet ConvertToPseudoJet(TLorentzVector v);

  /// Convert a vector of fastjet pseudojets  to a vector of TLorentzVectors
  vector<TLorentzVector> ConvertTo4Vector(vector<fastjet::PseudoJet> v);

  /// Convert a fastjet pseudojets to a TLorentzVectors
  TLorentzVector ConvertTo4Vector(fastjet::PseudoJet v);

  /// Compute the razor variable MR
  double CalcMR(TLorentzVector j1,  TLorentzVector j2);
  
  /// Compute the razor variable MRT
  double CalcMRT(TLorentzVector j1,  TLorentzVector j2, TLorentzVector MET);
  
  /// Compute the razor variable R
  double CalcR(TLorentzVector j1,  TLorentzVector j2, TLorentzVector MET);
  
  /// Compute the razor variable Rsq
  double CalcRsq(TLorentzVector j1,  TLorentzVector j2, TLorentzVector MET);
  
  /// Delta phi in radians in between two angles.
  template <typename T> 
  T DeltaPhi(T phi1, T phi2); /// Delta phi in radians in between two angles.

  /// Delta R in between two pairs (eta,phi).
  template <typename T>
  T DeltaR(T eta1, T phi1, T eta2, T phi2); 

  // Compute HT from a collection of TLorentzVector
  double HT(vector<TLorentzVector> jets);

  // Compute HT from a collection of fastjet::PseudoJet
  double HT(vector<fastjet::PseudoJet> jets);

  /// Sorts a collection of TLorentzVector by Pt
  vector<TLorentzVector> SortByPt(vector<TLorentzVector> v);

  /// Sorts a collection of fastjet::PseudoJet by Pt
  vector<fastjet::PseudoJet> SortByPt(vector<fastjet::PseudoJet> v);
  
  /// Sorts a collection of TLorentzVector by Et.
  vector<TLorentzVector> SortByEt(vector<TLorentzVector> v);

  /// Sorts a collection of fastjet::PseudoJet by Et.
  vector<fastjet::PseudoJet> SortByEt(vector<fastjet::PseudoJet> v);

  /// Apply acceptance requirement to a collection of TLorentzVector
  vector<TLorentzVector> SelectByAcceptance(vector<TLorentzVector> v, double pT=0., double eta=10);

  /// Apply acceptance requirement to a collection of fastjet::PseudoJet
  vector<fastjet::PseudoJet> SelectByAcceptance(vector<fastjet::PseudoJet> v, double pT=0., double eta=10);

  /// Check if a Particle has em charge
  bool IsCharged(int pdgid);

  /// find highest pT Lorentx Vector in a list
  int HighestPt(vector<fastjet::PseudoJet> jets, int otherJet);


};

template <typename T>
T RecoTools::DeltaPhi(T phi1, T phi2) {
  T result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

template <typename T>
T RecoTools::DeltaR(T eta1, T phi1, T eta2, T phi2) {
  T dphi = DeltaPhi(phi1,phi2);
  T result = sqrt((eta1-eta2)*(eta1-eta2)+dphi*dphi);
  return result;
}

#endif
