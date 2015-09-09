#include "RecoTools.hh"

RecoTools::RecoTools() { 
}

RecoTools::~RecoTools() {
}

fastjet::ClusterSequence RecoTools::JetMaker(vector<fastjet::PseudoJet> particles, fastjet::JetDefinition jet_def) {
  // run the clustering, extract the jets
  return fastjet::ClusterSequence(particles, jet_def);
}

fastjet::PseudoJet RecoTools::ConvertToPseudoJet(TLorentzVector v) {
  return fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E());
}

vector<fastjet::PseudoJet> RecoTools::ConvertToPseudoJet(vector<TLorentzVector> v) {
  vector<fastjet::PseudoJet> pseudojets;
  for(int i=0; i<v.size(); i++) {
    pseudojets.push_back(ConvertToPseudoJet(v[i]));
  }
  return pseudojets;
}

TLorentzVector RecoTools::ConvertTo4Vector(fastjet::PseudoJet v) {
  return TLorentzVector(v.px(), v.py(), v.pz(), v.e());
}

vector<TLorentzVector> RecoTools::ConvertTo4Vector(vector<fastjet::PseudoJet> v) {
  vector<TLorentzVector> fourvec;
  for(int i=0; i<v.size(); i++) {
    fourvec.push_back(ConvertTo4Vector(v[i]));
  }
  return fourvec;
}

vector<fastjet::PseudoJet> RecoTools::SelectByAcceptance(vector<fastjet::PseudoJet> v, double pT, double eta){
  vector<fastjet::PseudoJet> vOUT;
  for(int i = 0; i < v.size(); i++){
    if(v[i].pt() > pT && fabs(v[i].eta()) < eta){
      vOUT.push_back(v[i]);
    }
  }
  return vOUT;
}

vector<TLorentzVector> RecoTools::SelectByAcceptance(vector<TLorentzVector> v, double pT, double eta){
  vector<TLorentzVector> vOUT;
  for(int i = 0; i < v.size(); i++){
    if(v[i].Pt() > pT && fabs(v[i].Eta()) < eta){
      vOUT.push_back(v[i]);
    }
  }
  return vOUT;
}

vector<fastjet::PseudoJet> RecoTools::SortByPt(vector<fastjet::PseudoJet> v){
  vector<fastjet::PseudoJet> sorted;
  vector<pair<double,int> > pT;
  int N = v.size();
  for(int i=0; i<N; i++) {
    double pt = v[i].pt();
    pT.push_back(std::make_pair(pt,i));
  }
  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
  for(int i=0; i<pT.size(); i++) 
    sorted.push_back(v[pT.at(i).second]);
  return sorted;
}

vector<TLorentzVector> RecoTools::SortByPt(vector<TLorentzVector> v){
  vector<TLorentzVector> sorted;
  vector<pair<double,int> > pT;
  int N = v.size();
  for(int i=0; i<N; i++) {
    double pt = v[i].Pt();
    pT.push_back(std::make_pair(pt,i));
  }
  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
  for(int i=0; i<pT.size(); i++) 
    sorted.push_back(v[pT.at(i).second]);
  return sorted;
}

vector<fastjet::PseudoJet> RecoTools::SortByEt(vector<fastjet::PseudoJet> v){
  return ConvertToPseudoJet(SortByEt(ConvertTo4Vector(v)));
}

vector<TLorentzVector> RecoTools::SortByEt(vector<TLorentzVector> v){
  vector<TLorentzVector> sorted;
  vector<pair<double,int> > pT;
  int N = v.size();
  for(int i=0; i<N; i++) {
    double pt = v[i].Et(); // No need to change variable names.
    pT.push_back(std::make_pair(pt,i));
  }
  // sort from smallest to larger, then reverse.
  std::sort(pT.begin(), pT.end());
  std::reverse(pT.begin(), pT.end());
  for(int i=0; i<pT.size(); i++) 
    sorted.push_back(v[pT.at(i).second]);
  return sorted;
}

bool RecoTools::IsCharged(int pdgid){
  int id = abs(pdgid);
  if(id == 12 || // nu_e
     id == 14 || // nu_mu
     id == 16 || // nu_tau
     id == 22 || // photon
     id ==  130 || // KL
     id ==  2112 || // n
     id == 1000021 || // gluino
     id == 1000022 || // chi0_1
     id == 1000039 || // gravitino
     id ==  39 || // graviton
     id ==  111 // pi0
     ) return false;
  // the others should be charged or unstable
  return true;
}

fastjet::PseudoJet RecoTools::AngularMatching(vector<fastjet::PseudoJet> j, fastjet::PseudoJet V) {
  double DRmin = 99999999999.;
  int iGOOD = -99;
  for(int i=0; i< j.size(); i++) {
    if(j[i].delta_R(V)<DRmin) {
      DRmin = j[i].delta_R(V);
      iGOOD = i;
    }
  }
  return j[iGOOD];
}

TLorentzVector RecoTools::AngularMatching(vector<TLorentzVector> j, TLorentzVector V) {
  return ConvertTo4Vector(AngularMatching(ConvertToPseudoJet(j), ConvertToPseudoJet(V)));
}

double RecoTools::CalcMR(TLorentzVector ja,  TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();
  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
  double mygamma = 1./sqrt(1.-mybeta*mybeta);
  //gamma times MRstar                                                                                                                                                         \
  temp *= mygamma;
  return temp;
}

double RecoTools::CalcMRT(TLorentzVector ja,  TLorentzVector jb, TLorentzVector met){
  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Vect().Dot(ja.Vect()+jb.Vect());
  temp /= 2.;
  temp = sqrt(temp);
  return temp;
}

double RecoTools::CalcR(TLorentzVector j1,  TLorentzVector j2, TLorentzVector MET){
  return CalcMRT(j1, j2, MET)/CalcMR(j1,j2);
}

double RecoTools::CalcRsq(TLorentzVector j1,  TLorentzVector j2, TLorentzVector MET){
  double R = CalcR(j1,j2,MET);
  return R*R;
}

int RecoTools::HighestPt(vector<fastjet::PseudoJet> jets, int otherJet) {
  double maxPt = 0.;
  int higestPt = -99;
  for(int i=0; i<jets.size(); i++) {
    if(i == otherJet) continue;
    if(jets[i].pt()>maxPt) {
      maxPt = jets[i].pt();
      higestPt = i;
    }
  }
  return higestPt;
}


