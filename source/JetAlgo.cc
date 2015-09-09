#include "JetAlgo.hh"

/*! \fn JetAlgo    
* \brief   constructtor
*/
JetAlgo::JetAlgo() 
{
  /* initialization variable with default value  */
  rFlag = false;
  jjFlag = false;
  wjFlag = false;

}

/*! \fn ~jetAlgo     
* \brief  destructor
*/
JetAlgo::~JetAlgo() 
{
  delete[] nParticles;
  delete[] nParam;
  for(unsigned short i = 0; i < nType; i++)
    delete[] leafs[i];
  
  tree->Delete();
  file->Close();
  file->Delete();
}

/*! \fn setTree     
* \brief Function for upload TTree with data for TLorenzVector.
*
* \param  string inputFile - path to the file
* \param  string inputTree - name of tree in ROOT file
* \param  vector<vector<string>> namesLeafs - array with names Tleaf
* type namesLeafs: namesLeafs[optional][optional]
* Example:
*  leafs[0][0] = "ChargedHadron.PT";
*  leafs[0][1] = "ChargedHadron.Eta";
*  leafs[0][2] = "ChargedHadron.Phi";
*
*  leafs[1][0] = "Electron.PT";
*  leafs[1][1] = "Electron.Eta";
*  leafs[1][2] = "Electron.Phi";
*  leafs[1][3] = "Electron.M";
* 
* \param  string simulationMode = "Rec" or "Gen"
*  name of mode for data simulation: 
*                                   Rec - reconstructed
*                                   Gen - generation
 * \return int - check tag 
*/
int JetAlgo::setInputTree(string inputFile, string inputTree,
        vector<vector<string>> namesLeafs, string simulationMode) 
{
  mode = simulationMode;
  nType = namesLeafs.size(); // nType = numbersType, get numbersType particles.
                            
  nParticles = new int[nType]; // numbersParticles
  nParam = new unsigned short[nType]; // nParam = numbersParameters.
    if(!nParam) return 0;
  leafs = new TLeaf** [nType]; // array for leafs.
    if(!leafs) return 0;

  /* initialization TTree */
  file = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFile.c_str());
  if(!file)
    file = new TFile(inputFile.c_str());

  tree = (TTree*)file->Get(inputTree.c_str());
  if(!tree) { 
    cout << "Can't open input tree or file" << endl;
    return 0;
  }

  /* initialization array TLeaf */
  for(unsigned short i = 0; i < nType; i++) {

    nParam[i] = namesLeafs[i].size();
    leafs[i] = new TLeaf* [nParam[i]];
      if(!leafs[i]) return 0;
    for(unsigned short j = 0; j < nParam[i]; j++) {
      leafs[i][j] = tree->GetLeaf(namesLeafs[i][j].c_str());
        if(!leafs[i][j]) {
          cout << "Can't open input leaf:" << namesLeafs[i][j].c_str() << endl;
          return 0;
        }
    }
  }
  return 1;
}

/*! \fn AllParticles (not universal)      
* \brief Function for get data of Event form TTree(setTree).
                          unsigned long nEntrie
* \param  vector<fastjet::PseudoJet> *particles - return data  
* \param  unsigned long nEntrie - number Event for read data from TTree
* 
* \return int - check tag 
*/
int JetAlgo::AllParticles(vector<fastjet::PseudoJet> *particles, unsigned long nEntrie) 
{

  TLorentzVector thisP; // thisParticle

  for(unsigned short i = 0; i < nType; i++) {

    Float_t **leafValue = new Float_t* [nParam[i]]; // dynamic array for data
      if(!leafValue) return 0;

    for(unsigned short j = 0; j < nParam[i]; j++) {
      leafs[i][j]->GetBranch()->GetEntry(nEntrie);

      nParticles[i] = leafs[i][j]->GetLen(); // work only for same Len 
      leafValue[j] = new Float_t[nParticles[i]];
        if(!leafValue[j]) return 0;

      for(unsigned long n = 0; n < nParticles[i]; n++)
        leafValue[j][n] = (Float_t)leafs[i][j]->GetValue(n); // get data

    }
    unsigned long n = nParticles[i];
    for(unsigned long j = 0; j < n; j++) {
      if(nParam[i] == 3) { 
        thisP.SetPtEtaPhiM(leafValue[0][j], leafValue[1][j],
                           leafValue[2][j], 0.);
        particles->push_back(ConvertToPseudoJet(thisP));
      }
    }
    for(unsigned short j = 0; j < nParam[i]; j++)
	    delete[] leafValue[j];
  }
  return 1;
}

/*! \fn AllParticles (not universal)     
* \brief Function for get data of Event form TTree(setTree).
                          unsigned long nEntrie
* \param  vector<fastjet::PseudoJet> *particles - return data  
* \param  unsigned long nEntrie - number Event for read data from TTree
* 
* \return int - check tag 
*/
int JetAlgo::AllGenParticles(vector<fastjet::PseudoJet> *particles, unsigned long nEntrie) 
{

  TLorentzVector thisP; // thisParticle
  for(unsigned short i = 0; i < nType; i++) {

    Float_t **leafValue = new Float_t* [nParam[i]]; // work only for same Len 

      if(!leafValue) return 0;

    for(unsigned short j = 0; j < nParam[i]; j++) {
      leafs[i][j]->GetBranch()->GetEntry(nEntrie);

      nParticles[i] = leafs[i][j]->GetLen(); // обработать не совпадения размерностей
      leafValue[j] = new Float_t[nParticles[i]];
        if(!leafValue[j]) return 0;

      for(unsigned long n = 0; n < nParticles[i]; n++)
        leafValue[j][n] = (Float_t)leafs[i][j]->GetValue(n); // get data

    }
    
    unsigned long n = nParticles[i];
    for(unsigned long j = 0; j < n; j++) {
      // cut is in mm
      if(nParam[i] == 6) {
	      
        // px, py, pz, E, Status, PID
	      if(leafValue[4][j] == 1) {
	        if(fabs(leafValue[5][j]) == 12 ||
	           fabs(leafValue[5][j]) == 14 ||
	           fabs(leafValue[5][j]) == 16) continue;
          thisP.SetPxPyPzE(leafValue[0][j], leafValue[1][j],
		                       leafValue[2][j], leafValue[3][j]);
	        particles->push_back(ConvertToPseudoJet(thisP));
        }
      }
    }
    for(unsigned short j = 0; j < nParam[i]; j++)
	    delete[] leafValue[j];
  }
  return 1;
}

/*! \fn setVerbose    
* \brief  Function for set verbose param 
* 
* \param bool verbose
*/
void JetAlgo::setVerbose(bool vb)
{
  verbose  = vb;
}

/*! \fn setAlgorithm    
* \brief  Function for set basic algorithm "antiKt"; "kt"; "cam" 
* 
* \param string algorithm
*                         "antiKt"
*                         "kt"
*                         "cam"
*/
void JetAlgo::setAlgorithm(string inputAlgorithm)
{
  algorithm = inputAlgorithm;
}

/*! \fn setRadius    
* \brief  Function for set parameters for jRadius Algorithm
* 
* \param float radius - radius for calculate 
* \param int PT 
* \param int Eta 
*/
void JetAlgo::setRadius(float radius, double PT, double Eta) 
{
  rFlag = true;
  rRadius = radius;
  rPT = PT; rEta = Eta; 
}

/*! \fn setJetOfJet    
* \brief  Function for set parameters for jRadius Algorithm
* 
* \param float radius - first radius 
* \param float lastRadius - last radius (result radius) 
* \param int PT 
* \param int Eta 
*/
void JetAlgo::setJetOfJet(float radius, float lastRadius,
                          double PT, double Eta) 
{
  jjFlag = true;
  jjRadius = radius; jjLastRadius = lastRadius;
  jjPT = PT; jjEta = Eta;
}

/*! \fn setWideJet    
* \brief  Function for set parameters for diJet Algorithm
* 
* \param float radius - first radius 
* \param float lastRadius - last radius (result radius) 
* \param int PT 
* \param int Eta 
*/
void JetAlgo::setWideJet(float radius, float lastRadius,
                         double PT, double Eta) 
{
  wjFlag = true;
  wjRadius = radius; wjLastRadius = lastRadius;
  wjPT = PT; wjEta = Eta;
}

/*! \fn jRadius    
* \brief  Function for jet algorithm (anntiKt, Kt, cambrige)
* 
* \param vector<fastjet::PseudoJet> particles) - input data 
*
* \return float - invariant mass for two biggest jets
*/
float JetAlgo::jRadius(vector<fastjet::PseudoJet> particles)
{
  float massJet = -99.;
  float radius = rRadius; 

  fastjet::JetDefinition jetDef;
  vector<fastjet::PseudoJet> resultJets;

  if(strcmp ("antiKt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
  if(strcmp ("kt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::kt_algorithm, radius);
  if(strcmp ("cam", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::cambridge_algorithm, radius);

  // cluster jets
  fastjet::ClusterSequence jetCS(particles, jetDef);
  resultJets = jetCS.inclusive_jets();
  resultJets = SelectByAcceptance(fastjet::sorted_by_pt(resultJets), rPT, rEta);
  if(resultJets.size() > 1)
    massJet  = float((resultJets[0]+resultJets[1]).m());
  
  return massJet;
}

/*! \fn jetOfJet    
* \brief  Function for jetOfJet.
* 
* \param vector<fastjet::PseudoJet> particles) - input data 
* 
* \return float - invariant mass for two biggest jets
*/
float JetAlgo::jetOfJet(vector<fastjet::PseudoJet> particles) 
{
  float massJet = -99.;
  float radius = jjRadius; 

  fastjet::JetDefinition jetDef;
  vector<fastjet::PseudoJet> resultJets;
  fastjet::ClusterSequence jetCSLast;

  if(strcmp ("antiKt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
  if(strcmp ("kt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::kt_algorithm, radius);
  if(strcmp ("cam", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::cambridge_algorithm, radius);

  fastjet::ClusterSequence jetCS(particles, jetDef);
  resultJets = jetCS.inclusive_jets();
  resultJets = SelectByAcceptance(fastjet::sorted_by_pt(resultJets), jjPT, jjEta);

  /* JetOfJet */
  if(resultJets.size() > 1) {
      radius = jjLastRadius;

  if(strcmp ("antiKt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
  if(strcmp ("kt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::kt_algorithm, radius);
  if(strcmp ("cam", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::cambridge_algorithm, radius);
      jetCSLast = JetMaker(resultJets, jetDef);
      resultJets = SelectByAcceptance(fastjet::sorted_by_pt(jetCSLast.inclusive_jets()), jjPT, jjEta);
    }
    if(resultJets.size() > 1)
    massJet = float((resultJets[0]+resultJets[1]).m());
  return massJet;
}

/*! \fn diJet    
* \brief  Function for diJet.
* 
* \param vector<fastjet::PseudoJet> particles) - input data 
*
* \return float - invariant mass for two biggest jets
*/
float JetAlgo::wideJet(vector<fastjet::PseudoJet> particles)
{
  float massJet = -99.;
  float radius = wjRadius; 
  float lastRadius = wjLastRadius;
  
  fastjet::JetDefinition jetDef;
  vector<fastjet::PseudoJet> resultJet;
  fastjet::AreaDefinition area;
  fastjet::PseudoJet WideJet[2];
  bool inWideJet[2]; inWideJet[0] = false; inWideJet[1] = false;

  if(strcmp ("antiKt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
  if(strcmp ("kt", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::kt_algorithm, radius);
  if(strcmp ("cam", algorithm.c_str()) == 0)
    jetDef = fastjet::JetDefinition(fastjet::cambridge_algorithm, radius);

  fastjet::ClusterSequence jetCS(particles, jetDef);
  resultJet = jetCS.inclusive_jets();
  
  /* jet with default radius */
  resultJet = SelectByAcceptance(fastjet::sorted_by_pt(resultJet), 
                                     wjPT, wjEta);
    if(resultJet.size() > 1) {

        radius = wjLastRadius;
        WideJet[0] = resultJet[0];
        WideJet[1] = resultJet[1];

        for(unsigned short j = 2; j < resultJet.size(); j++) {
          inWideJet[0] = false;
          inWideJet[1] = false;
          if(resultJet[j].pt() < wjPT)
              continue;
          if(fabs(resultJet[j].eta()) > wjEta)
              continue;
          if(resultJet[0].delta_R(resultJet[j]) < lastRadius)
              inWideJet[0] = true;
          if(resultJet[1].delta_R(resultJet[j]) < lastRadius)
              inWideJet[1] = true;
          if(inWideJet[0] && inWideJet[1]) {
            if(resultJet[0].delta_R(resultJet[j]) < resultJet[1].delta_R(resultJet[j]))
                inWideJet[1] = false;
            else
                inWideJet[0] = false;
          }
          if(inWideJet[0])
              WideJet[0] = WideJet[0] + resultJet[j];
          if(inWideJet[1])
              WideJet[1] = WideJet[1] + resultJet[j];
        }
        massJet = float((WideJet[0]+WideJet[1]).m());
    }
    return massJet;
}

/*! \fn invMass    
* \brief  Function for read input Tree and get invariant mass with algorithm.
* 
* \param TH1F *histResult - return data (invariant mass) in hist 
*
*/
int JetAlgo::loop(string outFileName, unsigned short numberSteps,
                   float stepR, double stepPT, double stepEta) 
{
  vector<fastjet::PseudoJet> pfcands;
  
  float trRadius = rRadius;
  float trPT = rPT;
  float trEta = rEta;
  
  float tjjLastRadius = jjLastRadius;
  float tjjPT = jjPT;
  float tjjEta = jjEta;

  float twjLastRadius = wjLastRadius;
  float twjPT = wjPT;
  float twjEta = wjEta;
      
  TFile *outFile = new TFile(outFileName.c_str(), "UPDATE");
  if(!outFile) {
    cout << "Can not create output file" << endl;
    return 0;
  }

  TTree* outTree = new TTree(mode.c_str(), mode.c_str());
  if(!outTree) {
    cout << "Can not create output tree" << endl;
    return 0;
  }

  float massJRadius[numberSteps];
  float massJetOfJet[numberSteps];
  float massWideJet[numberSteps];
  for(unsigned short i = 0; i < numberSteps; i++) {
    if(rFlag == true) {
      outTree->Branch(Form("%s_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), rRadius, rPT, rEta),
                      &massJRadius[i],
                      Form("%s_%.2fR_%.2fPT_%.2fEta/F", algorithm.c_str(), rRadius, rPT, rEta));
      rRadius = rRadius + stepR;
      rPT = rPT + stepPT;
      rEta = rEta + stepEta;
    }
    if(jjFlag == true) {
      outTree->Branch(Form("jetOfjet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), jjRadius, jjLastRadius, jjPT, jjEta),
                      &massJetOfJet[i],
                      Form("jetOfJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta/F", algorithm.c_str(), jjRadius, jjLastRadius, jjPT, jjEta));
      jjLastRadius = jjLastRadius + stepR;
      jjPT = jjPT + stepPT;
      jjEta = jjEta + stepEta;
    }
    if(wjFlag == true) {
      outTree->Branch(Form("wideJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), wjRadius, wjLastRadius, wjPT, wjEta),
                      &massWideJet[i],
                      Form("wideJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta/F", algorithm.c_str(), wjRadius, wjLastRadius, wjPT, wjEta));
      wjLastRadius = wjLastRadius + stepR;
      wjPT = wjPT + stepPT;
      wjEta = wjEta + stepEta;
    }
  }
    rRadius = trRadius; jjLastRadius = tjjLastRadius; wjLastRadius = twjLastRadius;
    rPT = trPT; jjPT = tjjPT; wjPT = twjPT;
    rEta = trEta; jjEta = tjjEta; wjEta = twjEta;
  
  // loop over entries
  Long64_t nentries = tree->GetEntries();
  
  if(verbose == true)
    cout <<  "Number of entries = " << nentries << endl;
  
  for (Long64_t jEntry = 0; jEntry < nentries; jEntry++) {
    if(nentries < 0) {
      cout << "Number of entries  = " << nentries 
           << "It's not correct" << endl;
      return 0;
    }
    
  if(verbose == true) {
    if(jEntry%50 == 0)
      std::cout << ">>> Processing event # " << jEntry << std::endl;
  }

    // get event
  if(strcmp ("Rec", mode.c_str()) == 0)
    AllParticles(&pfcands, jEntry); 
  if(strcmp ("Gen", mode.c_str()) == 0)
    AllGenParticles(&pfcands, jEntry); 
 
    for(unsigned short i = 0; i < numberSteps; i++) {
      if(rFlag == true) {
        rRadius = rRadius + stepR;
        rPT = rPT + stepPT;
        rEta = rEta + stepEta;
        massJRadius[i] = jRadius(pfcands);
      }
      if(jjFlag == true) {
        jjLastRadius = jjLastRadius + stepR;
        jjPT = jjPT + stepPT;
        jjEta = jjEta + stepEta;
        massJetOfJet[i] = jetOfJet(pfcands);
      }
      if(wjFlag == true) {
        wjLastRadius = wjLastRadius + stepR;
        wjPT = wjPT + stepPT;
        wjEta = wjEta + stepEta;
        massWideJet[i] = wideJet(pfcands);
      }
    } 
    rRadius = trRadius; jjLastRadius = tjjLastRadius; wjLastRadius = twjLastRadius;
    rPT = trPT; jjPT = tjjPT; wjPT = twjPT;
    rEta = trEta; jjEta = tjjEta; wjEta = twjEta;

    outTree->Fill();    
    pfcands.clear();
  }
  outTree->Write();    
  outTree->Delete();
  
  outFile->Close();
  outFile->Delete();
  return 1;
}

//______________________________________________________________________________
// ANALYSIS PART

/*! \fn integraDraw    
* \brief  Function draw result from integra function.
* 
* \param TH1F *hist - data hist 
* \param double first - first value of interval 
* \param double last - last value of interval 
* \param double mean - mean value 
* \param double xmin - min value for xAxis
* \param double bin - bin size 
* \param double percent - percentage of the number of events for integral   
* \param TH1F *histResult - return data (invariant mass) in hist 

* \param (return) TCanvas *canvas - merge hist and result integral 
*/
void JetAlgo::integraDraw(TH1F *hist, TCanvas *canvas, double first, double last, double mean, 
                          double xmin, double bin, float percent)
{
  int nPointX = 1000; 
  int nPointY = 1000; 
  double x = 0;
  Float_t y = 0;

  TLegend *leg = new TLegend(0.6,0.6,0.89,0.89);
  TGraph *gr = new TGraph(nPointX*nPointY);
   gr->SetMarkerColor(46);
   gr->SetLineColor(46);
   gr->SetFillColor(46);

  for(int i = 0; i <= nPointX; i++) {
      x = first + (last-first)/nPointX*i;
      y = hist->GetBinContent(((x-xmin)/bin)+1);
      for(int j = 0; j <= nPointY; j++) {
        gr->SetPoint(i*nPointX + j, x, y/nPointY*j);
      }
   }
   leg->AddEntry(hist, "hist" ,"f");
   leg->AddEntry((TObject*)0, Form("mean=%.2f, GeV", hist->GetMean()), "");
   leg->AddEntry((TObject*)0, Form("RMS=%.2f, GeV", hist->GetRMS()), "");
 
   leg->AddEntry(gr, Form("integral_percent=%.2f", percent) ,"f");
   leg->AddEntry((TObject*)0, Form("meanCenter=%.2f, GeV", hist->GetBinCenter(hist->GetMaximumBin())), "");
   leg->AddEntry((TObject*)0, Form("range=%.2f, GeV", last - first), "");
   
   hist->Draw();
   gr->Draw("P");
   leg->Draw();

}

/*! \fn integra    
* \brief Function for find min diapason for percentage of the number of events.
* 
* \param TH1F *hist - data hist 
* \param double binWidth - bin size 
* \param double percent - percentage of the number of events for integral   

* \param double *first - (return) first value of interval 
* \param double *last - (return) last value of interval 
* \param double *centerMaxBin - (return) centerMaxBin 
*/
void JetAlgo::integra(TH1F *hist, float binWidth, float percent, 
                      double *first, double *last, double *centerMaxBin) 
{
  double center = hist->GetBinCenter(hist->GetMaximumBin());
  *centerMaxBin = center;
  
  int sumEntries = 1;
  int entries = (int)(hist->GetEntries()/100*percent);
  int indexBinLeft = hist->GetMaximumBin();
  int indexBinRight = hist->GetMaximumBin();
  int stepsLeft = 0;
  int stepsRight = 0;
  
  double carriageRight = hist->GetBinCenter(indexBinLeft);
  double carriageLeft = hist->GetBinCenter(indexBinRight);
  double centerBinLeft = hist->GetBinCenter(indexBinLeft);
  double centerBinRight = hist->GetBinCenter(indexBinRight);
  
  int numberEventsBinLeft = hist->GetBinContent(indexBinLeft);
  int numberEventsBinRight = hist->GetBinContent(indexBinRight);
  
  double stepWidthLeft = binWidth/numberEventsBinLeft;
  double stepWidthRight = binWidth/numberEventsBinRight;
  
  while(sumEntries <= entries) {
    if((carriageLeft <= (centerBinLeft - binWidth/2))) {
      indexBinLeft = indexBinLeft - 1;
      centerBinLeft = hist->GetBinCenter(indexBinLeft);
      numberEventsBinLeft = hist->GetBinContent(indexBinLeft);
      stepWidthLeft = binWidth/numberEventsBinLeft;
    }
    if((carriageRight >= (centerBinRight + binWidth/2))) {
      indexBinRight = indexBinRight + 1;
      centerBinRight = hist->GetBinCenter(indexBinRight);
      numberEventsBinRight = hist->GetBinContent(indexBinRight);
      stepWidthRight = binWidth/numberEventsBinRight;
    }

    if(stepWidthLeft == stepWidthRight) {
      if(stepsLeft <= stepsRight) {
        stepsRight = stepsRight + 1;
        carriageRight = carriageRight + stepWidthRight; 
        sumEntries = sumEntries + 1;
      }
      if(stepsLeft > stepsRight) {
        stepsLeft = stepsLeft + 1;
        carriageLeft = carriageLeft - stepWidthLeft; 
        sumEntries = sumEntries + 1;
      }
    }
    if(stepWidthLeft < stepWidthRight) {
        stepsLeft = stepsLeft + 1;
        carriageLeft = carriageLeft - stepWidthLeft; 
        sumEntries = sumEntries + 1;
    }
    if(stepWidthLeft > stepWidthRight) {
        stepsRight = stepsRight + 1;
        carriageRight = carriageRight + stepWidthRight; 
        sumEntries = sumEntries + 1;
    }
  }
  *first = carriageLeft;
  *last = carriageRight;
}

/*! \fn integra    
* \brief Function for find min diapason for percentage of the number of events.
* 
* \param vector<double> *data - data for calculate 
* \param double percent - percentage of the number of events for integral   

* \param double *first - (return) first value of interval 
* \param double *last - (return) last value of interval 
*/
void JetAlgo::integra(vector<double> *data, float percent, 
                      double *first, double *last) 
{
  
  double deltaEntries = numeric_limits<double>::max();
  double deltaEntriesT = 0;
  
  int window = (int)(data->size()/100*percent);
  sort(data->begin(), data->end());
  double *array = data->data();
 
  int Min;
  int Max;
  for(int i = 0; i < data->size() - window/2; i++) {
    if(i - window/2 >= 0) {
      Min = i - window/2; 
      Max = i + window/2; 
    }
    if(i - window/2 < 0) {
      Min = 0; 
      Max = i+window/2 + abs(i-window/2); 
    }
      deltaEntriesT =  abs(array[Max] - array[Min]);
    if(deltaEntries >= deltaEntriesT) {
      deltaEntries = deltaEntriesT;
      *first = array[Min];
      *last = array[Max];
    }
  }
}

/*! \fn meanWeight    
* \brief Function for find mean between three max bin.
* 
* \param TH1F *hist - data hist 

* \param double *mean - (return) mean value
*/
void JetAlgo::meanWeight(TH1F *hist, double *mean)
{
  double weight = hist->GetBinContent(hist->GetMaximumBin());
  double weightLow = hist->GetBinContent(hist->GetMaximumBin()-1);
  double weightHigh = hist->GetBinContent(hist->GetMaximumBin()+1);

  double center = hist->GetBinCenter(hist->GetMaximumBin());
  double centerLow = hist->GetBinCenter(hist->GetMaximumBin()-1);
  double centerHigh = hist->GetBinCenter(hist->GetMaximumBin()+1);
  
  *mean = (weight*center + weightLow*centerLow + weightHigh*center)/(weight+weightLow+weightHigh);
}

/*! \fn analysis    

 * \brief Function analysis jet's algorithm from loop function
* 
* \param string fileName - name for create file name 
* \param unsigned short numberSreps - number of steps (from loop fun) 
* \param float stepR - value fro radius step 
* \param float stepPT - value fro PT step 
* \param float stepEta - value fro Eta step 
* \param float percent - percent fro integral 

*/
int JetAlgo::analysis(string fileName, unsigned short numberSteps,
                   float stepR, double stepPT, double stepEta, float percent) 
{
  vector<fastjet::PseudoJet> pfcands;
 
  TFile *file = new TFile(fileName.c_str(), "UPDATE");
  if(!file) {
    cout << "Can not open file" << endl;
    return 0;
  }
  
  TTree *inputTreeRec = (TTree*)file->Get("Rec");
  if(!inputTreeRec) {
    cout << "Can not open inputput tree" << endl;
    return 0;
  }
  TTree *inputTreeGen = (TTree*)file->Get("Gen");
  if(!inputTreeGen) {
    cout << "Can not open inputput tree" << endl;
    return 0;
  }

  unsigned long nentries = 0; 
  if(inputTreeGen->GetEntries() == inputTreeRec->GetEntries())
    nentries = inputTreeGen->GetEntries();
  else {
    cout << "Not correct data for analysis" << endl;
    return 0;
  }
  
  int M = 0;
  if(rFlag == true) M = M + 1;
  if(jjFlag == true) M = M + 1;
  if(wjFlag == true) M = M + 1;
  
  string name;
  float radius[M]; float PT[M]; float Eta[M];

  TLeaf *leafGen[M];
  TLeaf *leafRec[M];

  vector<vector<double>> dataGen;
  dataGen.resize(M);
  vector<vector<double>> dataRec;
  dataRec.resize(M);

  TH1F *histGen[numberSteps][M];
  float xMinGen[numberSteps][M]; float xMaxGen[numberSteps][M];
  TH1F *histRec[numberSteps][M];
  float xMinRec[numberSteps][M]; float xMaxRec[numberSteps][M];
 
  TH1F *histDiffRG[numberSteps][M];
  
  TCanvas *canvasResultGen[numberSteps][M];
  TCanvas *canvasResultRec[numberSteps][M];

  TGraph *resolutionGen[M];
  TGraph *maxBinGen[M];
  TGraph *resolutionRec[M];
  TGraph *maxBinRec[M];

  int binGen[numberSteps][M]; int binRec[numberSteps][M]; double binWidth[numberSteps][M];
  double firstGen[numberSteps][M]; double lastGen[numberSteps][M]; double maxGen[numberSteps][M];
  double firstRec[numberSteps][M]; double lastRec[numberSteps][M]; double maxRec[numberSteps][M];
  double lastDiffRG[numberSteps][M]; double firstDiffRG[numberSteps][M];

  for(unsigned short m = 0; m < M; m++) {
    if(m == 0) { radius[m] = rRadius; PT[m] = rPT; Eta[m] = rEta; }
    if(m == 1) { radius[m] = jjLastRadius; PT[m] = jjPT; Eta[m] = jjEta; }
    if(m == 2) { radius[m] = wjLastRadius; PT[m] = wjPT; Eta[m] = wjEta; }
    for(unsigned short i = 0; i < numberSteps; i++) {
    	if(m == 0)
        name = Form("%s_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), radius[m], PT[m], Eta[m]); 
      if(m == 1) 
        name = Form("jetOfJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), jjRadius, radius[m], PT[m], Eta[m]);
      if(m == 2) 
        name = Form("wideJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), wjRadius, radius[m], PT[m], Eta[m]);
      
      xMinGen[i][m] = numeric_limits<float>::max(); xMaxGen[i][m] = numeric_limits<float>::min();
      xMinRec[i][m] = numeric_limits<float>::max(); xMaxRec[i][m] = numeric_limits<float>::min();
      
      leafGen[m] = inputTreeGen->GetLeaf(name.c_str());
      leafRec[m] = inputTreeRec->GetLeaf(name.c_str());
    
      histDiffRG[i][m] = new TH1F(Form("diffRecGen_%s", name.c_str()), 
                                  Form("diffRecGen_%s; invMassRec-invMassGen, GeV; Nentries", name.c_str()),                               
                                  100, 1, -1);
     
      for(unsigned long j = 0; j < nentries; j++) {
        
        leafGen[m]->GetBranch()->GetEntry(j);
        
        dataGen[m].push_back(leafGen[m]->GetValue());

        if(xMinGen[i][m] > leafGen[m]->GetValue())
          xMinGen[i][m] = leafGen[m]->GetValue();
        if(xMaxGen[i][m] < leafGen[m]->GetValue())
          xMaxGen[i][m] = leafGen[m]->GetValue();
       
        leafRec[m]->GetBranch()->GetEntry(j);
        
        dataRec[m].push_back(leafRec[m]->GetValue());
        
        if(xMinRec[i][m] > leafRec[m]->GetValue())
          xMinRec[i][m] = leafRec[m]->GetValue();
        if(xMaxRec[i][m] < leafRec[m]->GetValue())
          xMaxRec[i][m] = leafRec[m]->GetValue();
        
        histDiffRG[i][m]->Fill(leafRec[m]->GetValue()-leafGen[m]->GetValue());
      }
      
        integra(&dataGen[m], percent, &firstGen[i][m], &lastGen[i][m]);
        dataGen[m].clear();

        integra(&dataRec[m], percent, &firstRec[i][m], &lastRec[i][m]);
        dataRec[m].clear();
        
        binWidth[i][m] = (lastRec[i][m]-firstRec[i][m])-(lastGen[i][m]-firstGen[i][m]);
        
      radius[m] = radius[m] + stepR; PT[m] = PT[m] + stepPT; Eta[m] = Eta[m] + stepEta;
    }
  }
  
  for(unsigned short m = 0; m < M; m++) {
    if(m == 0) { radius[m] = rRadius; PT[m] = rPT; Eta[m] = rEta; }
    if(m == 1) { radius[m] = jjLastRadius; PT[m] = jjPT; Eta[m] = jjEta; }
    if(m == 2) { radius[m] = wjLastRadius; PT[m] = wjPT; Eta[m] = wjEta; }
    for(unsigned short i = 0; i < numberSteps; i++) {
      if(m == 0)
        name = Form("%s_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), radius[m], PT[m], Eta[m]); 
      if(m == 1) 
        name = Form("jetOfJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), jjRadius, radius[m], PT[m], Eta[m]);
      if(m == 2) 
        name = Form("wideJet_%s_%.2fR_%.2fR_%.2fPT_%.2fEta", algorithm.c_str(), wjRadius, radius[m], PT[m], Eta[m]);
      
      leafGen[m] = inputTreeGen->GetLeaf(name.c_str());
      leafRec[m] = inputTreeRec->GetLeaf(name.c_str());
        
      binGen[i][m] =  (int)((xMaxGen[i][m]-xMinGen[i][m])/binWidth[i][m]);

      histGen[i][m] = new TH1F(Form("gen_%s", name.c_str()), 
                               Form("gen_%s; invMass, GeV; hist(invMass)/(bin*Nentries), 1/GeV", name.c_str()),                               
                               binGen[i][m], xMinGen[i][m], xMaxGen[i][m]);
        
      binRec[i][m] = (int)((xMaxRec[i][m]-xMinRec[i][m])/binWidth[i][m]);

      histRec[i][m] = new TH1F(Form("rec_%s", name.c_str()), 
                               Form("rec_%s; invMass, GeV; hist(invMass)/(bin*Nentries), 1/GeV", name.c_str()),                               
                               binRec[i][m], xMinRec[i][m], xMaxRec[i][m]);
      
      for(unsigned long j = 0; j < nentries; j++) {
          leafGen[m]->GetBranch()->GetEntry(j);
          histGen[i][m]->Fill(leafGen[m]->GetValue());
          leafRec[m]->GetBranch()->GetEntry(j);
          histRec[i][m]->Fill(leafRec[m]->GetValue());
      }
      
      double scaleGen = ((xMaxGen[i][m]-xMinGen[i][m])/binGen[i][m])*nentries;
      histGen[i][m]->Scale(1/scaleGen);
      meanWeight(histGen[i][m], &maxGen[i][m]);
      double scaleRec = ((xMaxRec[i][m]-xMinRec[i][m])/binRec[i][m])*nentries;
      histRec[i][m]->Scale(1/scaleRec);
      meanWeight(histRec[i][m], &maxRec[i][m]);

      canvasResultGen[i][m] = new TCanvas(Form("integral_gen_%s", name.c_str()),
                                          Form("integral_gen_%s", name.c_str()),
                                          200,10,700,500);
      integraDraw(histGen[i][m], canvasResultGen[i][m], firstGen[i][m], lastGen[i][m], maxGen[i][m], xMinGen[i][m], (xMaxGen[i][m]-xMinGen[i][m])/binGen[i][m], percent);

      canvasResultRec[i][m] = new TCanvas(Form("integral_rec_%s", name.c_str()),
                                          Form("integral_rec_%s", name.c_str()),
                                          200,10,700,500);
      integraDraw(histRec[i][m], canvasResultRec[i][m], firstRec[i][m], lastRec[i][m], maxRec[i][m], xMinRec[i][m], (xMaxRec[i][m]-xMinRec[i][m])/binRec[i][m], percent);
      
      radius[m] = radius[m] + stepR; PT[m] = PT[m] + stepPT; Eta[m] = Eta[m] + stepEta;
    }
  }
  
  file->mkdir("hist_differenceRecGen"); 
  file->mkdir("hist_Gen"); 
  file->mkdir("hist_Rec"); 
  file->mkdir("integral_Gen"); 
  file->mkdir("integral_Rec"); 
  file->mkdir("resolution_meanCenter_Gen"); 
  file->mkdir("resolution_meanCenter_Rec"); 

  for(unsigned short m = 0; m < M; m++) {
    if(m == 0) { name = ""; radius[m] = rRadius; PT[m] = rPT; Eta[m] = rEta; }
    if(m == 1) { name = "jetOfJet"; radius[m] = jjLastRadius; PT[m] = jjPT; Eta[m] = jjEta; }
    if(m == 2) { name = "wideJet"; radius[m] = wjLastRadius; PT[m] = wjPT; Eta[m] = wjEta; }
    
    resolutionGen[m] = new TGraph(numberSteps);
    resolutionGen[m]->SetMarkerStyle(8);
    resolutionGen[m]->SetName(Form("gen_%s_%s_resolution; radius; resolution, GeV", name.c_str(), algorithm.c_str()));
    resolutionGen[m]->SetTitle(Form("gen_%s_%s_resolution; radius; resolution, GeV;", name.c_str(), algorithm.c_str()));
    
    maxBinGen[m] = new TGraph(numberSteps);
    maxBinGen[m]->SetMarkerStyle(8);
    maxBinGen[m]->SetName(Form("gen_%s_%s_meanCenter; radius; maxBinCenter, GeV", name.c_str(), algorithm.c_str()));
    maxBinGen[m]->SetTitle(Form("gen_%s_%s_meanCenter; radius; maxBinCenter, GeV", name.c_str(), algorithm.c_str()));
    
    resolutionRec[m] = new TGraph(numberSteps);
    resolutionRec[m]->SetMarkerStyle(8);
    resolutionRec[m]->SetName(Form("rec_%s_%s_resolution; radius; resolution, GeV", name.c_str(), algorithm.c_str()));
    resolutionRec[m]->SetTitle(Form("rec_%s_%s_resolution; radius; resolution, GeV", name.c_str(), algorithm.c_str()));
    
    maxBinRec[m] = new TGraph(numberSteps);
    maxBinRec[m]->SetMarkerStyle(8);
    maxBinRec[m]->SetName(Form("rec_%s_%s_meanCenter; radius; maxBinCenter, GeV", name.c_str(), algorithm.c_str()));
    maxBinRec[m]->SetTitle(Form("rec_%s_%s_meanCenter; radius; maxBinCenter, GeV", name.c_str(), algorithm.c_str()));

    float xMinRangeMeanGen = numeric_limits<float>::max(); float xMaxRangeMeanGen = numeric_limits<float>::min();
    float xMinRangeMeanRec = numeric_limits<float>::max(); float xMaxRangeMeanRec = numeric_limits<float>::min();
    for(unsigned short i = 0; i < numberSteps; i++) {

      resolutionGen[m]->SetPoint(i, radius[m], lastGen[i][m]-firstGen[i][m]);
      maxBinGen[m]->SetPoint(i, radius[m], maxGen[i][m]);
        if(xMinRangeMeanGen > maxGen[i][m])
          xMinRangeMeanGen = maxGen[i][m];
        if(xMaxRangeMeanGen < maxGen[i][m])
          xMaxRangeMeanGen = maxGen[i][m];
      resolutionRec[m]->SetPoint(i, radius[m], lastRec[i][m]-firstRec[i][m]);
      maxBinRec[m]->SetPoint(i, radius[m], maxRec[i][m]);
        if(xMinRangeMeanRec > maxRec[i][m])
          xMinRangeMeanRec = maxRec[i][m];
        if(xMaxRangeMeanRec < maxRec[i][m])
          xMaxRangeMeanRec = maxRec[i][m];

       file->cd("hist_Gen");
       histGen[i][m]->Write();
       
       file->cd("hist_Rec");
       histRec[i][m]->Write();
       
       file->cd("hist_differenceRecGen"); 
       histDiffRG[i][m]->Write();
       
       file->cd("integral_Gen"); 
       histGen[i][m]->SetStats(0); 
       canvasResultGen[i][m]->Write();
       
       file->cd("integral_Rec"); 
       histRec[i][m]->SetStats(0); 
       canvasResultRec[i][m]->Write();
      
       radius[m] = radius[m] + stepR; PT[m] = PT[m] + stepPT; Eta[m] = Eta[m] + stepEta;
    }
   
    maxBinGen[m]->SetMaximum(xMaxRangeMeanGen + 1000);
    maxBinGen[m]->SetMinimum(xMinRangeMeanGen - 1000);
    maxBinRec[m]->SetMaximum(xMaxRangeMeanRec + 1000);
    maxBinRec[m]->SetMinimum(xMinRangeMeanRec - 1000);
    
    file->cd("resolution_meanCenter_Gen"); 
    resolutionGen[m]->Write();
    maxBinGen[m]->Write();
    file->cd("resolution_meanCenter_Rec"); 
    resolutionRec[m]->Write();
    maxBinRec[m]->Write();
  }

  for(unsigned short m = 0; m < M; m++) {
    leafGen[m]->Delete();
    leafRec[m]->Delete();
    resolutionGen[m]->Delete();
    maxBinGen[m]->Delete();
    resolutionRec[m]->Delete();
    maxBinRec[m]->Delete();
    for(unsigned short i = 0; i < numberSteps; i++) {
      histGen[i][m]->Delete();
      histRec[i][m]->Delete();
      histDiffRG[i][m]->Delete();
      canvasResultGen[i][m]->Close();
      canvasResultRec[i][m]->Close();
    
    }
  }
  inputTreeRec->Delete();
  inputTreeGen->Delete();
  file->Close();
  file->Delete();
 
  return 1;
}
