//-------------------------------------------------------
// Description: main file for analysis algorithms
// Authors:  Shchablo, Shchablo@gmail.com
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

// Analyses class
#include "JetAlgo.hh"

using namespace std;

// Main function that runs the analysis algorithm on the
// specified input files
int main(int argc, char* argv[]) 
{
  cout <<"#-------------------------------------------------------------------------" << endl;
  cout <<"#-------------------------------------------------------------------------" << endl;
  cout <<"#                     JetAlgo, Analysis jet's algorithms"                   << endl;
  cout <<"#-------------------------------------------------------------------------" << endl;

  //
  /* Algorithm */
  string algorithm = "antiKt"; // "antiKt"; "kt"; "cam".
  int numberSteps = 13; // numbers of circle
  float stepR = 0.1;
  double PT = 30; double stepPT = 0;
  double Eta = 3; double stepEta = 0;

  /* just Algorithm (antiKt)  */
  float percent = 25;
  float rRadius = 0.4;
  double rPT; double rEta;
  /* JetOfJet */
  float jjRadius = 0.4; float jjLastRadius = 0.4;
  double jjPT; double jjEta;
  /* WideJet */
  float wjRadius = 0.4; float wjLastRadius = 0.4;
  float wjPT; float wjEta;

  rPT = PT; jjPT = PT; wjPT = PT;
  rEta = Eta; jjEta = Eta; wjEta = Eta;
  //

  if ( argc < 2 ) {
    cout << "To run CMSApp please specify the input file" << endl; 
    cout << "Example:  ./jetAlgo inputFile.root" << endl;
    cout << "OPTIONS:        " << endl;
    cout << "-output             Name of the output file" << endl;
    cout << "--tree=[name]       Choose name of tree for read" << endl;
    
    cout << "--verbose           Increase verbosity level" << endl;

    cout << "--algorithm=[algo]  Choose basic algorithm. algo: antiKt, kt, cam" << endl;
    cout << "--gen               Run genetarion level Analysis" << endl;
    cout << "--rec               Run reconstruction level Analysis" << endl;
    cout << "--analysis          Run Analysis data" << endl;
   
    cout << "--allJet            Run all jet's algorithm" << endl;
    cout << "--r                 Turn on just jet algorithm" << endl;
    cout << "--jj                Turn on jet of jet algorithm" << endl;
    cout << "--wj                Turn on wide jet algorithm" << endl;

    cout << "--ns=[value]        Number of steps" << endl;
    cout << "--stepR=[value]     Set value for radius step" << endl;
    cout << "--stepPT=[value]    Set value for pt step" << endl;
    cout << "--stepEta=[value]   Set value for eta step" << endl;
    cout << "--percent[value]    Set percent number of events for integral" << endl;
    return 1;
  }

  char inputBsmFile[150];
  string inputTreeName = "Delphes";
  strcpy(inputBsmFile, argv[1]);

  char outFileName[150] = "none";

  bool verbose  = false;
  bool writeOut = false;
  bool rJet = false;
  bool jjJet = false;
  bool wjJet = false;
  bool gen = false;
  bool rec = false;
  bool analysis = false;

  for(int i = 1; i < argc; i++) {
    if(strncmp(argv[i], "-output", 7) == 0) {
      sscanf(argv[i], "-output=%s", outFileName);
      writeOut = true;
    }
    if(strncmp(argv[i],"--tree=", 7) == 0) {
      char tmp[128];
      sscanf(argv[i], "--tree=%s", tmp);
      inputTreeName = tmp;
    }
    if(strncmp(argv[i], "--verbose", 9) == 0)
      verbose = true;
    if(strncmp(argv[i], "--allJet", 8) == 0) {
       rJet = true;
       jjJet = true;
       wjJet = true;
    }
    if(strncmp(argv[i], "--r", 3) == 0)
      rJet = true;
    if(strncmp(argv[i], "--jj", 4) == 0)
      jjJet = true;
    if(strncmp(argv[i], "--wj", 4) == 0)
      wjJet = true;
    if(strncmp(argv[i], "--gen", 5) == 0)
      gen = true;
    if(strncmp(argv[i], "--rec", 5) == 0)
      rec = true;
    if(strncmp(argv[i], "--analysis", 10) == 0)
      analysis = true;
    if(strncmp(argv[i], "--ns=", 5) == 0)
      sscanf(argv[i],"--ns=%d", &numberSteps);
    if(strncmp(argv[i], "--algorithm=", 12) == 0) {
      char tmp[128];
      sscanf(argv[i], "--algorithm=%s", tmp);
      algorithm = tmp;
    }
    if(strncmp(argv[i], "--stepR=", 8) == 0)
      sscanf(argv[i], "--stepR=%f", &stepR);
    if(strncmp(argv[i], "--stepPT=", 9) == 0)
      sscanf(argv[i], "--stepPT=%lf", &stepPT);
    if(strncmp(argv[i], "--stepEta=", 10) == 0)
      sscanf(argv[i], "--stepEta=%lf", &stepEta);
    if(strncmp(argv[i],"--percent=", 10) == 0)
      sscanf(argv[i], "--percent=%f", &percent);
  }
   
  if(strncmp(inputBsmFile, "none", 4) != 0) {
    
    if((rJet == true || jjJet == true || wjJet == true) 
        && (rec == true || gen == true)) {

      /* create class */
      JetAlgo jetAlgo;
      jetAlgo.setVerbose(verbose);

      jetAlgo.setAlgorithm(algorithm);
      if(rJet == true)
        jetAlgo.setRadius(rRadius, rPT, rEta);
      if(jjJet == true)
        jetAlgo.setJetOfJet(jjRadius, jjLastRadius, jjPT, jjEta);
      if(wjJet == true)
        jetAlgo.setWideJet(wjRadius, wjLastRadius, wjPT, wjEta);

      if(rec == true) {
        /* names of leaf for read from TTree */
        vector<vector<string>> leafs;
        leafs.resize(3);
        for(int i = 0; i < leafs.size(); i++) {
            leafs[i].resize(3);
        }
        leafs[0][0] = "ChargedHadron.PT";
        leafs[0][1] = "ChargedHadron.Eta";
        leafs[0][2] = "ChargedHadron.Phi";

        leafs[1][0] = "NeutralHadron.ET";
        leafs[1][1] = "NeutralHadron.Eta";
        leafs[1][2] = "NeutralHadron.Phi";

        leafs[2][0] = "Photon.PT";
        leafs[2][1] = "Photon.Eta";
        leafs[2][2] = "Photon.Phi";

          cout <<"#-------------------------------------------------------------------------" << endl;
          cout << "# Reconstruction level" << endl;
          cout << "# InputFile: " << inputBsmFile << endl;
          cout << "# InputTree: " << "Delphes" << endl;
          for(unsigned short i = 0; i < leafs.size(); i++) {
            cout << "# leafs:";
            for(unsigned short j = 0; j < leafs[i].size(); j++)
              cout << " " << leafs[i][j];
            cout << endl;
          }
          cout << "# OutputFile: " << outFileName << endl;
          cout << "# OutputTree: " << "Rec" << endl;
          cout << "# Parameters: " 
               << "numberSteps=" << numberSteps
               << " stepR=" << stepR
               << " stepPT=" << stepPT
               << " stepEta=" << stepEta << endl;
          if(rJet == true) {
            cout << "# justAlgorithm, firstParam: " 
                 << "radius=" << rRadius 
                 << " PT=" << rPT 
                 << " Eta=" << rEta << endl;
          }
          if(jjJet == true) {
            cout << "# jetOfJet, firstParam: " 
                 << "firstRadius=" << jjRadius 
                 << " lastRadius=" << jjLastRadius 
                 << " PT=" << jjPT 
                 << " Eta=" << jjEta << endl;
          }
          if(wjJet == true) {
            cout << "# wideJet, firstParam: " 
                 << "firstRadius=" << wjRadius 
                 << " lastRadius=" << wjLastRadius 
                 << " PT=" << wjPT 
                 << " Eta=" << wjEta << endl;
          }
          cout <<"#-------------------------------------------------------------------------" << endl;
        int flag = jetAlgo.setInputTree(inputBsmFile, inputTreeName.c_str(), leafs, "Rec");
        if(flag) {
          int algoFlag = jetAlgo.loop(outFileName, numberSteps, stepR, stepPT, stepEta);
          if(!algoFlag)
            return 1;
        }
        else
          return 1;
      }
      if(gen == true) {
        /* names of leaf for read from TTree */
        vector<vector<string>> leafs;
        leafs.resize(1);
        leafs[0].resize(6);

        leafs[0][0] = "Particle.Px";
        leafs[0][1] = "Particle.Py";
        leafs[0][2] = "Particle.Pz";
        leafs[0][3] = "Particle.E";
        leafs[0][4] = "Particle.Status";
        leafs[0][5] = "Particle.PID";
        
          cout <<"#-------------------------------------------------------------------------" << endl;
          cout << "# Generation level" << endl;
          cout << "# InputFile: " << inputBsmFile << endl;
          cout << "# InputTree: " << "Delphes" << endl;
          for(unsigned short i = 0; i < leafs.size(); i++) {
            cout << "# leafs:";
            for(unsigned short j = 0; j < leafs[i].size(); j++)
              cout << " " << leafs[i][j];
            cout << endl;
          }
          cout << "# OutputFile: " << outFileName << endl;
          cout << "# OutputTree: " << "Gen" << endl;
          cout << "# Parameters: " 
               << "numbersteps=" << numberSteps
               << " stepR=" << stepR
               << " stepPT=" << stepPT
               << " stepEta=" << stepEta << endl;
          if(rJet == true) {
            cout << "# justAlgorithm, firstParam: " 
                 << "radius=" << rRadius 
                 << " PT=" << rPT 
                 << " Eta=" << rEta << endl;
          }
          if(jjJet == true) {
            cout << "# jetOfJet, firstParam: " 
                 << "firstRadius=" << jjRadius 
                 << " lastRadius=" << jjLastRadius 
                 << " PT=" << jjPT 
                 << " Eta=" << jjEta << endl;
          }
          if(wjJet == true) {
            cout << "# wideJet, firstParam: " 
                 << "firstRadius=" << wjRadius 
                 << " lastRadius=" << wjLastRadius 
                 << " PT=" << wjPT 
                 << " Eta=" << wjEta << endl;
          }
          cout <<"#-------------------------------------------------------------------------" << endl;
      
        int flag = jetAlgo.setInputTree(inputBsmFile, inputTreeName.c_str(), leafs, "Gen");
        if(flag) {
          int algoFlag = jetAlgo.loop(outFileName, numberSteps, stepR, stepPT, stepEta);
          if(!algoFlag)
            return 1;
        }
        else
          return 1;
      }
        if(analysis == true && rec == true && gen == true) {
          int algoFlag = jetAlgo.analysis(outFileName, numberSteps, stepR, stepPT, stepEta, percent);
          if(!algoFlag)
            return 1;
        }
        else
          cout << "Can't run analysis" << endl;
    }
  }
  cout << "The End" << endl;
}
