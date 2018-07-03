#include <iostream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TF1.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include <TROOT.h>
#include "TSystem.h"
#include "TGaxis.h"
#include <cstdlib>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TH2D.h" 

void Enable_and_Set_Branches(TTree* & tree);

// Setting parameters //////////////////////////////////////////////////

// range to the number of electrons per cluster
  int emin = 1300; int emax = 17000; 
  
  int ohdu_numer = 4;
  
  std::string MCInputFile = "Output.root";
 
////////////////////////////////////////////////////////////////////////

  int Entries_mc = 1;
  int Entries_exp = 1;
  int xmin = 0; // range for histograms
  int xmax =20; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;
  
  int nbins = xmax+1;
   
////////////////////////////////////////////////////////////////////////

  int runID; int ohdu; int expoStart; 
  int nSat; int flag; 
  int xMin; int xMax;
  int yMin; int yMax;
  Float_t e; Float_t n;
  Float_t xBary; Float_t yBary;
  Float_t xVar; Float_t yVar;
  int  nSavedPix;
  int xPix; int yPix;
  Float_t ePix;
    
using namespace std;

void MakePlots(){  
     
// Get input files//////////////////////////////////////////////////////

// Experimental Data ///////////////////////////////////////////////////
    TFile * f_exp = TFile::Open("55Fe_exp.root");
    if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
    TTree * texp = (TTree*) f_exp->Get("hitSumm");
    TH1D * h_exp_n  =  new TH1D("h_exp_n"    , "Distribucion de tamanio de clusters", nbins, xmin, xmax);
    h_exp_n -> Sumw2();

	int Entries_exp = texp -> GetEntries();
	cout<<"Entries in experimental data file: "<<Entries_exp<<endl;
	
	Enable_and_Set_Branches(texp); 
	
// Monte Carlo Data ////////////////////////////////////////////////////
    
    TFile * f_mc = TFile::Open("Out_N0=200_DC=500_A=180_B=13.root");
    
    if (!f_mc->IsOpen()) {std::cerr << "ERROR: cannot open the root file with MC data" << std::endl;}   
    TTree * tmc    = (TTree*) f_mc->Get("hitSumm");
    TH1D * h_mc_n      =  new TH1D("h_mc_n"        , "Distribucion de tamanio de clusters", nbins, xmin, xmax);
    h_mc_n -> Sumw2();
    
    int Entries_mc = tmc -> GetEntries();
    cout<<"Enties in MC file: "<<Entries_mc<<endl;

	Enable_and_Set_Branches(tmc);
	
// Get information from trees///////////////////////////////////////////

// Experimental Data ///////////////////////////////////////////////////	
    for(int i_event=0;i_event<Entries_exp; i_event++){
   
    texp->GetEntry(i_event);
	if (ohdu == ohdu_numer) {
		if (e>emin & e<emax){  // Condition over the number of electrons
			if (xBary>10 & xBary<490){  // Condition over position
				if (yBary>5 & yBary<45){  // Condition over position
					// Fill the histogram with the variable n
					h_exp_n -> Fill(n); 
				}
			}
		}	
	}
    }
        
// MonteCarlo Data /////////////////////////////////////////////////////           
    for(int i_event=0;i_event<Entries_mc; i_event++){
    
    tmc->GetEntry(i_event);

		if (e>emin & e<emax){  // Condition over the number of electrons
			if (xBary>10 & xBary<490){  // Condition over position
				if (yBary>5 & yBary<45){  // Condition over position
					// Fill the histogram with the variable n
					h_mc_n -> Fill(n); 
				}
			}
		}	
    }
            
// Plot histograms /////////////////////////////////////////////////////
    TCanvas*  cc_s   = new TCanvas("cc_s","cc_s",1000,500);
    cc_s->cd();
       
    h_exp_n -> SetLineColor(kRed +3);         h_mc_n -> SetLineColor(kBlue);
	h_exp_n -> SetMarkerColor(kRed +3);       h_mc_n -> SetMarkerColor(kBlue);
	h_exp_n -> SetMarkerSize(0.8);              h_mc_n -> SetMarkerSize(0.8);
	h_exp_n -> SetMarkerStyle(22);              h_mc_n -> SetMarkerStyle(23);
    
    int norm=1; //Normalization    
    cc_s->SetLogy(1);
    Double_t scale_exp = norm/h_exp_n->Integral();
	h_exp_n->Scale(scale_exp);
    h_exp_n ->Draw("HIST E1");
    //cc_s->SaveAs( "/home/dario/CCD/toy/figures/Dist_n_exp.png"); 
    
    Double_t scale_mc = norm/h_mc_n->Integral();
	h_mc_n->Scale(scale_mc);
    h_mc_n ->Draw("HIST E1 same");
    //cc_s->SaveAs( "/home/dario/CCD/toy/figures/Dist_n_MC.png"); 

    //delete cc_s;
    
}

void Enable_and_Set_Branches(TTree* & tree){

  tree->SetBranchStatus("*",1); //enable all branches

  tree->SetBranchAddress ("runID",&runID);
  tree->SetBranchAddress ("ohdu",&ohdu);
  tree->SetBranchAddress ("expoStart",&expoStart);
  tree->SetBranchAddress ("nSat",&nSat);
  tree->SetBranchAddress ("flag",&flag);
  tree->SetBranchAddress ("xMin",&xMin);
  tree->SetBranchAddress ("xMax",&xMax);
  tree->SetBranchAddress ("yMin",&yMin);
  tree->SetBranchAddress ("yMax",&yMax);
  tree->SetBranchAddress ("e",&e);
  tree->SetBranchAddress ("n",&n);
  tree->SetBranchAddress ("xBary",&xBary);
  tree->SetBranchAddress ("yBary",&yBary);
  tree->SetBranchAddress ("xVar",&xVar);
  tree->SetBranchAddress ("yVar",&yVar);
  tree->SetBranchAddress ("nSavedPix",&nSavedPix);
  tree->SetBranchAddress ("xPix",&xPix);
  tree->SetBranchAddress ("yPix",&yPix);
  tree->SetBranchAddress ("ePix",&ePix);
}
