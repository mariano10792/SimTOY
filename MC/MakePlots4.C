#include <iostream>
#include <vector>
#include <string>
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
using namespace std;

void Enable_and_Set_Branches(TTree* & tree);

// Setting parameters //////////////////////////////////////////////////

// range for the number of electrons per cluster
  int emin =1500 ; int emax = 1680;
  int ohdu_numer = 2;
//number of bins to take into account for chi2
  int   bines = 100;
  float minePix = 0.5; // will be clasified as 1e-

  ////////////////////////////////////////////////////////////////////////


  int Entries_exp = 1;
  int xmin = 0; // range for histograms
  int xmax =10; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;

  int nbins = xmax*1+1;

  const int maxClusterSize = 50000;

////////////////////////////////////////////////////////////////////////

  int runID; int ohdu; int expoStart;
  int nSat; int flag;
  int xMin; int xMax;
  int yMin; int yMax;
  Float_t e; Float_t n;
  Float_t xBary; Float_t yBary;
  Float_t xVar; Float_t yVar;
  int  nSavedPix;
  int xPix[maxClusterSize];
  int yPix[maxClusterSize];
  Float_t ePix[maxClusterSize];

using namespace std;

////////////////////////////////////////////////////////////////////////
//  Older compilers
////////////////////////////////////////////////////////////////////////
string itos(const int i){
	ostringstream ossAux;
	ossAux << i;
	return ossAux.str();
}



void MakePlots3(){

// Experimental Data ///////////////////////////////////////////////////
// Get input files//////////////////////////////////////////////////////

cout<<"min ePix"<< minePix<<endl;

					TFile * f_exp = TFile::Open("/home/mariano/MEGAsync/images_from_mkids/21Oct2018_prueba_a2,5/T100/afterskipper2root/clustered/output_2_21Oct2018_prueba_a2,5_T100.root");
					if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
					TTree * texp = (TTree*) f_exp->Get("hitSumm");
					TH1D * h_exp_n  =  new TH1D("h_exp_n", "Distribucion de tamano de clusters", nbins, emin, emax);
					h_exp_n -> Sumw2();

					int Entries_exp = texp -> GetEntries();
					cout<<"Entries in experimental data file: "<<Entries_exp<<endl;

					Enable_and_Set_Branches(texp);

				// Get information from trees///////////////////////////////////////////

					for(int i_event=0;i_event<Entries_exp; i_event++){

					texp->GetEntry(i_event);

						if (ohdu == ohdu_numer) {
							if (e>emin && e<emax){  // number of electrons

								// Check if one of the pixels in the cluster is smaller that minePix
								bool noLowPixInCluster = true;
								for (int p = 0; p < nSavedPix; ++p){
									if(ePix[p]<minePix){
										noLowPixInCluster = false;
										break;
									}
								}
								bool noBadTransferInCluster = true;
								// for (int p = 0; p < nSavedPix; ++p){
								// 	if(xPix[p]==305){
								// 	noBadTransferInCluster = false;
								// 	break;
								// }
							}

								if (noLowPixInCluster){
									if (noBadTransferInCluster){
									if (xBary>10 && xBary<290){
										if (yBary>2 && yBary<28){
											if (xPix[0]>0){

											Fill the histogram with the variable ePix
											for (int p = 0; p < nSavedPix; ++p){
											h_exp_n -> Fill(ePix[p]);
												}
											}
										}
									}
								}
							}
						}
					}
				}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

				// Plot histograms /////////////////////////////////////////////////////

          TCanvas*  clusters   = new TCanvas("clusters","clusters",1200,1200);
					clusters->cd();
					gStyle->SetOptStat(0);
					//gPad->DrawFrame(0.5,0.85,3.5,1.1,"Summary of fit results;;Fit Bias");

					h_exp_n -> SetLineColor(kRed);
          h_exp_n -> SetLineWidth(2.0);
					h_exp_n -> SetMarkerColor(kRed);
					h_exp_n -> SetMarkerSize(2.0);
					h_exp_n -> SetMarkerStyle(24);

					//int norm=1; //Normalization
					//clusters->SetLogy(1);
					//Double_t scale_exp = norm/h_exp_n->Integral();
					//h_exp_n->Scale(scale_exp);
					h_exp_n->SetMaximum(3500);
					h_exp_n ->Draw("HIST E1");
					//clusters->SaveAs( "./figures/Dist_n_exp.png");

					vector<double> EXP(bines);
					vector<double> EXP_error(bines);

					for(int i_bin=1;i_bin<bines; i_bin++){
						EXP[i_bin]= (h_exp_n->GetBinContent(i_bin));
						EXP_error[i_bin] = h_exp_n->GetBinError(i_bin);
						cout<<EXP[i_bin]<<" "<<EXP_error[i_bin]<<endl;
					}


					TLegend *lg = new TLegend(0.35,0.7,0.9,0.9,NULL,"brNDC");
					lg->SetBorderSize(1);
					lg->SetTextFont(50);
					lg->SetTextSize(0.03146);
					lg->SetFillColor(10);
					lg->SetFillStyle(1001);
					//lg->AddEntry(h_mc_n,("MC N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)).c_str(),"l");
					//lg->AddEntry(h_mc_n,("emin="+itos(emin)+"_emax="+itos(emax)+"_minePix="+itos(minePix)+"_emin="+itos(emin)+"_emax="+itos(emax)).c_str(),"l");
          lg->AddEntry(h_exp_n,"Datos experimentales","l");

					lg->Draw("same");

					gPad->Update();
					//getchar();


				//}
			//}
		//}
	//}
}

////////////////////////////////////////////////////////////////////////

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
