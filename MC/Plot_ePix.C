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
  int emin =0 ; int emax = 200;
  int ohdu_numer = 4;
//number of bins to take into account for chi2
  int   bines = 7;
  float minePix = 10; // will be clasified as 1e-

  ////////////////////////////////////////////////////////////////////////



  int Entries_mc = 1;
  int Entries_exp = 1;
  int xmin = emin; // range for histograms
  int xmax =emax; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;

  int nbins = xmax*10+1;

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
  int J=0;
  

using namespace std;

////////////////////////////////////////////////////////////////////////
//  Older compilers
////////////////////////////////////////////////////////////////////////
string itos(const int i){
	ostringstream ossAux;
	ossAux << i;
	return ossAux.str();
}



void Plot_ePix(){

// Experimental Data ///////////////////////////////////////////////////
// Get input files//////////////////////////////////////////////////////

cout<<"min ePix"<< minePix<<endl;

					TFile * f_exp = TFile::Open("55Fe_exp.root");
					if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
					TTree * texp = (TTree*) f_exp->Get("hitSumm");
					TH1D * h_exp_n  =  new TH1D("h_exp_n", "Distribucion de tamanio de clusters", nbins, xmin, xmax);
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
								bool LowPixInCluster = false;
								for (int p = 0; p < nSavedPix; ++p){
									if(ePix[p]<minePix){
										LowPixInCluster = true;
										break;
									}
								}
								bool noBadTransferInCluster = true;
								for (int p = 0; p < nSavedPix; ++p){
									if(xPix[p]==305){
									noBadTransferInCluster = false;
									break;
								}
							}

								if (LowPixInCluster){
									if (noBadTransferInCluster){
									if (xBary>10 && xBary<490){
										if (yBary>5 && yBary<45){
											if (xPix[0]>250){
												if(n>0){
											//Fill the histogram with the variable ePix
											for (int p = 0; p < nSavedPix; ++p){
											h_exp_n -> Fill(ePix[p]);
											
											
											J++;
											}
											// Fill the histogram with the variable n
											//h_exp_n -> Fill(ePix);
											
												}
											}
										}
									}
								}
							}
						}
					}
				}

									cout<<"J "<<J<<endl;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

	//for(int N=1000;N<1001; N++){
		//for(int DC=1000;DC<1001; DC=DC++){
			//for(int A=2500;A<2501; A=A++){
				//for(int B=30;B<31; B=B++){
				int N=1;
				int DC=540;
				int A=3500;
				int B=10;

				// Monte Carlo Data ////////////////////////////////////////////////////
				// Get input files//////////////////////////////////////////////////////


					TFile * f_mc = TFile::Open(("./Grilla DC A=3500 B=10/output_N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)+".root").c_str());

					if (!f_mc->IsOpen()) {std::cerr << "ERROR: cannot open the root file with MC data" << std::endl;}
					TTree * tmc    = (TTree*) f_mc->Get("hitSumm");
					TH1D * h_mc_n      =  new TH1D("h_mc_n", ("Distribucion de tamanio de clusters -- N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)).c_str(), nbins, xmin, xmax);
					h_mc_n -> Sumw2();

					int Entries_mc = tmc -> GetEntries();
					cout<<"Entries in MC file: "<<Entries_mc<<endl;

					Enable_and_Set_Branches(tmc);

				// Get information from trees///////////////////////////////////////////

					for(int i_event=0;i_event<Entries_mc; i_event++){

					tmc->GetEntry(i_event);

						if (e>emin && e<emax){  // number of electrons

							// Check if one of the pixels in the cluster is smaller that minePix
							bool noLowPixInCluster = true;
							for (int p = 0; p < nSavedPix; ++p){
								if(ePix[p]<minePix){
									noLowPixInCluster = false;
									break;
								}
							}

							if (noLowPixInCluster){
									

										// Fill the histogram with the variable n
											h_exp_n -> Fill(e);
					

								    	// Fill the histogram with the variable ePix
											//for (int p = 0; p < nSavedPix; ++p){
											//h_mc_n -> Fill(ePix[p]);
											//}
							}
						}
					}




				// Plot histograms /////////////////////////////////////////////////////
					TCanvas*  clusters   = new TCanvas("clusters","clusters",1200,1200);
					clusters->cd();
					gStyle->SetOptStat(0);
					//gPad->DrawFrame(0.5,0.85,3.5,1.1,"Summary of fit results;;Fit Bias");

					h_exp_n -> SetLineColor(kRed +3);
					h_mc_n -> SetLineColor(kBlue);

					h_exp_n -> SetMarkerColor(kRed +3);
					h_mc_n -> SetMarkerColor(kBlue);

					h_exp_n -> SetMarkerSize(1.0);
					h_mc_n -> SetMarkerSize(1.0);

					h_exp_n -> SetMarkerStyle(24);
					h_mc_n -> SetMarkerStyle(24);

					//int norm=1; //Normalization
					//clusters->SetLogy(1);
					//Double_t scale_exp = norm/h_exp_n->Integral();
					//h_exp_n->Scale(scale_exp);
					h_exp_n->SetMaximum(30000);
					h_exp_n ->Draw("HIST E1");
					//clusters->SaveAs( "./figures/Dist_n_exp.png");

					//Double_t scale_mc = norm/h_mc_n->Integral();
					//h_mc_n->Scale(scale_mc);
					//h_mc_n->SetMaximum(norm);
					h_mc_n ->Draw("HIST E1 same");

					vector<double> EXP(bines);
					vector<double> EXP_error(bines);

					for(int i_bin=1;i_bin<bines; i_bin++){
						EXP[i_bin]= (h_exp_n->GetBinContent(i_bin));
						EXP_error[i_bin] = h_exp_n->GetBinError(i_bin);
						cout<<EXP[i_bin]<<" "<<EXP_error[i_bin]<<endl;
					}

					vector<double> MC(bines);
					vector<double> MC_error(bines);

					for(int i_bin=1;i_bin<bines; i_bin++){
						MC[i_bin]= h_mc_n->GetBinContent(i_bin);
						MC_error[i_bin] = h_mc_n->GetBinError(i_bin);
						cout<<MC[i_bin]<<" "<<MC_error[i_bin]<<endl;
					}

					Double_t chi2=0;
					for(int i_bin=1;i_bin<bines; i_bin++){
						chi2=pow((EXP[i_bin]-MC[i_bin]),2)/MC_error[i_bin];
					}
						cout<<endl;
						cout<<"Chi2: "<<chi2<<endl;
						
						cout<<"Nro de eventos: "<<J<<endl;

					TLegend *lg = new TLegend(0.35,0.7,0.9,0.9,NULL,"brNDC");
					lg->SetBorderSize(1);
					lg->SetTextFont(50);
					lg->SetTextSize(0.03146);
					lg->SetFillColor(10);
					lg->SetFillStyle(1001);
					lg->AddEntry(h_mc_n,("MC N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)).c_str(),"l");
					lg->AddEntry(h_mc_n,("emin="+itos(emin)+"_emax="+itos(emax)+"_minePix="+itos(minePix)+"_emin="+itos(emin)+"_emax="+itos(emax)).c_str(),"l");
					lg->AddEntry(h_mc_n,("Cond: x_yBary ; ePix =! 305"),"l");
					
					/*
					lg->AddEntry(h_exp_n,("Datos -             Chi2: "+itos(chi2)).c_str(),"l");
					lg->AddEntry(h_mc_n,("MC N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)).c_str(),"l");
					lg->AddEntry(h_mc_n,(,"l");
					*/
					//lg->AddEntry(chi2,"Chi cuadrado","l");
					lg->Draw("same");
					
					clusters->SaveAs(("Dist_n_MC_N="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)+"_minePix="+itos(minePix)+"_emin="+itos(emin)+"_emax="+itos(emax)+"_LOG.png").c_str());

					gPad->Update();
					//getchar();
					if (chi2 < 50) {
						//gPad->WaitPrimitive();
					}


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
