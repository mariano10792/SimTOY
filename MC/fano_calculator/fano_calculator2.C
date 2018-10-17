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
  int emin =1505 ; int emax = 1585;
  int ohdu_numer = 4;
//number of bins to take into account for chi2
  int bines = 100;
  float minePix = 0; // will be clasified as 1e-

  ////////////////////////////////////////////////////////////////////////

  int Entries_mc = 1;
  int Entries_exp = 1;
  //int xmin = 0; // range for histograms
  //int xmax =0; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;

  int nbins = 5;

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



void fano_calculator2(){

// Experimental Data ///////////////////////////////////////////////////
// Get input files//////////////////////////////////////////////////////

cout<<"min ePix: "<< minePix<<endl;
					int M=10;
					
					ofstream myfile;
					myfile.open ("../mejor/figures/example.txt");	
							
					for(int m=1;m<M+1;++m){
					
					
					TFile * f_exp = TFile::Open(("../mejor/merged_"+itos(m)+".root").c_str());
					if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
					TTree * texp = (TTree*) f_exp->Get("hitSumm");
					
					
					// esto no se que hace...
					TH1D * h_exp_e[nbins+1];
					TString histname_tmp_exp;
					
					//for each n
					for(int ystarbin=0;ystarbin<nbins+1;++ystarbin){
					histname_tmp_exp  = "h_exp_e";
					histname_tmp_exp += ystarbin;
					h_exp_e[ystarbin] = new TH1D(histname_tmp_exp.Data(),"",bines,emin,emax);
					h_exp_e[ystarbin]->Sumw2();
					}

					int Entries_exp = texp -> GetEntries();
					cout<<"Entries in experimental data file: "<<Entries_exp<<endl;

					Enable_and_Set_Branches(texp);

				// Get information from trees///////////////////////////////////////////

				vector<double> mean_exp_fit(nbins+1);
				vector<double> sigma_exp_fit(nbins+1);
				vector<double> fano_exp_fit(nbins+1);
				
				vector<double> mean_exp_fit_error(nbins+1);
				vector<double> sigma_exp_fit_error(nbins+1);
				vector<double> fano_exp_fit_error(nbins+1);
				
				vector<double> events_exp_fit(nbins+1);


				for (int ene=1; ene<nbins+1; ene++){
					int i = 0;
					for(int i_event=0;i_event<Entries_exp; i_event++){
					texp->GetEntry(i_event);

				//			if (ohdu == ohdu_numer) {
								if (e>emin && e<emax){  // number of electrons
                  if (n==ene){

                    // Check if one of the pixels in the cluster is smaller that minePix
      							bool noLowPixInCluster = true;
      							for (int p = 0; p < nSavedPix; ++p){
      								if(ePix[p]<minePix){
      									noLowPixInCluster = false;
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

      							if (noLowPixInCluster){
      								if (noBadTransferInCluster){
      								if (xBary>250 && xBary<490){
      									if (yBary>3 && yBary<48){

                 // cout << "i= "<< i << endl;
									h_exp_e[ene]->Fill(e);
									h_exp_e_total->Fill(e);

									events_exp[i]=e;
									i++;
                          }
                        }
                      }
                    }
									}
								}
					//		}
						}
						// aca calculamos varianza y media

						h_exp_e[ene]->Fit("gaus");
						TF1 *fit1 = (TF1*)h_exp_e[ene]->GetFunction("gaus");
						
						TCanvas*  canvitas   = new TCanvas("clusters","clusters",800,500);
						canvitas->cd();
						fit1->Draw("HIST E1");
						h_exp_e[ene]->Draw("HIST E1 same");
						canvitas->SaveAs(("../mejor/figures/Fit_M="+itos(m)+"_N="+itos(ene)+".png").c_str());

												
						sigma_exp_fit[ene]= fit1->GetParameter(2);
						sigma_exp_fit_error[ene]= fit1->GetParError(2);
						mean_exp_fit[ene]= fit1->GetParameter(1);
						mean_exp_fit_error[ene]= fit1->GetParError(1);
						fano_exp_fit[ene]= pow(sigma_exp_fit[ene],2)/mean_exp_fit[ene];
						fano_exp_fit_error[ene]= pow(pow(2*sigma_exp_fit[ene]*sigma_exp_fit_error[ene]/mean_exp_fit[ene],2)+pow(mean_exp_fit_error[ene]*sigma_exp_fit[ene]*sigma_exp_fit[ene]/(mean_exp_fit[ene]*mean_exp_fit[ene]),2),0.5);
						events_exp_fit[ene]=h_exp_e[ene]->GetEntries();
					
						

					}
						// save to file
						myfile << "Mean" << "	" <<"Sigma" << "	" <<"Fano factor" << "	" <<"Fano error" << "	" << " # Events for M="<< m <<", increasing n" <<endl;
						for(int i = 1; i < nbins+1; i ++) {
						myfile  << mean_exp_fit[i]  << "	" << sigma_exp_fit[i] << "	" << fano_exp_fit[i] << "	" << fano_exp_fit_error[i]  << "	" <<  events_exp_fit[i] << endl;
						}			
						
					
	} //closes for 
	myfile.close();
						
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
