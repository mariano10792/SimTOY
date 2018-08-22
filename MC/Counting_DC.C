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
  int emin =1400 ; int emax = 1500;
  int ohdu_numer = 4;
//number of bins to take into account for chi2
  int bines = 7;
  float minePix = 1.5; // will be clasified as 1e-

  ////////////////////////////////////////////////////////////////////////

  int Entries_mc = 1;
  int Entries_exp = 1;
  int xmin = 0; // range for histograms
  int xmax =4; // range for histograms
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



void Counting_DC(){

// Experimental Data ///////////////////////////////////////////////////
// Get input files//////////////////////////////////////////////////////

					TH1D * h_exp_e[267];
					TString histname_tmp_exp;

					for(int ystarbin=0;ystarbin<nbins+1;++ystarbin){
					histname_tmp_exp  = "h_exp_e";
					histname_tmp_exp += ystarbin;
					h_exp_e[ystarbin] = new TH1D(histname_tmp_exp.Data(),"",30,emin,emax);
					h_exp_e[ystarbin]->Sumw2();
					}

					

				// Get information from trees///////////////////////////////////////////

				//vector<double> events_exp(50000);
				//vector<double> mean_exp(nbins+1);
				//vector<double> sigma_exp(nbins+1);
				//vector<double> fano_exp(nbins+1);
				int initfile = 73;
				int nfiles= 267;
				vector<double> ENE(nfiles-initfile); //194 espacios
				vector<double> I(nfiles-initfile); //194 espacios
				
		for(int ene=initfile; ene<nfiles; ene++){
					
					TFile * f_exp = TFile::Open(("../../caca/output_"+itos(ene)+".root").c_str());
					if ((!f_exp) || (f_exp->IsZombie())) { //deletes file if it doesnt find it
					    std::cout << "(no file || file is zombie) skip file: "
					              << f_exp << std::endl;
									ENE[ene-initfile]=ene;
									I[ene-initfile]=0;
					    delete f_exp; continue;
					  }
					
					if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
					TTree * texp = (TTree*) f_exp->Get("hitSumm");



					int Entries_exp = texp -> GetEntries();
					cout<<"Entries in experimental data file "<< ene <<" are "<<Entries_exp<<endl;

					Enable_and_Set_Branches(texp);
					int i = 0;
					for(int i_event=0;i_event<Entries_exp; i_event++){
					texp->GetEntry(i_event);

							if (ohdu == ohdu_numer) {
								
							

      							
      						for (int p = 0; p < nSavedPix; ++p){
									if (xPix[0]<250){
										if(ePix[p]<8 *minePix){
      									i++;
									}
								}
							}
						
				


                 // cout << "i= "<< i << endl;

		                      
						    }
						  }
						 cout<<"DC in experimental data file "<< ene <<" are "<<i<<endl;
						 ENE[ene-initfile]=ene;
						 I[ene-initfile]=i;
						 					  
					    }


					////matrix with charges
					//int **N1 = new int*[1];  int **N2 = new int*[2]; int **N3 = new int*[3]; int **N4 = new int*[4]; int **N5 = new int*[5]; int **N6 = new int*[6]; int **N7 = new int*[7];

					//for(int i = 0; i < cols; ++i) {
					//N1[i] = new int[events[1]]; N2[i] = new int[events[2]]; N3[i] = new int[events[3]]; N4[i] = new int[events[4]]; N5[i] = new int[events[5]]; N6[i] = new int[events[6]]; N7[i] = new int[events[7]];
					//}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

	//for(int N=1000;N<1001; N++){
		//for(int DC=1000;DC<1001; DC=DC++){
			//for(int A=2500;A<2501; A=A++){
				//for(int B=30;B<31; B=B++){
				int N=200;
				int DC=0;
				int A=3500;
				int B=10;

				// Monte Carlo Data ////////////////////////////////////////////////////
				// Get input files//////////////////////////////////////////////////////


					TFile * f_mc = TFile::Open(("output_N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)+".root").c_str());

					if (!f_mc->IsOpen()) {std::cerr << "ERROR: cannot open the root file with MC data" << std::endl;}
					TTree * tmc    = (TTree*) f_mc->Get("hitSumm");
					//TH1D * h_mc_n      =  new TH1D("h_mc_n", ("Distribucion de tamanio de clusters -- N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)).c_str(), nbins, xmin, xmax);
					TH1D * h_mc_e[nbins+1];
					TString histname_tmp;

					for(int J=0;J<nbins+1;++J){
					histname_tmp  = "h_mc_e";
					histname_tmp += J;
					h_mc_e[J] = new TH1D(histname_tmp.Data(),"",30,emin,emax);
					h_mc_e[J]->Sumw2();
					}

					int Entries_mc = tmc -> GetEntries();
					cout<<"Entries in MC file: "<<Entries_mc<<endl;

					Enable_and_Set_Branches(tmc);

				  vector<double> events_mc(50000);
				  vector<double> mean_mc(nbins+1);
				  vector<double> sigma_mc(nbins+1);


				for (int ene=1; ene<nbins+1; ene++){
          int i = 0;
					for(int i_event=0;i_event<Entries_mc; i_event++){

					//texp->GetEntry(i_event);



								if (e>emin && e<emax){  // number of electrons
                  if (n==ene){
									h_mc_e[ene]->Fill(e);

									sigma_mc[ene]= h_mc_e[ene]->GetRMS();
									mean_mc[ene]= h_mc_e[ene]->GetMean();

//									events_mc[i]=e;
									i++;
                  cout << "i= "<< i << endl;






								}
							}
						}
						// aca calculamos varianza y media






					}


					/*
				// Plot histograms /////////////////////////////////////////////////////
					TCanvas*  clusters   = new TCanvas("clusters","clusters",1000,500);
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

					int norm=1; //Normalization
					//clusters->SetLogy(1);
					Double_t scale_exp = norm/h_exp_n->Integral();
					h_exp_n->Scale(scale_exp);
					h_exp_n->SetMaximum(0.5);
					h_exp_n ->Draw("HIST E1");
					//clusters->SaveAs( "./figures/Dist_n_exp.png");

					Double_t scale_mc = norm/h_mc_n->Integral();
					h_mc_n->Scale(scale_mc);
					h_mc_n->SetMaximum(0.5);
					h_mc_n ->Draw("HIST E1 same");
					clusters->SaveAs(("Dist_n_MC_N="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)+".png").c_str());

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

					TLegend *lg = new TLegend(0.6,0.6,0.9,0.8,NULL,"brNDC");
					lg->SetBorderSize(1);
					lg->SetTextFont(62);
					lg->SetTextSize(0.03146);
					lg->SetFillColor(10);
					lg->SetFillStyle(1001);
					lg->AddEntry(h_exp_n,("Datos - Chi2: "+itos(chi2)).c_str(),"l");
					lg->AddEntry(h_mc_n,("MC N0="+itos(N)+"_DC="+itos(DC)+"_A="+itos(A)+"_B="+itos(B)).c_str(),"l");
					//lg->AddEntry(chi2,"Chi cuadrado","l");
					lg->Draw("same");
					*/



				//}
			//}
		//}
	//}
	
	
	
	/////////////
	// save file
	///////////
	ofstream myfile("./report_8e.txt");
	for(int ene=initfile; ene<nfiles; ene++){
	myfile << ENE[ene-initfile] <<"	"<< I[ene-initfile] <<"\n" ;
	
}	
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
