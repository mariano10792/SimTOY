#include <sstream>
#include <iomanip>
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
  int emin =1485 ; int emax = 6050;
  int ohdu_numer = 4;
//number of bins to take into account for chi2
  int bines = 2000;
  float minePix = 0; // will be clasified as 1e-

  ////////////////////////////////////////////////////////////////////////

  int Entries_mc = 1;
  int Entries_exp = 1;
  //int xmin = 0; // range for histograms
  //int xmax =0; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;

  int nbins = 6;

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



////////////////////////////////////////////////////////////////////////
// Cluster follower selects the 1828 clusters of the 1.fits, 6.fits etc 121.fits and makes a list.
// Later, it compares that list of clusters with the 1_2.fits, 6_7.fits etc and keeps only the clusters that exited in the original list.
// By doing this, we can evaluate the effect of superoposition of auger electrons with photoelectrons and with other auger electrons.
////////////////////////////////////////////////////////////////////////


//This version, on the other hand compares 1.fits with 1_2, 1_2_3 and 1_2_3_4 but doesnt jump to 6, it goes to 2, and then 3 allowing a
// mayor base of data but paying the correlation of this. It's complementary of cluster_follower

void cluster_follower2(){



cout<<"min ePix: "<< minePix<<endl;
					int M=5; //this will be done in order to look initially at fits 1+M*k k>=0 (1,6,11,16,etc)
          int I=0; //this integer is used for the fake run
          int K=100; //  because 24*M+1=121 and we got 127 fits (with this we go til 126).
          int i=0; //length of the m=1 events list (list_m1)

					//ofstream myfile;
					//myfile.open ("../mejor/figures/example.txt");


          TH1D * h_exp_e[M];
					TString histname_tmp_exp;

	        for(int ystarbin=0;ystarbin<M;++ystarbin){
					histname_tmp_exp  = "h_exp_e";
					histname_tmp_exp += ystarbin;
					h_exp_e[ystarbin] = new TH1D(histname_tmp_exp.Data(),"",bines,emin,emax);
					h_exp_e[ystarbin]->Sumw2();
					}

//As said, this part of the code creates a list of clusters with no merging (1.fits, 6.fits, etc).

          int of_list_m1[K+1];
          int list_m1[50000][2];
          int l=1;
          bool mybool =true; //this bool is need for the l increment loop to function.
               int j=0; //position in the list_m1 list
            for(int k=0;k<K+1;++k){

            if (mybool==false){l=l+1;}

          // open file and set TTree
          TFile * f_exp = TFile::Open(("../mejor/grupos_de_1/merge_"+itos(l)+"_"+itos(l)+".root").c_str());
          if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
          TTree * texp = (TTree*) f_exp->Get("hitSumm");
          int Entries_exp = texp -> GetEntries();
          cout<<"Entries in experimental data file: "<<Entries_exp<<endl;
          Enable_and_Set_Branches(texp);


          cout<<"l = "<<l<<endl;


          // Get information from trees///////////////////////////////////////////

          int t=0;

          for (int ene=0; ene<nbins+1; ene++){

          if (ene==0) {continue;}
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
                              t++;

                                i++;

                                of_list_m1[k]=i;

                                list_m1[j][0]=xPix[0];
                                list_m1[j][1]=yPix[0];
                                h_exp_e[0]->Fill(e);
                                j++;


                                // if (k==0) {
                                //   of_list_m1[k]=i;
                                // }else{
                                //   of_list_m1[k]=i;
                                // }


                              // if (m==2){
                              //   list_m1[j][0]=xPix[0];
                              //   list_m1[j][1]=yPix[0];
                              //   h_exp_e[0]->Fill(e);
                              //   j++;
                              //
                              //
                              // }
                              // cout << "t= " <<t<<endl;

                          }
                        }
                      }
                    }
                  }
                }
          	   }
              }
              if (l==1){mybool = false;}
             }

//the final state of i is the number of elements in the m1 cluster list
cout << i << endl;

// of_list_m1 create a CUMULATIVE list of elements (64,128,228,etc) and of_list_m1_bis creates a list of the
// number of elements in each fit (64,64,100,etc). Both are used in the compare loop downwards

int of_list_m1_bis[K+1];
for (int i = 0; i <K+1; i++) {
  if (i==0) {
    of_list_m1_bis[i]=of_list_m1[i];
  }else{
    of_list_m1_bis[i]=of_list_m1[i]-of_list_m1[i-1];
  }
}

//check the list

for (size_t j = 0; j < 64; j++) {
  cout << "list " << list_m1[j][1] << endl;
}




//this loop compares and saves in the TH1D.



            for(int m=2;m<5;++m){
              int l=1;
              bool mybool =true;
              int j=0;

              for(int k=0;k<K+1;++k){

              if (mybool==false){l=l+1;}

            // open file and set TTree
            TFile * f_exp = TFile::Open(("../mejor/grupos_de_"+itos(m+1)+"/merge_"+itos(l)+"_"+itos(l+m)+".root").c_str());
            if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
            TTree * texp = (TTree*) f_exp->Get("hitSumm");
            cout<<"l = "<<l<<endl;
            int Entries_exp = texp -> GetEntries();
            cout<<"Entries in experimental data file: "<<Entries_exp<<endl;

            Enable_and_Set_Branches(texp);

            // Get information from trees///////////////////////////////////////////


            int t=0;
            for (int ene=0; ene<nbins+1; ene++){

            if (ene==0) {continue;}
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
                            t++;


                              for (int p = 0; p < ene; ++p){

                                  int x_resta;
                                  int y_resta;

                                  //this loop
                                  for (int i=0; i<of_list_m1_bis[k];i++){

                                  //k=0 is the first fits so it needs a different code for the x(y)_resta calculation. No big deal.
                                  if (k==0) {
                                    x_resta=xPix[p]-list_m1[i][0];
                                    y_resta=yPix[p]-list_m1[i][1];
                                    }else{
                                      x_resta=xPix[p]-list_m1[of_list_m1[k-1]+i][0];
                                      y_resta=yPix[p]-list_m1[of_list_m1[k-1]+i][1];
                                  }


                          //          We found it! save it.
                                      if (x_resta==0 && y_resta==0){
                                      h_exp_e[m-1]->Fill(e);

//                                      cout << "m= " <<m << endl;
                                    //  cout << t << endl;
                                      break;
                                    }else{

                                      continue;}

                                  }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
            	   }
                }
                if (l==1){mybool = false;}
               }
              }


TCanvas *c1 = new TCanvas("c4","",200,10,1700,1000);
//c1->SetGrid();
c1->Divide(2,2);

for (size_t i = 0; i < 4; i++) {
//h_exp_e[i]->SetMarkerStyle(0);
//h_exp_e[i]->SetMarkerSize(2);
//h_exp_e[i] -> SetLineColor(kRed);
}

c1->cd(1);
h_exp_e[0]->Draw("HIST E1");
c1->cd(2);
h_exp_e[1]->Draw("HIST E1 same");
c1->cd(3);
h_exp_e[2]->Draw("HIST E1 same");
c1->cd(4);
h_exp_e[3]->Draw("HIST E1 same");








/*


            ///////////
						// aca calculamos varianza y media
            //////////


            TF1 *fit1 = new TF1("fit1","gaus",1510,1580);
            h_exp_e[ene]->Fit("fit1","Quiet","R");
						//h_exp_e[ene]->Fit("gaus","Quiet",1515,1575);
						//TF1 *fit1 = (TF1*)h_exp_e[ene]->GetFunction("gaus");

            //////////
            /// Use canvitas to save every fit.
            //////////

						// TCanvas*  canvitas   = new TCanvas("clusters","clusters",800,500);
						// canvitas->cd();
            // h_exp_e[ene]->Draw("HIST E1");
            // fit1->Draw("same");
            //
						// canvitas->SaveAs(("../mejor/figures/Fit_M="+itos(m)+"_N="+itos(ene)+"Fano4.png").c_str());


						sigma_exp_fit[ene][m]= fit1->GetParameter(2);
						sigma_exp_fit_error[ene][m]= fit1->GetParError(2);
						mean_exp_fit[ene][m]= fit1->GetParameter(1);
						mean_exp_fit_error[ene][m]= fit1->GetParError(1);
						fano_exp_fit[ene][m]= pow(sigma_exp_fit[ene][m],2)/mean_exp_fit[ene][m];
						fano_exp_fit_error[ene][m]= pow(pow(2*sigma_exp_fit[ene][m]*sigma_exp_fit_error[ene][m]/mean_exp_fit[ene][m],2)+pow(mean_exp_fit_error[ene][m]*sigma_exp_fit[ene][m]*sigma_exp_fit[ene][m]/(mean_exp_fit[ene][m]*mean_exp_fit[ene][m]),2),0.5);
						events_exp_fit[ene][m]=h_exp_e[ene]->GetEntries();
						events_exp_fit_error[ene][m]=0;

            if (ene==nbins){

              TF1 *fit2 = new TF1("fit2","gaus",1505,1585);
              h_exp_e_total->Fit("fit2","Quiet","R");

              sigma_exp_fit[ene][m]= fit2->GetParameter(2);
  						sigma_exp_fit_error[ene][m]= fit2->GetParError(2);
  						mean_exp_fit[ene][m]= fit2->GetParameter(1);
  						mean_exp_fit_error[ene][m]= fit2->GetParError(1);
  						fano_exp_fit[ene][m]= pow(sigma_exp_fit[ene][m],2)/mean_exp_fit[ene][m];
  						fano_exp_fit_error[ene][m]= pow(pow(2*sigma_exp_fit[ene][m]*sigma_exp_fit_error[ene][m]/mean_exp_fit[ene][m],2)+pow(mean_exp_fit_error[ene][m]*sigma_exp_fit[ene][m]*sigma_exp_fit[ene][m]/(mean_exp_fit[ene][m]*mean_exp_fit[ene][m]),2),0.5);
  						events_exp_fit[ene][m]=h_exp_e_total->GetEntries(); //useless
   						events_exp_fit_error[ene][m]=0; //useless



            }

					}
						// save to file --> dont' forget to uncomment upwards
					//	myfile << "Mean" << "	" <<"Sigma" << "	" <<"Fano factor" << "	" <<"Fano error" << "	" << " # Events for M="<< m <<", increasing n" <<endl;
					//	for(int i = 1; i < nbins+1; i ++) {
					//	myfile  << mean_exp_fit[i][m]  << "	" << sigma_exp_fit[i][m] << "	" << fano_exp_fit[i][m] << "	" << fano_exp_fit_error[i][m]  << "	" <<  events_exp_fit[i][m] << endl;

            if (m==1){
            m=0;
            I++; //
          } //this creates the run fake. To function it needs another instance above that only initiates with I=0. This instance must just count the number of events in the fits 1+M*k k>=0
        }
      }


	} //closes for
	myfile.close();




  TCanvas *c4 = new TCanvas("c4","",200,10,1700,1000);
  c4->SetFillColor(0);
  c4->SetGrid();
  c4->Divide(3,2);

for (int j=0; j<nbins;j++){


  double x_er[M];
  double x[M];
  double y[M];
  double y_er[M];
  for (Int_t i=0;i<M;i++) {
   x[i] = i+1;
   x_er[i] =0;
   y[i] = events_exp_fit[j+1][i+1];
   y_er[i] = events_exp_fit_error[j+1][i+1];
  }
   c4->cd(j+1);

   TGraph *gr = new TGraph(M,x,y);

   gr->SetMinimum(0);
   gr->SetMaximum(3000);
   gr->SetLineColor(j+1);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle(("Eventos para n="+itos(j+1)+"").c_str());
   gr->GetXaxis()->SetTitle("M");
   gr->GetYaxis()->SetTitle("Eventos");
   gr->Draw("AP");

if (j==nbins-1) {gr->SetTitle(("Eventos para n =< "+itos(j)+"").c_str());}

}




c4->SaveAs(("../mejor/figures/EventosGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.png").c_str());
c4->SaveAs(("../mejor/figures/EventosGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.root").c_str());


		TCanvas *c3 = new TCanvas("c3","",200,10,1700,1000);
		c3->SetFillColor(42);
		c3->SetGrid();
		c3->Divide(3,2);

	for (int j=0; j<nbins;j++){


    double x_er[M];
		double x[M];
		double y[M];
		double y_er[M];
		for (Int_t i=0;i<M;i++) {
		 x[i] = i+1;
		 x_er[i] =0;
		 y[i] = mean_exp_fit[j+1][i+1];
		 y_er[i] = mean_exp_fit_error[j+1][i+1]/2;
		}
	   c3->cd(j+1);

	   double promedio = TMath::Mean(M,&y[0]);
	   double RMS = TMath::RMS(M,&y[0]);
	   double prom[M];



	   cout << "promedio para n =" << j+1 <<"	" << promedio <<"+-" << RMS << endl;

	   //TGraph *gr = new TGraph(M,x,y);
	   TGraphErrors *gr = new TGraphErrors(M,x,y,x_er,y_er);


	   auto promedio_str = std::to_string(promedio);
	   auto rms_str = std::to_string(RMS);

	   gr->SetMinimum(1540);
	   gr->SetMaximum(1560);
	   gr->SetLineColor(j+1);
	   gr->SetLineWidth(4);
	   gr->SetMarkerColor(4);
	   gr->SetMarkerStyle(21);
	   gr->SetTitle(("n="+itos(j+1)+"_prom="+promedio_str+"+-"+rms_str+"").c_str());
	   gr->GetXaxis()->SetTitle("M");
	   gr->GetYaxis()->SetTitle("Mean|");
		gr->Draw("AP");
    if (j==nbins-1) {gr->SetTitle(("n =< "+itos(j)+"_prom="+promedio_str+"+-"+rms_str+"").c_str());}

	}

	c3->SaveAs(("../mejor/figures/MeanGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.png").c_str());
  c3->SaveAs(("../mejor/figures/MeanGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.root").c_str());


		TCanvas *c2 = new TCanvas("c2","",200,10,1700,1000);
		c2->SetFillColor(42);
		c2->SetGrid();
		c2->Divide(3,2);

	for (int j=0; j<nbins;j++){


		double x[M];
		double y[M];
		double x_er[M];
		double y_er[M];
		for (Int_t i=0;i<M;i++) {
		 x[i] = i+1;
		 x_er[i] =0;
		 y[i] = sigma_exp_fit[j+1][i+1];
		 y_er[i] = sigma_exp_fit_error[j+1][i+1]/2;
		}
	   c2->cd(j+1);


	   //TGraph *gr = new TGraph(M,x,y);
	   TGraphErrors *gr = new TGraphErrors(M,x,y,x_er,y_er);

	   gr->SetMinimum(10);
	   gr->SetMaximum(20);
	   gr->SetLineColor(j+1);
	   gr->SetLineWidth(4);
	   gr->SetMarkerColor(4);
	   gr->SetMarkerStyle(21);
	   gr->SetTitle(("n="+itos(j+1)+"").c_str());
	   gr->GetXaxis()->SetTitle("M");
	   gr->GetYaxis()->SetTitle("Sigma");
	   gr->Draw("AP");
     if (j==nbins-1) {gr->SetTitle(("n =< "+itos(j)+"").c_str());}
	}


		c2->SaveAs(("../mejor/figures/SigmaGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.png").c_str());
    c2->SaveAs(("../mejor/figures/SigmaGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.root").c_str());









		TCanvas *c1 = new TCanvas("c1","",200,10,1700,1000);
		c1->SetFillColor(0);
		c1->SetGrid();
		c1->Divide(3,2);

	for (int j=0; j<nbins;j++){


		double x[M];
		double y[M];
		double x_er[M];
		double y_er[M];
		for (Int_t i=0;i<M;i++) {
		 x[i] = i+1;
		 x_er[i] =0;
		 y[i] = fano_exp_fit[j+1][i+1];
		 y_er[i] = fano_exp_fit_error[j+1][i+1]/2;
		}
	   c1->cd(j+1);

	   double promedio = TMath::Mean(M,&y[0]);
     double promedio_rounded = (int)(promedio * 1000.0)/1000.0;
	   double RMS = TMath::RMS(M,&y[0]);
     double RMS_rounded = (int)(RMS * 1000.0)/1000.0;
	   double prom[M];


	   cout << "promedio para n =" << j+1 <<"	" << promedio_rounded <<"+-" << RMS_rounded << endl;

	   //TGraph *gr = new TGraph(M,x,y);
	   TGraphErrors *gr = new TGraphErrors(M,x,y,x_er,y_er);

     //if(j==0) {cout << "error  fano para n =" << j+1 <<"	" << y_er[0] << endl;}


	   string promedio_str = std::to_string(promedio_rounded);
	   string rms_str = std::to_string(RMS_rounded);



	   gr->SetMinimum(0.095);
	   gr->SetMaximum(0.2);
	   gr->SetLineColor(j+1);
	   gr->SetLineWidth(4);
	   gr->SetMarkerColor(1);
	   gr->SetMarkerStyle(21);
	   gr->SetTitle(("n="+itos(j+1)+" -> Promedio = "+promedio_str+"+-"+rms_str+"").c_str());
	   gr->GetXaxis()->SetTitle("M");
	   gr->GetYaxis()->SetTitle("Fano Factor");
		 gr->Draw("AP");
     if (j==nbins-1) {gr->SetTitle(("n =< "+itos(j)+" \n Promedio="+promedio_str+"+-"+rms_str+"").c_str());}

	}


	c1->SaveAs(("../mejor/figures/FanoGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.png").c_str());
  c1->SaveAs(("../mejor/figures/FanoGraph_N="+itos(nbins)+"_hasta M="+itos(M)+"_minePix="+itos(minePix)+"Fano4.root").c_str());

		// TCanvas *c0 = new TCanvas("c0","",200,10,1600,1000);
		// c0->cd();
		// h_exp_e_total->Draw();



*/

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



/*
void graph() {
   //Draw a simple graph
   // To see the output of this macro, click begin_html <a href="gif/graph.gif">here</a>. end_html
   //Author: Rene Brun

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

   c1->SetFillColor(42);
   c1->SetGrid();

   const Int_t n = 10;
   Double_t x[n], y[n];
   for (Int_t i=0;i<M;i++) {
     x[i] = i*+1;
     y[i] = fano_exp_fit[1][i+1];
     printf(" i %i %f %f \n",i,x[i],y[i]);
   }

   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("a simple graph");
   gr->GetXaxis()->SetTitle("X title");
   gr->GetYaxis()->SetTitle("Y title");
   gr->Draw("ACP");

   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
*/
