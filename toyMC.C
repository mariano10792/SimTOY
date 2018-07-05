// Version 2/07/2018 - 15:30 hs

// En esta version el parametro B no juega ningun rol.
/*
Código escrito para simular:

- Coordenada (xx,yy) donde tuvo lugar la interacción del foton X
- Profundidad zz donde tuvo lugar la interacción del foton X
- Número de cargas generadas en cada interacción
- Distribución de cargas que llega a la superficie del CCD
- Numero de cargas por pixel.
- Pixeles generador por corriente oscura
- Distribucion de los pixeles que conforman cada cluster

*/

#include <random>
#include <ctime>
#include <cstdlib>
#include "fitsio.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <chrono>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
using namespace std;


// Internal Parameters /////////////////////////////////////////////////
// CCD size. Real size: 4126 x 866.
int nx = 250;           // Number of pixels in x-direction
int ny = 50;            // Number of pixels in y-direction
int pixSize= 15;        // Pixel size side in microns
int xSize = nx*pixSize; // x CCD size in microns
int ySize = ny*pixSize; // y CCD size in microns
int sizeArray =nx*ny;   // Total number of pixels in the CCD

double tau_Si = 58;     // reference value 58 um

long fpixel[2] = {1,1};

// This loop runs over each X-ray interaction //////////////////////////
// void interaction(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7,  TH2F *h2p_int,  TH2F *h2p_TOTAL){
void interaction(TH2F *h2p_int,  TH2F *h2p_TOTAL, int N0, float A, float B, TH1D *h5){

vector<double> emeannumber(N0); //emeannumber: mean number of electrons generates by a X-ray //////////
vector<double> tau(N0); // skin depth

for (int i = 0; i < N0; ++i){
//cout<<i<<"  \n ";

TRandom3 kalfa(0);
double alfa= kalfa.Uniform(0,1.92);

//cout<<alfa<<"  \n ";

// Cuidado que emeannumber no es la energia sino el numero medio de
// electrone que produce cada X. Se parecen solo porque la ganancia y
// la energia de ionizacion son similares (3.7) y se cancelan mutuamente.

if (alfa<0.51){
	emeannumber[i]=5887.65; // [eV]
	tau[i]=tau_Si;
	//i++;
	//cout<<alfa<<" A  \n ";
}
if (0.51<alfa && alfa<1.51){
	emeannumber[i]=5898.75; // [eV]
	tau[i]=tau_Si;
	//i++;
	//cout<<alfa<<" B  \n ";
}
if (1.51<alfa && alfa<1.715){
	emeannumber[i]=6490.45; // [eV]
	tau[i]=tau_Si*(6.5/5.9); // a first order correction with energy ...
	//i++;
	//cout<<alfa<<" C  \n ";
}
if (1.715<alfa && alfa<1.92){
	//i++;
	emeannumber[i]=6535.2; // [eV]
	tau[i]=tau_Si*(6.5/5.9); // a first order correction with energy ...
	//cout<<alfa<<" D  \n ";
}
//cout<<i<<"  \n ";
}


for (int j = 0; j < N0; ++j){

	// Variables for random numbers ////////////////////////////////////////
	vector<double> xx(N0);
	vector<double> yy(N0);
	vector<double> zz(N0);

/*
	Para yN = 250 µm, ND < 6 × 1011 cm−3, y para yN = 675 µm, ND < 1;9 × 1011 cm−3 donde
	donde yN es el espesor de la CCD e ND es el dopaje de dadores (Miguel)
	permitividad relativa silicio 11.68 (wikipedia)
	permitivdad silicio 11.68 * 8.85×10^−12 F/m (wikipedia)
	carga electron  −1,6× 10^−19 C (also wikipedia)

	a1= (( 6 × 10^5 m−3 * −1,6× 10^−19 C ))  / (( 11.68 * 8.85×10^−12 F/m  ))	 (para yN=250um)
	unidades de a1= V/m^2

	diffusion coeficient silicon = k T µe / q  donde µe(uh) es la movilidad de los electrones (huecos)
	unidades de difusion m^2/t
	k= 1.38064852(79)×10−23 J K^-1
	T= 120 K
	q= carga electron = −1,6× 10^−19 C (also wikipedia)
	µh = 1,35 × 108 (120K)−2,20 * cm^2 (V s)^-1 (miguel fuente 51)

	D= (( 1.38064852(79)×10−23 J K^-1 * 120 K * 1,35 × 108 (120K)−2,20 * cm^2 (V s)^-1  ) /  1,6× 10^−19 C  ))
	nu es efectivamente uh

	E(y_w) habria que relativamente estimarlo con un V fijo.

	vector<double> AA(N0); //constante 1 formula sigma cuadrado = 2D/a1*nu
	vector<double> BB(N0); //constante 2 formula sigma cuadrado a1/E(y_w)

*/

	vector<double> electrons(N0);
    vector<double> sigma(N0);

	// Generate a pair (xx, yy) of uniform random numbers //////////////
	// (xx,yy) represents the coordinates of interaction in the CCD ////

    TRandom3 rx(0); //  seed=0  ->  different numbers every time
    //h1->Fill(xx[j]=rx.Uniform(0,xSize));
    xx[j]=rx.Uniform(0,xSize);

    TRandom3 ry(0); //  seed=0  ->  different numbers every time
    //h2->Fill(yy[j]=ry.Uniform(0,ySize));
    yy[j]=ry.Uniform(0,ySize);

	////////////////////////////////////////////////////////////////////
	// Generate zz: depth in the CCD
	// An exponential random number with coefficient tau

    TRandom3 rz(0); //  seed=0  ->  different numbers every time
    //h3->Fill(zz[j] = rz.Exp(tau));
	zz[j] = rz.Exp(tau[j]);
	if (zz[j]>200){
		continue;
	}
	//cout << "zeta = "<< zz[j] << endl;

	//cout << endl;
	//cout << "x = "<< xx[j] << " , " << "y = "<< yy[j] << " , " << "z = "<< zz[j] << endl;

    ////////////////////////////////////////////////////////////////////
    // emeannumber is the mean number of electrons generated by the primary electron.
    // Primary electron here is who carries the energy released
    // by the x-ray after photoelectric interaction.

	// Generate a Poisson random number with mu = emeannumber //////////

	// Fe55 decay


	TRandom3 re(0); //seed=0  ->  different numbers every time
	//h4->Fill(electrons[j] = re.Poisson(emeannumber[j]));
	electrons[j] = re.Poisson(emeannumber[j]);

	//cout << "#e = "<< electrons[j] << endl;













	////////////////////////////////////////////////////////////////////
    // sigma of charge distribution on the CCD surface in microns
    // proportional to the square root of zz depth
    // AA y BB variables from Moroni 2015.

    h5->Fill(sigma[j] = pow(-A/100*log((-B/10000)*zz[j]+1),0.5));
	//cout << "sigma = "<< sigma[j] << endl;
    // Reference:
    //sigma[j] = pow(-A*log(abs((B/10000)*zz[j]-1)),0.5);

















    //cout << "sigma = "<< sigma[j] << endl;

    ////////////////////////////////////////////////////////////////////
    // charge distribution on the CCD surface
    // It generates electrons[j] random numbers with
    // bivariate multinormal distribution

    // Coordinates x and y of each charge colected /////////////////////
    vector<double> chargex(electrons[j]);
    vector<double> chargey(electrons[j]);

    for (int k = 0; k < electrons[j]; ++k) {

		////////////////////////////////////////////////////////////////
		// rcx is the number of excited electrons after difussion
		// Generate a gaussian random number with mean xx[j]
		TRandom3 rcx(0); //seed=0  ->  different numbers every time
		//h6->Fill(chargex[i] = rcx.Gaus(xx[j],sigma[j]));
		chargex[k] = rcx.Gaus(xx[j],sigma[j]);
		//cout << endl;
		//cout << "x on CCD = "<< chargex[i] << endl;

		////////////////////////////////////////////////////////////////
		// rcy is the number of excited electrons after difussion
		// Generate a gaussian random number with mean yy[j]
		TRandom3 rcy(0); //seed=0  ->  different numbers every time
		//h7->Fill(chargey[i] = rcy.Gaus(yy[j],sigma[j]));
		chargey[k] = rcy.Gaus(yy[j],sigma[j]);
		//cout << "y on CCD = "<< chargey[j] << endl;
	}

for (int k = 0; k < electrons[j]; ++k){

	 h2p_int->Fill(chargex[k], chargey[k]); //Only interactions
	 h2p_TOTAL->Fill(chargex[k], chargey[k]); // Only interactions, up to now ...

}
	//cout << "#e = "<< j << endl;
}   // End loop over x-rays interactions



}

// Add dark Current ////////////////////////////////////////////////////
void addDC(TH2F *h2p_DC, TH2F *h2p_TOTAL, int darkC){

	//for (int i = 0; i < darkC; ++i){
		//TRandom3 XDC(0); //  seed=0  ->  different numbers every time
		//int xdc=XDC.Uniform(0,xSize);
		//TRandom3 YDC(0); //  seed=0  ->  different numbers every time
		//int ydc=YDC.Uniform(0,ySize);

	    //double xdc =  rand() % xSize;
        //double ydc =  rand() % ySize;
        //h2p_DC->Fill(xdc,ydc);
        //h2p_TOTAL->Fill(xdc,ydc);
	//}

		vector<int> xdc(darkC);
		vector<int> ydc(darkC);

		for (int i = 0; i < darkC; ++i){

		TRandom3 YDC(0); //  seed=0  ->  different numbers every time
		ydc[i]=YDC.Uniform(0,ySize);
	}

		for (int i = 0; i < darkC; ++i){

		TRandom3 XDC(0); //  seed=0  ->  different numbers every time
		xdc[i]=XDC.Uniform(0,xSize);

	}

		for (int i = 0; i < darkC; ++i){

		h2p_DC->Fill(xdc[i],ydc[i]);
        h2p_TOTAL->Fill(xdc[i],ydc[i]);

	}
}

// Save the content of each histogram into pix variables ///////////////
void pix(double*  pix_int, double*  pix_dc, double*  pix_total, TH2F *h2p_int, TH2F *h2p_DC, TH2F *h2p_TOTAL){
for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {

		   pix_int[j*nx+i]= h2p_int->GetBinContent(i+1,j+1);
		   //cout<<"bin: "<<i*nx+j<<" Interactions: "<<pix_int[i*nx+j]<<endl;

		   pix_dc[j*nx+i]= h2p_DC->GetBinContent(i+1,j+1);
		   //cout<<"bin: "<<i*nx+j<<" Dark Current: "<<pix_dc[i*nx+j]<<endl;

		   pix_total[j*nx+i]= h2p_TOTAL->GetBinContent(i+1,j+1);
		   //cout<<"bin: "<<i*nx+j<<" Interac + DC: "<<pix_total[i*nx+j]<<endl;
		   //cout<<endl;
		}
}
}

// Save simulation as fits files
void save(int nx, int ny, int SizeArray, double*  pix_int, double*  pix_dc, double*  pix_total, int N0, int darkC, int A, int B, int R){

cout<< "Starting to save fits files ..."<<endl;

// Save histograms as .fits files //////////////////////////////////////

	int naxis = 2;
    int status = 0;
    long naxesOut[9] = {nx, ny, 1, 1, 1, 1, 1, 1, 1}; //va desde 0 hasta donde marques

///////////////////////Real Interactions ///////////////////////////////
/*
	std::string outMeanFitsFile1 = "interactions_A=";
    fitsfile *outClusterptr1;
    fits_create_file(&outClusterptr1, (outMeanFitsFile1+ std::to_string(A)+"_B="+ std::to_string(B)+".fits").c_str(), &status);
	fits_create_img(outClusterptr1, -32, naxis, naxesOut, &status);
    fits_write_pix(outClusterptr1, TDOUBLE, fpixel, sizeArray, pix_int, &status);
    fits_close_file(outClusterptr1,  &status);
    cout<< "real_interactions.fits saved "<<endl;

////////////////////////// Dark Current ////////////////////////////////

	std::string outMeanFitsFile2 = "i_dark_current_A=";
    fitsfile *outClusterptr2;
    fits_create_file(&outClusterptr2, (outMeanFitsFile2+ std::to_string(A)+"_B="+ std::to_string(B)+".fits").c_str(), &status);
	fits_create_img(outClusterptr2, -32, naxis, naxesOut, &status);
    fits_write_pix(outClusterptr2, TDOUBLE, fpixel, sizeArray, pix_dc, &status);
    fits_close_file(outClusterptr2,  &status);
    cout<< "dc.fits saved "<<endl;

*/
///////////////  Real Interactions + Dark Current///////////////////////

	std::string outMeanFitsFile3 = "./MC/MC_N0=";
    fitsfile *outClusterptr3;
    fits_create_file(&outClusterptr3, (outMeanFitsFile3+std::to_string(N0)+"_DC="+std::to_string(darkC)+"_A="+ std::to_string(A)+"_B="+ std::to_string(B)+"_R="+ std::to_string(R)+".fits").c_str(), &status);
	fits_create_img(outClusterptr3, -32, naxis, naxesOut, &status);
    fits_write_pix(outClusterptr3, TDOUBLE, fpixel, sizeArray, pix_total, &status);
    fits_close_file(outClusterptr3,  &status);
    cout<< "real_interactions+dc.fits saved "<<endl;

////////////////////////////////////////////////////////////////////////
}

// Plot histograms /////////////////////////////////////////////////////
/*void hist(TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7){
if (N0>99){
	TCanvas* canvas = new TCanvas("canvas","canvas",10,10,700,800);
	canvas->Divide(3,3);

	canvas->cd(1);
	h1->Draw(""); // xx

	canvas->cd(2);
	h2->Draw(""); // yy

	canvas->cd(3);
	h3->Draw(""); // zz

	canvas->cd(4);
	gStyle->SetOptStat(10);
	h4->Draw(""); // #e

	canvas->cd(5);
	h5->Draw(""); // sigma

	canvas->cd(6);
	h6->Draw(""); // x

	canvas->cd(7);
	h7->Draw(""); // y
}


}*/

int main(int argc, char* argv[]){

TStopwatch t; // time counter
t.Start();

int N0 = atoi(argv[1]);    // N0 represents the number of X-rays that interact with the CCD ////////
int darkC = atoi(argv[2]); // Dark Current total events

float A = atoi(argv[3]); // first parameter to fit (2D/a1*nu)
float B = atoi(argv[4]); // second parameter to fit (a1/E(y_w))


int R = atoi(argv[5]); // RUN number // only to label simulations with the same parameters

double*  pix_int = new double[sizeArray];   // pixels array with real interactions
double*  pix_dc = new double[sizeArray];    // pixels array with dark current
double*  pix_total = new double[sizeArray]; // pixels array with all information

cout<<"Variables definition done"<<endl;

// Initialization
for (int i = 0; i < sizeArray; ++i) pix_int[i]=0;
for (int i = 0; i < sizeArray; ++i) pix_dc[i]=0;
for (int i = 0; i < sizeArray; ++i) pix_total[i]=0;

cout<<"Pix variables initialization done"<<endl;

 //Removed in order to include emeannumber into interaction /////////////////
//// Book 1D-histograms //////////////////////////////////////////////////
	//int binnumber = N0/5;

	//auto h1 = new TH1D("h1","x-coordinate",binnumber,-xSize/4,1.25*xSize);
	//auto h2 = new TH1D("h2","y-coordinate",binnumber,-ySize/4,1.25*ySize);
	//auto h3 = new TH1D("h3","z-coordinate",binnumber,-tau/5,5*tau);

	//double esigma= pow(emeannumber,0.5); //raiz de emeannumber
	//auto h4 = new TH1D("h4","number of electrons",binnumber,emeannumber-5*esigma,emeannumber+5*esigma);
		auto *h5 = new TH1D("h5","sigma",200,0,40);
	//auto h6 = new TH1D("h6","spot x",binnumber,-xSize/2,3*xSize/2);
	//auto h7 = new TH1D("h7","spot y",binnumber,-ySize/2,3*ySize/2);


// CCD - Interactions Only
	TCanvas *ch2p2 = new TCanvas("ch2p2","ch2p2",1000,1000);

	TH2F *h2p_int = new TH2F("h2p_int","",nx,-0,xSize,ny,0,ySize); //este es el histograma 2D. xbins, xmin, xmax y luego y idem.
	h2p_int->SetName("Pixels");
	h2p_int->SetTitle("CCD Interactions");

// CCD - Dark Current Only
	TH2F *h2p_DC = new TH2F("h2p_DC","",nx,-0,xSize,ny,0,ySize); //este es el histograma 2D. xbins, xmin, xmax y luego y idem.
	h2p_DC->SetName("Pixels");
	h2p_DC->SetTitle("CCD Dark Current");

// CCD - Interactions + Dark Current
	TH2F *h2p_TOTAL = new TH2F("h2p_TOTAL","",nx,-0,xSize,ny,0,ySize); //este es el histograma 2D. xbins, xmin, xmax y luego y idem.
	h2p_TOTAL->SetName("Pixels");
	h2p_TOTAL->SetTitle("CCD Interactions + Dark Current");

cout<<"Histograms already booked"<<endl;

cout<<"Monte Carlo simulations runing ..."<<endl;

// This loop runs over each X-ray interaction //////////////////////////
//interaction(h1,h2,h3,h4,h5,h6,h7, h2p_int, h2p_TOTAL);
interaction(h2p_int, h2p_TOTAL,N0, A, B, h5);

cout<<endl;
cout<<"Interactions done"<<endl;

// Add dark Current ////////////////////////////////////////////////////
addDC(h2p_DC, h2p_TOTAL, darkC);
cout<<"Dark current simulation done"<<endl;

// Save the content of each histogram into pix variables ///////////////
pix(pix_int,pix_dc,pix_total, h2p_int, h2p_DC, h2p_TOTAL);
cout<<"Content of each histogram saved into pix variables"<<endl;

// Save into fits
save(nx, ny, sizeArray, pix_int, pix_dc, pix_total, N0, darkC, A ,B, R);
cout<<"Content of pix variables saved into fits files"<<endl;

// Histograms (not checked)
// hist(h1,h2,h3,h4,h5,h6,h7);

// Show CCD 2D plot ////////////////////////////////////////////////////
	ch2p2->Divide(2,2);
	ch2p2->cd(1);
	h2p_int->Draw("COLZ"); // Interactions
	ch2p2->cd(2);
	h2p_DC->Draw("COLZ"); // Dark Current
	ch2p2->cd(3);
	h2p_TOTAL->Draw("COLZ"); // Interactions + Dark Current
	ch2p2->cd(4);
	h5->Draw("COLZ"); // Interactions + Dark Current

// Print TCanvas into pdf
	TString ps = "./MC/CCD_N0"+std::to_string(N0)+"DC"+std::to_string(darkC)+"A"+ std::to_string(A)+"B"+ std::to_string(B)+"R"+ std::to_string(R)+".pdf";
	ch2p2->Print(ps+"[");
	ch2p2->Print(ps);
	ch2p2->Print(ps+"]");

	t.Stop();
	t.Print();

return 0;
}    // end
