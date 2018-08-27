#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>
#include <sys/resource.h>

#include "globalConstants.h"


#include "TFile.h"
#include "TNtuple.h"
#include "TObject.h"
#include "tinyxml2.h"
#include "gConfig.h"


using namespace std;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program extracts hits and tracks, computes their relevant parameters\n";
    cout << "and saves them in a root file.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << yellow;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> \n\n";
  cout << "\nOptions:\n";
  cout << "  -q for quiet (no screen output)\n";
  cout << "  -m <mask file> to provide a bad pixels mask\n";
  cout << "  -s <HDU number> for processing a single HDU \n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

string bitpix2TypeName(int bitpix){
  
  string typeName;
  switch(bitpix) {
      case BYTE_IMG:
          typeName = "BYTE(8 bits)";
          break;
      case SHORT_IMG:
          typeName = "SHORT(16 bits)";
          break;
      case LONG_IMG:
          typeName = "INT(32 bits)";
          break;
      case FLOAT_IMG:
          typeName = "FLOAT(32 bits)";
          break;
      case DOUBLE_IMG:
          typeName = "DOUBLE(64 bits)";
          break;
      default:
          typeName = "UNKNOWN";
  }
  return typeName;
}

struct track_t{
//   TNtuple &nt;
  vector<Int_t>    xPix;
  vector<Int_t>    yPix;
  vector<Float_t>  adc;
  
  Int_t flag;
  Int_t nSat;
  Int_t id;
  
  Int_t xMin;
  Int_t xMax;
  Int_t yMin;
  Int_t yMax;
  
  track_t() : flag(0), nSat(0), xMin(0), xMax(0), yMin(0), yMax(0) {};
  void fill(const Int_t &x , const Int_t &y, const Float_t &adcVal, const Int_t &l){ xPix.push_back(x); yPix.push_back(y); adc.push_back(adcVal); };
  
  void reset(){ xPix.clear(); yPix.clear(); adc.clear(); 
                flag=0; nSat=0; xMin=0; xMax=0; yMin=0; yMax=0; };
};

void extractTrack(double* outArray, const int &i, const int &nX, const int &nY, track_t &hit, const char* mask, const double &kCal, int recLevel = 0){
  
  if(recLevel>gHitMaxSize)
    return;
  
  int hitX = i%nX;
  int hitY = i/nX;
  const double &adcPixi = outArray[i];
  hit.flag = hit.flag|mask[i];

  if(adcPixi>kSat){
    hit.flag = hit.flag|kSatFlag;
    hit.nSat += 1;
  }
  hit.fill(hitX, hitY, adcPixi/kCal, 0);
  
  outArray[i] = kExtractedMask-hit.id;
  
  const double kAddThr = kCal/2;
  //West
  if(hitX>0){
    const double &pixAdcVal = outArray[i-1];
    if(pixAdcVal>kAddThr){
      extractTrack(outArray, i-1, nX, nY, hit, mask, kCal, recLevel+1);
    }
  }
  
  //South
  if(hitY>0){
    const double &pixAdcVal = outArray[i-nX];
    if(pixAdcVal>kAddThr){
      extractTrack(outArray, i-nX, nX, nY, hit, mask, kCal, recLevel+1);
    }
  }
  
  //East
  if(hitX<nX-1){
    const double &pixAdcVal = outArray[i+1];
    if(pixAdcVal>kAddThr){
      extractTrack(outArray, i+1, nX, nY, hit, mask, kCal, recLevel+1);
    }
  }
  
  //North
  if(hitY<nY-1){
    const double &pixAdcVal = outArray[i+nX];
    if(pixAdcVal>kAddThr){
      extractTrack(outArray, i+nX, nX, nY, hit, mask, kCal, recLevel+1);
    }
  }
  
  return;
}

void computeHitParameters(const track_t &hit, const double &kCal, Float_t *hitParam){
  
  double xSum   = 0;
  double ySum   = 0;
  double wSum   = 0;
  double eSum = 0;
  int    nPix   = 0;
  for(unsigned int i=0;i<hit.xPix.size();++i){
    const double &adc = hit.adc[i]; 
    if(adc>kExtractedMask){
      if(adc>0){
      	xSum += hit.xPix[i]*adc;
      	ySum += hit.yPix[i]*adc;
      	wSum += adc;
      }
      if(adc<kSat) eSum += std::round(adc);
      ++nPix;
    }
  }
  float xBary = xSum/wSum;
  float yBary = ySum/wSum;
  
  double x2Sum = 0;
  double y2Sum = 0;
  for(unsigned int i=0;i<hit.xPix.size();++i){
    const double &adc = hit.adc[i];
    if(adc>kExtractedMask){
      if(adc>0){
      	double dx = (hit.xPix[i] - xBary);
      	double dy = (hit.yPix[i] - yBary);
      	x2Sum   += dx*dx*adc;
      	y2Sum   += dy*dy*adc;
      }
    }
  }
  float xVar = x2Sum/wSum;
  float yVar = y2Sum/wSum;
  
  hitParam[0] = eSum;
  hitParam[1] = nPix;
  hitParam[2] = xBary;
  hitParam[3] = yBary;
  hitParam[4] = xVar;
  hitParam[5] = yVar;
}

struct hitTreeEntry_t{
  track_t hit;
  Int_t   runID;
  Int_t   ohdu;
  Int_t   expoStart;
  Int_t   nSavedPix;
  
  Float_t *hitParam;
  
  hitTreeEntry_t(const Int_t nPar): runID(-1), ohdu(-1), expoStart(-1),nSavedPix(0){ hitParam = new Float_t[nPar]; };
  ~hitTreeEntry_t(){ delete[] hitParam; };
};

void initHitTree(TTree &hitSumm, hitTreeEntry_t &evt ){
  hitSumm.Branch("runID",    &(evt.runID),    "runID/I");
  hitSumm.Branch("ohdu",     &(evt.ohdu),     "ohdu/I");
  hitSumm.Branch("expoStart",     &(evt.expoStart),     "expoStart/I");
  
  hitSumm.Branch("nSat", &(evt.hit.nSat), "nSat/I");
  hitSumm.Branch("flag", &(evt.hit.flag), "flag/I");
  hitSumm.Branch("xMin", &(evt.hit.xMin), "xMin/I");
  hitSumm.Branch("xMax", &(evt.hit.xMax), "xMax/I");
  hitSumm.Branch("yMin", &(evt.hit.yMin), "yMin/I");
  hitSumm.Branch("yMax", &(evt.hit.yMax), "yMax/I");
  
  gConfig &gc = gConfig::getInstance();
  const int kHitMaxSize = gc.getHitMaxSize();
  
  for(int n=0;n<gNExtraTNtupleVars;++n){
    hitSumm.Branch(gExtraTNtupleVars[n],  &(evt.hitParam[n]),  (string(gExtraTNtupleVars[n])+"/F").c_str());
  }
  
  hitSumm.Branch("nSavedPix", &(evt.nSavedPix), "nSavedPix/I");
  evt.hit.xPix.reserve(kHitMaxSize+1);
  hitSumm.Branch("xPix",  &(evt.hit.xPix[0]), "xPix[nSavedPix]/I");
  evt.hit.yPix.reserve(kHitMaxSize+1);
  hitSumm.Branch("yPix",  &(evt.hit.yPix[0]), "yPix[nSavedPix]/I");
  evt.hit.adc.reserve(kHitMaxSize+1);
  hitSumm.Branch("ePix", &(evt.hit.adc[0]), "ePix[nSavedPix]/F");
  
}

void refreshTreeAddresses(TTree &hitSumm, hitTreeEntry_t &evt)
{
  hitSumm.SetBranchAddress("xPix",  &(evt.hit.xPix[0]));
  hitSumm.SetBranchAddress("yPix",  &(evt.hit.yPix[0]));
  hitSumm.SetBranchAddress("ePix", &(evt.hit.adc[0]));
}


std::string trim(const std::string& str, const std::string& whitespace = " \t\'"){ // removes leading and trailing spaces
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos) return ""; // no content
    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

int getExpoInfoFromHeader(fitsfile *fptr, Int_t &tShut, Int_t &tExpo){
	int status = 0;   /*  CFITSIO status value MUST be initialized to zero!  */

	int hdutype;
	char keyValue[1024] = "";
	char comment[1024]  = "";

	Int_t tDate = -1;
	tShut = -1;
	tExpo = -1;
	// Get general data from ext=0 (hdu=1) header
	fits_movabs_hdu(fptr, 1, &hdutype, &status);
	if(status!=0) return -2;

	fits_read_keyword(fptr, "DATE", keyValue, comment, &status);
	if(status==0){ // key exist
		std::tm tm = {};
		strptime(trim(keyValue).c_str(), "%Y-%m-%dT%H:%M:%S", &tm);
		time_t t = std::mktime(&tm);
		tDate = t;	    
	}

	fits_read_keyword(fptr, "UTSHUT", keyValue, comment, &status);
		if(status==0){ // key exist
		std::tm tm = {};
		strptime(trim(keyValue).c_str(), "%Y-%m-%dT%H:%M:%S", &tm);
		time_t t = std::mktime(&tm);
		tShut = t;	    
	}
	tExpo = tDate - tShut;

	return 0;
}


int searchForTracks(TFile *outF, TTree &hitSumm, hitTreeEntry_t &evt, double* outArray, const int runID, const int ohdu, const int expoStart, const long totpix, const int nX, const int nY, char* mask){
  
  gConfig &gc = gConfig::getInstance();
  const double kCal       = gc.getExtCal(ohdu);
  const double kSeedThr   = kCal/2;
  const bool  kSaveTracks = gc.getSaveTracks();
  const char *kTrackCuts  = gc.getTracksCuts().c_str();
  
  const int kNVars = gNBaseTNtupleVars + gNExtraTNtupleVars;
  Float_t *ntVars = new Float_t[kNVars];
  if(mask == 0){
  	mask = new char[nX*nY]();
  } 
  
  outF->cd();
  
  unsigned int hitN = hitSumm.GetEntries();
  
  evt.runID = runID;
  evt.ohdu = ohdu;
  evt.expoStart = expoStart;
  track_t &hit = evt.hit;
  
  TTree hitSummAux("hitSummAux","hitSummAux");
  if(kSaveTracks){
    initHitTree(hitSummAux, evt);
    hitSummAux.SetCircular(1);
  }
  
  
  if(gVerbosity){
    cout << "\nProcessing runID " << runID << " ohdu " << ohdu << " -> sigma: " << gc.getExtSigma(ohdu) << ":\n";
  }
  for(long i=0;i<totpix;++i){
    
    if(outArray[i]>kSeedThr){
      
      evt.nSavedPix = 0;
      hit.reset();
      hit.id = hitN;
      ++hitN;
      extractTrack(outArray, i, nX, nY, hit, mask, kCal);
      computeHitParameters( hit, kCal, evt.hitParam );
      
      if(kSaveTracks){
        hitSummAux.Fill();
        if( hitSummAux.GetEntries(kTrackCuts) == 1 ){
          evt.nSavedPix = hit.xPix.size();
          refreshTreeAddresses(hitSumm, evt);
        }
      }
      
      hitSumm.Fill();
    }
    if(gVerbosity){
      if(i%1000 == 0) showProgress(i,totpix);
    }
  }
  delete[] ntVars;
  // delete[] mask;
  outF->cd();
  hitSumm.Write(hitSumm.GetName(),TObject::kOverwrite);
  if(gVerbosity){
    showProgress(1,1);
  }
  
  return 0;
}


bool readCardValue(fitsfile  *infptr, const char *keyName, double &value){
  
  int status = 0;
  char record[1024] = "";
  fits_read_card(infptr, keyName, record, &status);
  if(status==KEY_NO_EXIST){
    status=0;
    return false;
  }
  else{
    string sRec(record);
    size_t tPosEq = sRec.find("=");
    size_t tPosSl = sRec.find("/");
    istringstream recISS( sRec.substr(tPosEq+1, tPosSl-tPosEq-1) );
    recISS >> value;
    return true;
  }

}

int readMask(const char* maskName, vector <char*> &masks, const vector<int> &singleHdu){
  int status = 0;
  int nhdu = 0;
  double nulval = 0.;
  int anynul = 0;

  
  fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
  fits_open_file(&infptr, maskName, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  fits_get_num_hdus(infptr, &nhdu, &status);
  if (status != 0) return(status);
  
  /* check the extensions to process*/
  for(unsigned int i=0;i<singleHdu.size();++i){
    if(singleHdu[i] > nhdu){
      fits_close_file(infptr,  &status);
      cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      return -1000;
    }
  }
  
  vector<int> useHdu(singleHdu);
  if(singleHdu.size() == 0){
    for(int i=0;i<nhdu;++i){
      useHdu.push_back(i+1);
    }
  }
  const unsigned int nUseHdu=useHdu.size();
  
  
    
  for (unsigned int eN=0; eN<nUseHdu; ++eN)  /* Main loop through each extension */
  {
    const int n = useHdu[eN];
    
    /* get input image dimensions and total number of pixels in image */
    int hdutype, bitpix, naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    for (int i = 0; i < 9; ++i) naxes[i] = 1;
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
    long totpix = naxes[0] * naxes[1];
//       double bzero;
//       ffgky(infptr, TBYTE, "BZERO", &bzero, NULL, &status);
//       if (status){
// 	status = 0;
// 	bzero = 0.0;
//       }
    
    /* Don't try to process data if the hdu is empty */    
//       cout << (hdutype != IMAGE_HDU) << (naxis == 0) << (totpix == 0) << endl;
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      masks.push_back(0);
      continue;
    }
    
    char* maskArray = new char[totpix];
//       for(int i=0;i<totpix;++i) outArray[i] = 0;
    
    /* Open the input file */
    fits_movabs_hdu(infptr, n, &hdutype, &status);
    if (status != 0) return(status);
    int xMin=1;
    int xMax=naxes[0];
    int yMin=1;
    int yMax=naxes[1];
    
    /* Read the images as doubles, regardless of actual datatype. */
    long fpixel[2]={xMin,yMin};
    long lpixel[2]={xMax,yMax};
    long inc[2]={1,1};
    fits_read_subset(infptr, TBYTE, fpixel, lpixel, inc, &nulval, maskArray, &anynul, &status);
    if (status != 0){
      fits_report_error(stderr, status);
      return(status);
    }
    masks.push_back(maskArray);
  }
  fits_close_file(infptr,  &status);
  return status;
}

void writeConfigTree(TFile *outF){
  
  gConfig &gc = gConfig::getInstance();
  
  outF->cd();
  TTree configTree("config","config");
  
  Float_t kSigma[101];
  for(int i=0;i<=100;++i){
    kSigma[i] = gc.getExtSigma(i);
  }
  configTree.Branch("sigma", &kSigma, "sigma[101]/F");
  
  Float_t kCal[101];
  for(int i=0;i<=100;++i){
    kCal[i] = gc.getExtCal(i);
  }
  configTree.Branch("sigma", &kCal, "cal[101]/F");
  
  Float_t kSeedThr = gc.getSeedThr();
  configTree.Branch("seedThr", &kSeedThr, "seedThr/F");
  
  Float_t kAddThr = gc.getAddThr();
  configTree.Branch("addThr", &kAddThr, "addThr/F");
  
  Int_t kStackSize  = gc.getStackSize();
  configTree.Branch("stackSize", &kStackSize, "stackSize/I");
  
  Int_t kHitMaxSize  = gc.getHitMaxSize();
  configTree.Branch("hitMaxSize", &kHitMaxSize, "hitMaxSize/I");
 
  Bool_t kSaveTracks  = gc.getSaveTracks();
  configTree.Branch("saveTracks", &kSaveTracks, "saveTracks/B");
  
  TString kTrackCuts  = gc.getTracksCuts();
  configTree.Branch("trackCuts", &kTrackCuts);

  configTree.Fill();
  configTree.Write();
}

int computeImage(const vector<string> &inFileList,const char *maskName, const char *outFile, const vector<int> &singleHdu){
  int status = 0;
  double nulval = 0.;
  int anynul = 0;
  
  // Read mask
  bool noMask = false;
  vector <char*> masks;
  if(strlen(maskName) != 0) readMask(maskName, masks, singleHdu);
  else{
  	noMask = true;
  } 
  
  int nhdu = 0;
  const unsigned int nFiles  = inFileList.size();
  
 // gConfig &gc = gConfig::getInstance();
  
  TFile outRootFile(outFile, "RECREATE");
  
  writeConfigTree(&outRootFile);
  
  TTree hitSumm("hitSumm","hitSumm");
  hitTreeEntry_t evt(gNExtraTNtupleVars);
  initHitTree(hitSumm, evt);

  Int_t xSize = -1;
  Int_t ySize = -1;
  Int_t tShut = -1;
  Int_t tExpo = -1;
  
  for(unsigned int fn=0; fn < nFiles; ++fn){
    
    fitsfile  *infptr; /* FITS file pointers defined in fitsio.h */
    fits_open_file(&infptr, inFileList[fn].c_str(), READONLY, &status); /* Open the input file */
    if (status != 0) return(status);
    fits_get_num_hdus(infptr, &nhdu, &status);
    if (status != 0) return(status);

    getExpoInfoFromHeader(infptr, tShut, tExpo);
      
    /* check the extensions to process*/
    for(unsigned int i=0;i<singleHdu.size();++i){
      if(singleHdu[i] > nhdu){
      	fits_close_file(infptr,  &status);
      	cerr << red << "\nError: the file does not have the required HDU!\n\n" << normal;
      	return -1000;
      }
    }
    
    vector<int> useHdu(singleHdu);
    if(singleHdu.size() == 0){
      for(int i=0;i<nhdu;++i){
      	useHdu.push_back(i+1);
      }
    }
    const unsigned int nUseHdu=useHdu.size();
    
    
    for (unsigned int eN=0; eN<nUseHdu; ++eN)  /* Main loop through each extension */
    {
      
      const int n = useHdu[eN];

      /* get input image dimensions and total number of pixels in image */
      int hdutype, bitpix, naxis = 0;
      long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      fits_movabs_hdu(infptr, n, &hdutype, &status);
      for (int i = 0; i < 9; ++i) naxes[i] = 1;
      fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);
      long totpix = naxes[0] * naxes[1];
      double bzero;
      ffgky(infptr, TDOUBLE, "BZERO", &bzero, NULL, &status);
      if (status){
      	status = 0;
      	bzero = 0.0;
      }
      
      /* Don't try to process data if the hdu is empty */    
//       cout << (hdutype != IMAGE_HDU) << (naxis == 0) << (totpix == 0) << endl;
      if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      	continue;
      }
      cout << "\nImage size: " << naxes[0] << "x" << naxes[1] << endl;
      
      double* outArray = new double[totpix];
      for(int i=0;i<totpix;++i) outArray[i] = 0;
      
      
      /* Open the input file */
      fits_movabs_hdu(infptr, n, &hdutype, &status);
      if (status != 0) return(status);
      int xMin=1;
      int xMax=naxes[0];
      int yMin=1;
      int yMax=naxes[1];
      if(xSize<naxes[0]) xSize = naxes[0];
	  if(ySize<naxes[1]) ySize = naxes[1];
      
      /* Read the images as doubles, regardless of actual datatype. */
      long fpixel[2]={xMin,yMin};
      long lpixel[2]={xMax,yMax};
      long inc[2]={1,1};
      
      fits_read_subset(infptr, TDOUBLE, fpixel, lpixel, inc, &nulval, outArray, &anynul, &status);
      if (status != 0) return(status);
            
      double runID = 0;
      readCardValue(infptr, "RUNID", runID);
      double ext = n;
      readCardValue(infptr, "OHDU", ext);
      double expoStart = 0;
      readCardValue(infptr, "EXPSTART", expoStart);

      if(noMask) searchForTracks(&outRootFile, hitSumm, evt, outArray, (int)runID, (int)ext, (int)expoStart, totpix, naxes[0], naxes[1], 0);
      else       searchForTracks(&outRootFile, hitSumm, evt, outArray, (int)runID, (int)ext, (int)expoStart, totpix, naxes[0], naxes[1], masks[eN]);
      /* clean up */
      delete[] outArray;
    }

    /* Close the input file */
    fits_close_file(infptr,  &status);   
    
  }
  outRootFile.cd();
  TTree imgParTree("imgParTree", "imgParTree");
  imgParTree.Branch("xSize",    &xSize, "xSize/I");
  imgParTree.Branch("ySize",    &ySize, "ySize/I");
  imgParTree.Branch("tShut",    &tShut, "tShut/I");
  imgParTree.Branch("tExpo",    &tExpo, "tExpo/I");
  imgParTree.Fill();
  imgParTree.Write();
  outRootFile.Close();
  
  return status;
}


void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], 
                           vector<int> &singleHdu, vector<string> &inFileList, string &maskFile, string &outFile, string &confFile){
  
  if(argc == 1) return 1;
  inFileList.clear();
  singleHdu.clear();
  bool outFileFlag = false;
  bool maskFileFlag = false;
  int opt=0;
  while ( (opt = getopt(argc, argv, "i:m:o:s:c:qhH?")) != -1) {
    switch (opt) {
      case 'm':
        if(!maskFileFlag){
          maskFile = optarg;
          maskFileFlag = true;
        }
        else{
          cerr << red << "\nError, can not set more than one mask file!\n\n" << normal;
          return 2;
        }
      break;
      case 'o':
        if(!outFileFlag){
          outFile = optarg;
          outFileFlag = true;
        }
        else{
          cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
          return 2;
        }
        break;
      case 's':
        singleHdu.push_back(atoi(optarg));
        break;
      case 'c':
        confFile = optarg;
        break;
      case 'q':
        gVerbosity = 0;
        break;
      case 'h':
      case 'H':
      default: /* '?' */
        return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }
  
  if(!maskFileFlag){
    cerr << yellow << "\nMask filename missing. Will use empty mask\n" << normal;
    maskFile = "";
  }

  for(int i=optind; i<argc; ++i){
    inFileList.push_back(argv[i]);
    if(!fileExist(argv[i])){
      cout << red << "\nError reading input file: " << argv[i] <<"\nThe file doesn't exist!\n\n" << normal;
      return 1;
    }
  }
  
  if(inFileList.size()==0){
    cerr << red << "Error: no input file(s) provided!\n\n" << normal;
    return 1;
  }
  
  return 0;
}

#include <iostream>
#include <sstream>
#include <locale>
#include <iomanip>
int main(int argc, char *argv[])
{
  checkArch(); //Check the size of the double and float variables.
  
  time_t start,end;
  double dif;
  time (&start);
  
  string maskFile;
  string outFile;
  string confFile = "extractConfig.xml";
  vector<string> inFileList;
  vector<int> singleHdu;

  int returnCode = processCommandLineArgs( argc, argv, singleHdu, inFileList, maskFile, outFile, confFile);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }

  /* Create configuration singleton and read configuration file */
  gConfig &gc = gConfig::getInstance();
  if(gc.readConfFile(confFile.c_str()) == false){
    return 1;
  }
  if(gVerbosity){
    cout << "\nConfig file: " << confFile << endl;
    gc.printVariables();
  }
  gHitMaxSize = gc.getHitMaxSize();

  
  /* Increase the stack size to be able to use a deeply
   * nested recursive function */
  const rlim_t kStackSize = gc.getStackSize() * 1024 * 1024;   // min stack size = gc.getStackSize() MB
  struct rlimit rl;
  int result;
  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0){
      if (rl.rlim_cur < kStackSize){
          rl.rlim_cur = kStackSize;
          result = setrlimit(RLIMIT_STACK, &rl);
          if (result != 0){
	    cerr << "Error increasing the stack size: setrlimit returned result = " << result << endl;
          }
      }
  }
  
  
  
  if(gVerbosity){
    cout << bold << "\nWill read the following files:\n" << normal;
    for(unsigned int i=0; i<inFileList.size();++i) cout << "\t" << inFileList[i] << endl;
    if(singleHdu.size()>0){
      cout << bold << "\nAnd the following extension:" << normal << endl << "\t";
      for(unsigned int i=0; i<singleHdu.size();++i) cout << singleHdu[i] << ",";
      cout << "\b " << endl;
    }
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }
  
  int status = computeImage( inFileList, maskFile.c_str(),  outFile.c_str(), singleHdu);
  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "\nAll done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  return status;
}
