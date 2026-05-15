#ifndef FITHH
#define FITHH

#include <iostream> 
#include <vector>
#include <algorithm> 
#include <map>
#include <filesystem>
#include <regex>

#include "PositroniumFitter.h"

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TMinuit.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TF1.h"
#include "TString.h"
#include "TRatioPlot.h"
#include "TMarker.h"

namespace fs = std::filesystem;

// Function 
void LoadData();
void LoadMC(std::string);
void LoadTree();
void LoadPMTInfo();
void Plot(Params, std::string, int evt = 0);

// Variables
bool SaveGamma, HamaOnly;
bool SaveElectron;
bool uniformMC;
std::string data_type;
std::string juno_version;
std::string SimEnergy;
std::string suffix;
std::string reprod;


int start;
int nEntries;
int nbpIBD;
int nbdIBD;

int xmin, xmax, nbin;
int maxNFit;

//Fit parameters
double A1, A2, t0, dt, bkg;

TFile* outfile;

TH1F* hdata;
TH1F* hIoni, *hAni, *hFit;

TH1F* hPrompt;
TH1F* hDelay; 
TH1F* hDeltaT;
TH1F* hoPsDeltaT;
TH1F* hTrueoPsDeltaT;


TTree* FitTree;
int Entry; 
int istat; 
float edm;
float charge;
float energy;
float dataRMS;
float chi2;
float chi2_ndf;
float chi2_prescan;
float tfirst, dtime, aIoni, aoPs, constant;
float tfirst_err_low, dtime_err_low, aIoni_err_low, aoPs_err_low, constant_err_low;
float tfirst_err_up, dtime_err_up, aIoni_err_up, aoPs_err_up, constant_err_up;
double corr[4*4];
float Simx, Simy, Simz;

TTree* PairTree;
int TotdEntries;

TFile* simoutfile;
TFile* gammaoutfile;
TTree* SimFit;
TTree* SimGammaFit;
int SimEntry = 0;
double Simdt = 0.0;
double Fitdt = 0.0;
double RMS = 0.0;
double std_dev = 0.0;
double mean = 0.0;
double skew = 0.0;
double kurto = 0.0;
double SimChi2Ndf = 0.0;
TString SimTag;

TChain IBD("ibd");
TChain newIBD("event");
float PromptE;
float DelayE;
float Prompt_vtx[3];
TString* tag = nullptr;
std::vector<double>* PromptTOF = nullptr;
std::vector<double>* DelayTOF = nullptr;

TChain Sim("Sim");
int SimEvtID;
std::vector<int>* PDGID;
std::vector<double>* InitT;
std::vector<float>* SimX;
std::vector<float>* SimY;
std::vector<float>* SimZ;
double SimDecayTime;

TChain Calib("Calib");
int CalibEvtID;
double ChargeTotLPMT;
std::vector<int>* PmtIDCalib;
std::vector<double>* CalibTimeTOF;

TChain Reco("Reco");
float Recx;
float Recy;
float Recz;

TChain GammaCalib("Calib");
int GammaCalibEvtID;
std::vector<int>* GammaPmtIDCalib;
std::vector<double>* GammaCalibTimeTOF;


TChain PMT("PMTInfo");
std::map<int, std::string> PmtType;

std::map<int, double> TruedT;


#endif 

