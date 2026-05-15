#ifndef POSIFITTER_HH
#define POSIFITTER_HH

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector> 
#include <algorithm>
#include <map>
#include <unordered_map>
#include <regex>
#include <filesystem>
#include <chrono>


namespace fs = std::filesystem;

#include "TH1F.h" 
#include "TH2D.h" 
#include "TCanvas.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TGraph.h"

double fPDFWrapper(double* x, double* par);


struct Params{
    int istat;
    double t0, dt, aIoni, aoPs, constant;
    double t0_err[2], dt_err[2], aIoni_err[2], aoPs_err[2], constant_err[2];
    double chi2, chi2ndf, edm, chi2prescan;
    double corr[4*4]; //correlation matrix
    TH1F *hIoni, *hoPs, *hFit;
};

class PositroniumFitter
{
private:
    bool debug, setparams, hama;
    bool prescan;
    int xmin, xmax, nbin;
    double Coeff;
    std::string SimEnergy;
    std::string datatype;

    int loop = 0;
    double chi2_prescan;
    Params pInit;

    TMinuit* gMinuit;

    TH1F *hdata, *hPDF1, *hPDF2, *hIoni, *hoPs, *hFit;
    TH1F *hPrompt, *hDelay, *hDeltaT;
    TH1F *hTrueDeltaT;

    TGraph *gr;
    
private:

    void shiftHistogram(TH1F*, TH1F*, double, double);
    void LoadPDF();
    void SetupMinuit();
    double PreScanT0Dt(double t0min, double t0max, double dtmin, double dtmax, int nt0, int ndt);

    
public:
    
    PositroniumFitter(TH1F* h = nullptr, bool Debug = false, int firstbin = 0, int lastbin = 700, int nbins = 350, std::string type = "data");
    ~PositroniumFitter();
    
    static void FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
    
    void setType(std::string type){datatype = type;};
    void setHamaOnly(){hama = true;};
    void setEnergy(std::string Energy){SimEnergy = Energy;};
    void setParameters(Params);
    void setDataHist(TH1F* h){hdata = h;};
    // void SetPDF(std::string);
    void DoFit();
    double fPDF(double *x, double *par);
    Params FetchParams();

    TGraph* getContour(){return gr;};
    
};



#endif
