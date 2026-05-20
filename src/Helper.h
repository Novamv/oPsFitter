#ifndef HELPER_HH
#define HELPER_HH

#include <iostream>
#include <map>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"


class Helper
{
private:
    // ── Init ────────────────────────────────────────────────────────
    bool debug;
    bool hamaOnly;
    int xmin, xmax, nbin;
    std::string datatype;

    std::string energy = "3";
    std::string juno_version = "J24.1.2";

    // ── PMT Information ───────────────────────────────────────────────────
    void LoadPMTInfo();
    
    // ── Create PDF ────────────────────────────────────────────────────────
    std::string suffix;
    bool CreateGe68PDF(const std::string& fname);

    double shift(TH1F* src, TH1F* dst, double TargetTime = 210 /*ns*/);

    // ── Sim ────────────────────────────────────────────────────────
    TChain Sim;
    int SimEvtID;
    std::vector<int>* PDGID = nullptr;
    std::vector<double>* InitT = nullptr;
      
public:
    Helper(bool debug_, int min_, int max_, int nbin_, std::string datatype_);
    ~Helper();
    
    // ── Utilities ────────────────────────────────────────────────────────
    std::map<int, std::string> PmtType;
    std::map<int, double> oPsTrueDecay;
    TH1F* hGe68 = nullptr;

    void ProgressBar(Long64_t current, Long64_t total, int width = 50);

    // ── PDF ────────────────────────────────────────────────────────
    void SetHamaOnly() {hamaOnly = true;};
    bool GeneratePDF();

    // ── Data ────────────────────────────────────────────────────────
    TChain Calib; // MC
    TChain Reco;  // MC
    TChain IBD;   // data

    int nEntries, EntryNb;
    TString Tag;
    float RecEnergy, Recx, Recy, Recz;
    std::vector<int>* PmtIDCalib = nullptr;
    std::vector<double>* HitTimeTOF = nullptr;

    void LoadData(bool delay = false);
};


#endif