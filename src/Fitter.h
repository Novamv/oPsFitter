#pragma once

#include <string>
#include <array>
#include <memory>
#include <chrono>
#include <iostream>

// ROOT
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"

// Minuit2 — modern, thread-safe, no global state
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

// ─────────────────────────────────────────────────────────────────────────────
// Params  —  plain data struct returned by FetchParams()
// ─────────────────────────────────────────────────────────────────────────────
struct FitParameters {
    double t0    = 0, dt   = 0;
    double aIoni = 0, aoPs = 0, constant = 0;

    // Asymmetric errors: [0] = lower (negative), [1] = upper (positive)
    std::array<double,2> t0_err       = {0,0};
    std::array<double,2> dt_err       = {0,0};
    std::array<double,2> aIoni_err    = {0,0};
    std::array<double,2> aoPs_err     = {0,0};
    std::array<double,2> constant_err = {0,0};

    double chi2      = 0, chi2ndf   = 0, chi2prescan = 0;
    double edm       = 0;
    int    istat     = 0;   // Minuit2: 0=not conv, 1=cov forced pos-def, 2=hesse ok, 3=minos ok

    // Correlation matrix — 5×5 flattened (parameters: t0, dt, aoPs, const, aIoni)
    std::array<double,25> corr = {};

    // Histograms for plotting (owned by the fitter, valid until the next DoFit call)
    TH1F* hIoni = nullptr;
    TH1F* hoPs  = nullptr;
    TH1F* hFit  = nullptr;
};


class Fitter
{
public:
    // ── Parameters ─────────────────────────────────────────────────────────
    static constexpr int kT0    = 0;
    static constexpr int kDt    = 1;
    static constexpr int kAIoni = 2;
    static constexpr int kAoPs  = 3;
    static constexpr int kConst = 4;
    static constexpr int kNPar  = 5;
private:

    // ── Init ────────────────────────────────────────────────────────
    TH1F* hdata = nullptr;

    // Normalised templates (unit integral), owned by this object
    std::unique_ptr<TH1F> hPDF_ioni;   // ionisation template (Ge68 or gamma)
    std::unique_ptr<TH1F> hPDF_oPs;    // o-Ps template (Ge68 shifted)

    // Working histograms updated by shiftHistogram inside the NLL
    std::unique_ptr<TH1F> hIoni_work;
    std::unique_ptr<TH1F> hoPs_work;
    std::unique_ptr<TH1F> hFit_work;

    bool        debug, hamaOnly;
    double      xmin, xmax;
    int         nbin;
    std::string datatype;

    bool converged;

    std::string pdf_source = "Ge68";

    // oPs amplitude conversion factor
    double oPsCoeff     = 0.0; // set in LoadPDF
    double oPsCoeff_err = 0.0; // set in LoadPDF

    // Seed from previous fit (optional)
    bool hasSeed = false;
    FitParameters seed;

    // ── Fit infrastructure ──────────────────────────────────────────────────
    std::unique_ptr<ROOT::Minuit2::Minuit2Minimizer> migrad;
    
    // The objective function: extended Poisson log-likelihood
    // par = {t0, dt, aIoni, aoPs, constant}
    double NLL(const double* par);
    
    // ── 2D pre-scan ────────────────────────────────────────────────────────
    bool run_prescan = true;
    double chi2prescan = 0.0;
    struct ScanResult { double t0, dt, chi2; };
    ScanResult PreScanT0Dt(double* init_values, int nt0 = 20, int ndt = 20);

    // ── PDF manipulation ───────────────────────────────────────────────
    // Does h_out = intensity * PDF(t - shift)
    void shiftHistogram(const TH1F* src, TH1F* dst, double shift, double intensity);

    // ── PDF loading ─────────────────────────────────────────────────────────
    void LoadPDF();

    // ── PDF loading ─────────────────────────────────────────────────────────
    double GetInflectionPoint(TH1F* h);

    // ── Minuit2 parameter config ────────────────────────────────────────────
    struct ParCfg {
        const char* name;
        double      step;
        double      lo;    // relative lower bound factor (e.g. -0.6 means init-60%)
        double      hi;    // relative upper bound factor
    };
    static const ParCfg kParConfig[kNPar];

public:
    // ── Constructor / destructor ────────────────────────────────────────────
    Fitter(bool debug, bool hama, double min, double max, int bins, std::string type);
    ~Fitter();

    void setDataHist(TH1F* h) {hdata = h;};
    void setHamaOnly() {hamaOnly = true;}
    
    // Run one fit routine
    void DoFit();

    // For Ionisation PDF
    void UsePDFSource(const std::string& source) {pdf_source = source;};
    
    // Extract results after DoFit()
    FitParameters FetchParams();
};
