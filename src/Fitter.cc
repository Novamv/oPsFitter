#include "Fitter.h"


// ─────────────────────────────────────────────────────────────────────────────
// Constructor / Destructor
// ─────────────────────────────────────────────────────────────────────────────
Fitter::Fitter(bool debug, bool hama, double min, double max, int bins, std::string type)
: debug(debug), hamaOnly(hama), xmin(min), xmax(max), nbin(bins), datatype(type)
{
    // ── Load PDF templates ──────────────────────────────────────────────────
    auto t0 = std::chrono::high_resolution_clock::now();
    LoadPDF();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "PositroniumFitter: PDF loaded in "
    << std::chrono::duration<double>(t1-t0).count() << " s\n";
    
    // ── Allocate working histograms once (reused every event) ───────────────
    hIoni_work = std::make_unique<TH1F>("_hIoni", "", nbin, xmin, xmax);
    hoPs_work  = std::make_unique<TH1F>("_hoPs",  "", nbin, xmin, xmax);
    hFit_work  = std::make_unique<TH1F>("_hFit",  "", nbin, xmin, xmax);
    hIoni_work->SetDirectory(nullptr);
    hoPs_work ->SetDirectory(nullptr);
    hFit_work ->SetDirectory(nullptr);

    // ── Set up Minuit2 minimizer ────────────────────────────────────────────
    migrad = std::make_unique<ROOT::Minuit2::Minuit2Minimizer>(ROOT::Minuit2::kMigrad);

}

Fitter::~Fitter()
{
}

// ─────────────────────────────────────────────────────────────────────────────
// Parameter configuration table
// ─────────────────────────────────────────────────────────────────────────────
const Fitter::ParCfg Fitter::kParConfig[kNPar] = {
    // name     step    lo_factor  hi_factor
    { "t0",    0.005,  -0.10,      0.10  },   // kT0    — (bounds in %)
    { "dt",    0.005,   0.00,      30.0  },   // kDt    — hard coded bounds
    { "aIoni", 100.0,  -0.30,      0.30  },   // kAIoni — (bounds in %)
    { "aoPs",  100.0,  -0.30,      0.30  },   // kAoPs  — (bounds in %)
    { "const", 0.1,    -0.10,      0.10  },   // kConst — (bounds in %)
};


// ─────────────────────────────────────────────────────────────────────────────
// DoFit — run one complete fit on hdata
// ─────────────────────────────────────────────────────────────────────────────
double Fitter::GetInflectionPoint(TH1F* h){
    double peak_x   = h->GetBinCenter(h->GetMaximumBin());
    double peak_val = h->GetMaximum();
    double fit_lo = peak_x - 50.0;
    double fit_hi = peak_x + 5.0;

    // f(x) = A / (1 + exp(-k*(x - x0))) sigmoid
    // par[0] = A    (amplitude ~ peak value)
    // par[1] = k    (steepness of rise)
    // par[2] = x0   (inflection point — our alignment reference)
    TF1 fsig("_fsig", "[0]/(1+exp(-[1]*(x-[2])))", fit_lo, fit_hi);
    fsig.SetParameter(0, peak_val);
    fsig.SetParameter(1, 1.0);      // positive: rising left to right
    fsig.SetParameter(2, peak_x - 5.0);
    fsig.SetParLimits(1, 0.1, 10.0);  // enforce rising slope

    h->Fit(&fsig, "RQ0");

    return fsig.GetParameter(2);

}
void Fitter::DoFit()
{
    if(!hdata) throw std::runtime_error("Fitter::DoFit(): hdata is null");

    double tStart, timediff, Amp_Ioni, Amp_oPs, bgnd;

    if(hasSeed){
        tStart    = seed.t0;
        timediff  = seed.dt;
        Amp_Ioni  = seed.aIoni;
        Amp_oPs   = seed.aoPs;
        bgnd      = seed.constant;
        hasSeed   = false;
    } else {
        // oPs annihilation delay time
        timediff = 1.0;
        
        // background
        TF1 f_bkg("f_bkg", "[0]", 0, 150);
        hdata->Fit(&f_bkg, "RQ0");
        bgnd = f_bkg.GetParameter(0);
        
        // oPs amplitude: use oPsCoeff as first guess
        Amp_oPs = 1.022 * oPsCoeff;
    
        // Ionisation amplitude: whatever is left
        Amp_Ioni = hdata->Integral() -  Amp_oPs;

        // PDF initial shift
        double inflection_data = GetInflectionPoint(hdata);
        // for(int i = 1; i <= hPDF_ioni->GetNbinsX(); ++i){
        //     hIoni_work->SetBinContent(i, Amp_Ioni*hPDF_ioni->GetBinContent(i) + bgnd);
        // }        
        double inflection_ioni  = GetInflectionPoint(hPDF_ioni.get());
        tStart = inflection_data - inflection_ioni;
        // tStart = hdata->GetBinCenter(hdata->GetMaximumBin()) - hPDF_ioni->GetBinCenter(hPDF_ioni->GetMaximumBin());

        std::cout << "\n[DoFit] Inflection data: " << inflection_data << ", inflection Ioni: " << inflection_ioni << std::endl;
        
    }

    std::cout << "[DoFit] starting values: "
              << "  t0=" << tStart << "  dt=" << timediff
              << "  aIoni=" << Amp_Ioni << "  aoPs=" << Amp_oPs
              << "  bgnd=" << bgnd << std::endl;

    // ── Build initial parameter set ─────────────────────────────────────────
    double initVals[kNPar] = { tStart, timediff, Amp_Ioni, Amp_oPs, bgnd };
    
    // ── 2D pre-scan of (t0, dt) ─────────────────────────────────────────────
    if(debug){
        TCanvas* c = new TCanvas();
        shiftHistogram(hPDF_ioni.get(), hIoni_work.get(), initVals[kT0], initVals[kAIoni]);
        shiftHistogram(hPDF_oPs.get(), hoPs_work.get(), initVals[kT0]+initVals[kDt], initVals[kAoPs]);

        hdata->SetLineColor(kBlack);
        hdata->Draw();
        hIoni_work->SetLineColor(kRed+1);
        hIoni_work->Draw("HIST+SAME");
        hoPs_work->SetLineColor(kViolet+1);
        hoPs_work->Draw("HIST+SAME");

        c->SaveAs("PreScanFit.root");
    }

    // if(run_prescan){
    //     auto scan = PreScanT0Dt(initVals); // t0 and dt first guesses are modified here
    //     chi2prescan = scan.chi2;
    //     initVals[kT0] = scan.t0; initVals[kDt] = scan.dt;
    // }

    if(debug){
        TCanvas* c = new TCanvas();
        shiftHistogram(hPDF_ioni.get(), hIoni_work.get(), initVals[kT0], initVals[kAIoni]);
        shiftHistogram(hPDF_oPs.get(), hoPs_work.get(), initVals[kT0]+initVals[kDt], initVals[kAoPs]);

        hdata->SetLineColor(kBlack);
        hdata->Draw();
        hIoni_work->SetLineColor(kRed+1);
        hIoni_work->Draw("HIST+SAME");
        hoPs_work->SetLineColor(kViolet+1);
        hoPs_work->Draw("HIST+SAME");

        c->SaveAs("PostScanFit.root");
    }


    // Functor for minimizers
    ROOT::Math::Functor functor(this, &Fitter::NLL, 5);

    // ── Step 1: SIMPLEX (rough first search) ─────────────────────────────────────────────
    std::cout << "Info SIMPLEX \n";
    ROOT::Minuit2::Minuit2Minimizer simplex(ROOT::Minuit2::kSimplex);
    simplex.SetFunction(functor);
    simplex.SetStrategy(1);
    simplex.SetPrecision(1e-12);
    simplex.SetTolerance(0.1);
    simplex.SetMaxFunctionCalls(10000);
    simplex.SetPrintLevel(debug ? 2 : -1);


    for(int i = 0; i < kNPar; ++i){
        const auto& config = kParConfig[i];
        double lo, hi;
        if(i == kDt){                          // absolute bounds
            lo = config.lo;
            hi = config.hi;
        } else {                               // relative bounds
            lo = initVals[i] * (1.0 + config.lo);
            hi = initVals[i] * (1.0 + config.hi);
        }
        if(lo > hi) std::swap(lo, hi);        // guard for negative initVals

        // if(i == kConst) simplex.SetFixedVariable(i, config.name, initVals[i]);
        else simplex.SetLimitedVariable(i, config.name, initVals[i], config.step, lo, hi);
    }
    simplex.FixVariable(kConst);

    simplex.Minimize();
    const double* x_simplex = simplex.X();

    // ── Step 2: MIGRAD seeded from simplex (scan around minimum) ─────────────────────────────────────────────

    std::cout << "Info MIGRAD \n";

    const int N_ATTEMPTS = 5;
    const double tolerance_list[N_ATTEMPTS] = {0.01, 0.0025, 1e-6, 0.5, 1.0};
    const int maxFcn_list[N_ATTEMPTS] = {200000, 400000, 600000, 10000000, 20000};

    bool check = false;
    double best_chi2 = 1e99;
    double best_edm = 1e99;
    double best_par[kNPar] = {x_simplex[kT0], x_simplex[kDt], x_simplex[kAIoni], x_simplex[kAoPs], x_simplex[kConst]};

    migrad->SetFunction(functor);
    migrad->SetStrategy(2);
    migrad->SetErrorDef(1.0); // 1.0 for chi2, 0.5 for NLL
    migrad->SetPrintLevel(debug ? 2 : -1);

    
    for(int iter = 0; iter < N_ATTEMPTS && !check; ++iter){
        migrad->SetTolerance(tolerance_list[iter]);
        migrad->SetMaxFunctionCalls(maxFcn_list[iter]);
        migrad->SetMaxIterations(maxFcn_list[iter]);

        for(int i = 0; i < kNPar; ++i){
            const auto& config = kParConfig[i];
            double lo, hi;
            if(i == kDt){                          // absolute bounds
                lo = config.lo;
                hi = config.hi;
            } else {                               // relative bounds
                lo = best_par[i] * (1.0 + config.lo);
                hi = best_par[i] * (1.0 + config.hi);
            }
            if(lo > hi) std::swap(lo, hi);
    
            // if(i == kConst) migrad->SetFixedVariable(i, config.name, best_par[i]);
            else migrad->SetLimitedVariable(i, config.name, best_par[i], config.step, lo, hi);
        }
        
        converged = migrad->Minimize();
        std::cout << "convergence " << converged << std::endl;

        int status  = migrad->Status();
        double chi2 = migrad->MinValue();
        double edm  = migrad->Edm();

        const double* par = migrad->X();
        for(int i = 0; i < kNPar; ++i) best_par[i] = par[i];

        if (edm < best_edm){
            best_chi2 = chi2;
            best_edm = edm;
        }

        check = (status <= 1 && converged);
        std::cout << "Validity check = " << check << " for tolerance " << tolerance_list[iter] << "\n";
        if (check) break;
        
    }

    if(!converged && debug){
        std::cout << "[DoFit] MIGRAD did not converge sucessfully (status " << migrad->Status() << ")" << std::endl;
    }

    std::cout << "  [DoFit] istat=" << migrad->Status()
              << "  EDM=" << migrad->Edm()
              << "  chi2=" << best_chi2 << "\n";
}



// ─────────────────────────────────────────────────────────────────────────────
// FetchParams — extract results after DoFit()
// ─────────────────────────────────────────────────────────────────────────────
FitParameters Fitter::FetchParams()
{
    FitParameters out;

    if(migrad->Status() >= 2 && converged) migrad->Hesse();

    const double* x = migrad->X();
    const double* e = migrad->Errors();

    out.t0       = x[kT0];
    out.dt       = x[kDt];
    out.aIoni    = x[kAIoni];
    out.aoPs     = x[kAoPs];
    out.constant = x[kConst];

    // Symmetric errors from HESSE
    out.t0_err      [0] = out.t0_err      [1] = e[kT0];
    out.dt_err      [0] = out.dt_err      [1] = e[kDt];
    out.aIoni_err   [0] = out.aIoni_err   [1] = e[kAIoni];
    out.aoPs_err    [0] = out.aoPs_err    [1] = e[kAoPs];
    out.constant_err[0] = out.constant_err[1] = e[kConst];

    // Asymmetric MINOS errors — only attempted when HESSE suceeded
    // if(migrad->Status() == 3 || converged){
    //     double elo, ehi;
    //     for (int i = 0; i < kNPar; ++i){
    //         if(migrad->GetMinosError(i, elo, ehi)){
    //             auto setErr = [&](std::array<double,2>& arr){
    //                 arr[0] = (elo != 0) ? elo : -e[i];  // MINOS may return 0 on failure
    //                 arr[1] = (ehi != 0) ? ehi :  e[i];
    //             };
    //             if (i == kT0)    setErr(out.t0_err);
    //             if (i == kDt)    setErr(out.dt_err);
    //             if (i == kAIoni) setErr(out.aIoni_err);
    //             if (i == kAoPs)  setErr(out.aoPs_err);
    //             if (i == kConst) setErr(out.constant_err);
    //         }
    //     }
    // }

    
    // Quality-of-fit
    out.chi2        = migrad->MinValue();
    out.chi2ndf     = migrad->MinValue() / double(nbin - migrad->NFree());
    out.chi2prescan = chi2prescan;
    out.edm         = migrad->Edm();
    out.istat       = migrad->Status();

    // Correlation matrix (5×5 flattened)
    // Minuit2 exposes the full covariance via CovMatrix(i,j)
    if (migrad->Status() >= 2) {
        for (int i = 0; i < kNPar; ++i) {
            for (int j = 0; j < kNPar; ++j) {
                double cii = migrad->CovMatrix(i, i);
                double cjj = migrad->CovMatrix(j, j);
                double denom = std::sqrt(cii * cjj);
                out.corr[i * kNPar + j] = (denom > 0.0)
                    ? migrad->CovMatrix(i, j) / denom : 0.0;
            }
        }
    }

    // Histograms for plotting — fill from best-fit parameters
    shiftHistogram(hPDF_ioni.get(), hIoni_work.get(), out.t0,           out.aIoni);
    shiftHistogram(hPDF_oPs.get(),  hoPs_work.get(),  out.t0 + out.dt,  out.aoPs);
    hFit_work->Reset();
    hFit_work->Add(hIoni_work.get());
    hFit_work->Add(hoPs_work.get());
    for (int i = 1; i <= nbin; ++i)
        hFit_work->AddBinContent(i, out.constant);

    out.hIoni = hIoni_work.get();
    out.hoPs  = hoPs_work.get();
    out.hFit  = hFit_work.get();

    std::cout << "FetchParams: t0=" << out.t0 << " dt=" << out.dt
              << " aIoni=" << out.aIoni << " aoPs=" << out.aoPs
              << " const=" << out.constant
              << " istat=" << out.istat << " chi2/ndf=" << out.chi2ndf << std::endl;

    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// NLL — extended Poisson log-likelihood 
//
// par = { t0, dt, aIoni, aoPs, constant }
//
// L = 2 * sum_i [ mu_i - n_i + n_i * ln(n_i / mu_i) ]   (n_i > 0)
//   = 2 * sum_i [ mu_i ]                                (n_i = 0)
//
// plus a Gaussian penalty on aoPs to constrain it near the energy 
// expectation 1.022 * Coeff.
// ─────────────────────────────────────────────────────────────────────────────
double Fitter::NLL(const double* par)
{
    double t0    = par[kT0];
    double dt    = par[kDt];
    double aIoni = par[kAIoni];
    double aoPs  = par[kAoPs];
    double c     = par[kConst];

    double LogL = 0;

    shiftHistogram(hPDF_ioni.get(), hIoni_work.get(), t0,      aIoni);
    shiftHistogram(hPDF_oPs.get(),  hoPs_work.get(),  t0 + dt, aoPs);

    for (int i = 1; i <= nbin; ++i){
        const double n  = hdata->GetBinContent(i);
        const double mu = hIoni_work->GetBinContent(i) + hoPs_work->GetBinContent(i) + c;

        const double mu_safe = (mu > 0.0) ? mu : 1e-9;

        if(n > 0.0)
            LogL += 2.0 * (mu_safe - n + n*std::log(n / mu_safe));
        else
            LogL += 2.0 * mu_safe;
    }

    // Penalty: aoPs expected near 1.022 * Coeff, 40% tolerence (estimated using distribution of max)
    double expected_oPs = 1.022 * oPsCoeff; // Convert energy to amplitude
    double pull = (aoPs - expected_oPs) / (1.022 * oPsCoeff_err); // pearson chi2

    // double pull_ioni = (hdata->Integral() - aIoni - aoPs) / (oPsCoeff_err);
    // LogL += pull_ioni*pull_ioni;


    LogL += pull*pull;
    return LogL;
}
// ─────────────────────────────────────────────────────────────────────────────
// PreScanT0Dt — brute-force grid scan over (t0, dt)
// Returns the best (t0, dt, chi2) triplet.
// ─────────────────────────────────────────────────────────────────────────────
Fitter::ScanResult Fitter::PreScanT0Dt(double* init_values, int nt0, int ndt)
{
    const auto& conf_t0 = kParConfig[kT0];
    double t0min        = init_values[kT0] + init_values[kT0]*conf_t0.lo;
    double t0max        = init_values[kT0] + init_values[kT0]*conf_t0.hi;

    if(t0min > t0max) std::swap(t0min, t0max);

    const auto& conf_dt = kParConfig[kDt];
    double dtmin        = conf_dt.lo;
    double dtmax        = conf_dt.hi;

    std::cout << "[PreScan] grid " << nt0 << "×" << ndt
              << "  t0=[" << t0min << "," << t0max << "]"
              << "  dt=[" << dtmin << "," << dtmax << "]" << std::endl;
    
    ScanResult best { t0min, dtmin, 1e9 };

    for(int i = 0; i < nt0; ++i){
        init_values[kT0] = t0min + (t0max - t0min) * i / std::max(nt0 - 1, 1);
        for (int j = 0; j < ndt; ++j) {
            init_values[kDt] = dtmin + (dtmax - dtmin) * j / std::max(ndt - 1, 1);
            double chi2 = NLL(init_values);
            if (chi2 < best.chi2) { best = { init_values[kT0], init_values[kDt], chi2 }; }
        }
    }

    std::cout << "[PreScan] best  t0=" << best.t0
              << "  dt=" << best.dt << "  chi2=" << best.chi2 << std::endl;
    return best;
}
// ─────────────────────────────────────────────────────────────────────────────
// shiftHistogram — fill dst = intensity * src(t - shift), bin-by-bin
// dst must have the same axis as src (guaranteed by construction above).
// ─────────────────────────────────────────────────────────────────────────────
void Fitter::shiftHistogram(const TH1F* src, TH1F* dst, double shift, double intensity)
{
    dst->Reset("ICE");

    const int n      = src->GetNbinsX();
    const double min = src->GetXaxis()->GetXmin();
    const double max = src->GetXaxis()->GetXmax();

    for (int i = 1; i <= n; ++i){
        const double x  = src->GetBinCenter(i);
        const double xs = x - shift;

        double new_content = 0.0;
        // if(xs >= min && xs <= max) new_content = src->GetBinContent(src->GetXaxis()->FindBin(xs));
        if(xs >= min && xs <= max) new_content = src->Interpolate(xs);
        dst->SetBinContent(i, new_content);
    }
    dst->Scale(intensity);
}
// ─────────────────────────────────────────────────────────────────────────────
// LoadPDF — read templates from ROOT files
// Both templates are loaded from the same file (Ge68) for both PDF,
// ─────────────────────────────────────────────────────────────────────────────
void Fitter::LoadPDF()
{
    const std::string tag = (datatype == "MC") ? "MC" : "data";
    const std::string suffix = (hamaOnly) ? "_Hama" : "";
    const std::string foPsName = "/sps/juno/mlecocq/oPs/util/Ge68" + tag + "PDF_"+ std::to_string(nbin) +"bins" + suffix + ".root";
    // const std::string fIoniName = foPsName; // same source - change if needed
    const std::string fIoniName = "/sps/juno/mlecocq/oPs/util/"+ pdf_source + tag + "PDF_"+ std::to_string(nbin) +"bins" + suffix + ".root";


    // helper function that loads the PDFs for each component
    auto readHist = [&](const std::string& fname, const std::string& hname,
                        const std::string& label) -> std::unique_ptr<TH1F>
    {
        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie())
            throw std::runtime_error("LoadPDF: cannot open " + fname);

        TTree* tree = static_cast<TTree*>(f->Get("tree"));
        int nEvents = 0;
        if (tree) {
            tree->SetBranchAddress("nEvents", &nEvents);
            if(label == "oPs"){
                tree->SetBranchAddress("oPsScale", &oPsCoeff);
                tree->SetBranchAddress("oPsScale_err", &oPsCoeff_err);
            }
            tree->GetEntry(0);
        }

        TH1F* h = static_cast<TH1F*>(f->Get(hname.c_str()));
        if (!h) throw std::runtime_error("LoadPDF: histogram " + hname + " not found in " + fname);

        // Build a local re-binned copy on our axis
        auto out = std::make_unique<TH1F>(("_pdf_" + label).c_str(), "", nbin, xmin, xmax);
        out->SetDirectory(nullptr);

        for (int i = 1; i <= nbin; ++i) {
            int src_bin = h->FindBin(out->GetBinCenter(i));
            out->SetBinContent(i, h->GetBinContent(src_bin));
            out->SetBinError  (i, h->GetBinError  (src_bin));
        }

        const double integral = out->Integral();
        if (integral > 0.0) out->Scale(1.0 / integral); // normalised histogram to build PDF

        std::cout << "LoadPDF [" << label << "]: nEvents=" << nEvents
                  << "  integral (post-scale)=" << out->Integral() << "\n";

        // Coeff needed by DoFit() to estimate aoPs amplitude
        if (label == "oPs" && nEvents > 0) {
            // Integral to Energy conversion factor
            oPsCoeff = (integral / nEvents) / 1.022;
        }

        f->Close();
        return out;
    };

    hPDF_oPs = readHist(foPsName, "hTime_shift", "oPs");
    hPDF_ioni = readHist(fIoniName, "hTime_shift", "ioni");

    std::cout << "LoadPDF: Coeff=" << oPsCoeff << "\n";

    if (debug) {
        TCanvas c("_cPDF", "", 800, 500);
        hPDF_ioni->SetLineColor(kAzure - 4); hPDF_ioni->Draw("HIST");
        hPDF_oPs ->SetLineColor(kViolet - 5); hPDF_oPs->Draw("HIST SAMES");
        c.SaveAs((datatype + "_PDFs_check.root").c_str());
    }
}