#include "PositroniumFitter.h"

#include "TF1.h"
// Global pointer to the current fitter instance
static PositroniumFitter* instance = nullptr;

void FCNWrapper(int& npar, double* grad, double& fval, double* par, int flag) {
    if (instance)
        instance->FCN(npar, grad, fval, par, flag);  // Calls the real method
}

double fPDFWrapper(double *x, double* par){
    if(instance){
        return instance->fPDF(x, par);
    }
    else{
        return 0;
    }
}



PositroniumFitter::PositroniumFitter(TH1F* h, bool Debug, int firstbin, int lastbin, int nbins, std::string type)
    : hdata(h), debug(Debug), setparams(false), hama(false), xmin(firstbin), xmax(lastbin), nbin(nbins), datatype(type)
{
    // SetType("data");
    SetupMinuit();

    instance = this;
    std::cout << "Data type " << datatype << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    hPDF1 = new TH1F("hpdf1", "hpdf1", nbin, xmin, xmax);
    hPDF2 = new TH1F("hpdf2", "hpdf2", nbin, xmin, xmax);
    LoadPDF();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - begin;
    std::cout << "PDF Loaded in " << duration.count() << " seconds" <<std::endl;
    

    hIoni = (TH1F*)hPDF1->Clone("hIoni");
    hoPs = (TH1F*)hPDF2->Clone("hoPs");
    hFit = new TH1F("hFit", "Total Fit", nbin, xmin, xmax);
    std::cout << "hFit bins " << hFit->GetNbinsX() << std::endl;
}


PositroniumFitter::~PositroniumFitter()
{
    delete gMinuit;
}


void PositroniumFitter::SetupMinuit(){
    gMinuit = new TMinuit(5);
    if(debug) gMinuit->SetPrintLevel(1);
    else gMinuit->SetPrintLevel(-1);
    gMinuit->SetFCN(FCNWrapper);
}


double PositroniumFitter::PreScanT0Dt(double t0min, double t0max, double dtmin, double dtmax, int nt0, int ndt){
    std::cout << "\n[PreScan] 2D scan of (t0, dt)\n";

    prescan = true;

    double bestChi2 = 1e99;
    double bestT0 = 0.0;
    double bestDt = 0.0;

    // Fixed parameters during scan
    double par[5] = {0, 0, 0, 0, 0};
    double tmpErr, b1, b2;
    int ivar;
    TString name;

    gMinuit->mnpout(2, name, par[2], tmpErr, b1, b2, ivar);
    gMinuit->mnpout(3, name, par[3], tmpErr, b1, b2, ivar);
    gMinuit->mnpout(4, name, par[4], tmpErr, b1, b2, ivar);


    int npar = 5;
    double chi2;

    for (int i = 0; i < nt0; ++i) {
        double t0 = t0min + (t0max - t0min) * i / (nt0 - 1);

        for (int j = 0; j < ndt; ++j) {
            double dt = dtmin + (dtmax - dtmin) * j / (ndt - 1);

            par[0] = t0;
            par[1] = dt;

            FCNWrapper(npar, nullptr, chi2, par, 0);

            if (chi2 < bestChi2) {
                bestChi2 = chi2;
                bestT0 = t0;
                bestDt = dt;
            }
        }
    }

    std::cout << "[PreScan] Best point found:\n"
              << "  t0 = " << bestT0
              << ", dt = " << bestDt
              << ", chi2 = " << bestChi2 << std::endl;

    // Set starting values in Minuit
    int ierflg;
    double arglist[2];

    arglist[0] = 1;  // t0
    arglist[1] = bestT0;
    gMinuit->mnexcm("SET PAR", arglist, 2, ierflg);
    // gMinuit->FixParameter(0);

    arglist[0] = 2;  // dt
    arglist[1] = bestDt;
    gMinuit->mnexcm("SET PAR", arglist, 2, ierflg);

    prescan = false;

    return bestChi2;
}

void PositroniumFitter::setParameters(Params p){
    setparams = true;
    pInit = p;
}

void PositroniumFitter::DoFit(){

    double arglist[10];
    int ierflg;

    double Amp2 = 0.0, Amp1 = 0.0, tStart = 0.0, timediff = 0.0, bgnd = 0.0;

    
    if(!setparams){
        tStart = hdata->GetBinCenter(hdata->GetMaximumBin()) - hPDF1->GetBinCenter(hPDF1->GetMaximumBin());
        timediff = 2.0;
        // Estimate background from the first 150 ns
        TF1* f_bkg = new TF1("f_bkg", "[0]", 0, 150);
        hdata->Fit("f_bkg", "RQ0");
        bgnd = f_bkg->GetParameter(0);

        Amp2 = 1.022 * Coeff;  // oPs should correspond to 1.022 MeV
        Amp1 = (hdata->Integral() - Amp2 - bgnd*nbin) ; //Ioni is the rest
        
        // TSpectrum s(10);
        // int nFound = s.Search(hdata, 2, "nodraw", 0.7);
        // double xdata = s.GetPositionX()[0];

        // // nFound = s.Search(hPDF1, 2, "nodraw", 0.7);
        // tStart = xdata - s.GetPositionX()[0];

    }
    else{
        Amp2 = pInit.aoPs;  // Compute a renormalisation factor for oPs and Ioni
        Amp1 = pInit.aIoni ;
        tStart = pInit.t0;
        timediff = pInit.dt;
        bgnd = pInit.constant;
        setparams = false;
    }
    std::cout << "Info: Data Entries " << hdata->GetEntries() << " data integral " << hdata->Integral() << std::endl;
    std::cout << "Info: Starting values:" << " Amp 1 " << Amp1 << " Amp 2 " << Amp2 << " t0 " << tStart << " dt " << timediff << " bkg " << bgnd << std::endl;


    std::cout << "hdata nBins " << hdata->GetNbinsX() << " PDF nbins " << hPDF1->GetNbinsX() 
    << " first bin " << hPDF1->GetBinCenter(1) << " last bin " << hPDF1->GetBinCenter(hPDF1->GetNbinsX()) << std::endl;

    arglist[0] = 1.;
    gMinuit->mnexcm("SET ERR", arglist, 1,ierflg);

    gMinuit->mnparm(0, "t0", tStart, 0.005, tStart-0.3*tStart, tStart+0.3*tStart, ierflg); 
    gMinuit->mnparm(1, "dt", timediff, 0.005, -0.0, 30, ierflg);
    gMinuit->mnparm(2, "aoPs", Amp2, 100, Amp2 - 0.6*Amp2, Amp2 + 0.6*Amp2, ierflg);
    gMinuit->mnparm(3, "const", bgnd, 0.1, bgnd - 0.10*bgnd, bgnd + 0.10*bgnd, ierflg);
    gMinuit->mnparm(4, "aIoni", Amp1, 100, Amp1 - 0.6*Amp1, Amp1 + 0.6*Amp1, ierflg);


    // gMinuit->FixParameter(0); //t0
    // gMinuit->FixParameter(1); //dt
    // gMinuit->FixParameter(2); //aoPs
    // gMinuit->FixParameter(4); //aIoni
    gMinuit->FixParameter(3); //constant

// =============  2D Prescan of t0 and dt  =============

    double p[2];
    double tmpErr, b1, b2;
    int ivar;
    TString name;
    
    gMinuit->mnpout(0, name, p[0], tmpErr, b1, b2, ivar);
    gMinuit->mnpout(1, name, p[1], tmpErr, b1, b2, ivar);
    
    // First do a prescan of t0 and dt based on the chi2 value
    chi2_prescan = PreScanT0Dt(p[0]-0.3*p[0], p[0]+0.2*p[0], p[1]-0.9*p[1], p[1]+0.9*p[1], 20, 20);

// =============  start the fit  =============

    arglist[0] = 1;
    gMinuit->mnexcm("SET STR", arglist, 1, ierflg);

    gMinuit->mnsimp(); //Simplex
    
    arglist[0] = 500;
    arglist[1] = 1.0; //0.5 for NLL, 1.0 for Chi2
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    
    double amin, edm, errdef;
    int npari, nparx, istat;

    gMinuit->mnstat(amin, edm, errdef, npari, nparx, istat);

    if(istat == 3) gMinuit->mnexcm("MINOS", nullptr, 0, ierflg);

    std::cout << "Set Params " << setparams << std::endl; 
}

Params PositroniumFitter::FetchParams(){

    Params fitparameters; 

    // Get Chi2 and Chi2/NDF
    double amin, edm, errdef;
    int npari, nparx, istat;

    gMinuit->mnstat(amin, edm, errdef, npari, nparx, istat);

    fitparameters.istat = istat; 
    fitparameters.edm = edm;
    fitparameters.chi2 = amin;
    fitparameters.chi2ndf = amin / (nbin - npari); // 5 fitting parameters
    fitparameters.chi2prescan = chi2_prescan;

    // Get best fitting parameters

    TString name;
    double p[5];
    double ep[5];
    double bnd1[5];
    double bnd2[5];
    int varbl[5];

    gMinuit->mnpout(0, name, p[0], ep[0], bnd1[0], bnd2[0], varbl[0]); // t0
    gMinuit->mnpout(1, name, p[1], ep[1], bnd1[1], bnd2[1], varbl[1]); // dt
    gMinuit->mnpout(2, name, p[2], ep[3], bnd1[3], bnd2[3], varbl[3]); // aoPs
    gMinuit->mnpout(3, name, p[3], ep[4], bnd1[4], bnd2[4], varbl[4]); // c
    gMinuit->mnpout(4, name, p[4], ep[2], bnd1[2], bnd2[2], varbl[2]); // aIoni

    
    fitparameters.t0 = p[0]; 
    fitparameters.dt = p[1];
    fitparameters.aoPs = p[2];
    fitparameters.constant = p[3];
    fitparameters.aIoni = p[4]; //aIoni

    if(istat == 3){
        double errLow, errUp, errParab, gcc;

        gMinuit->mnerrs(0, errLow, errUp, errParab, gcc); //t0 error
        if(errUp == 0 && errLow == 0){
            fitparameters.t0_err[0] = errParab; 
            fitparameters.t0_err[1] = errParab;
        }
        else {
            fitparameters.t0_err[0] = errLow; 
            fitparameters.t0_err[1] = errUp;
        }
        gMinuit->mnerrs(1, errLow, errUp, errParab, gcc); //dt error
        if(errUp == 0 && errLow == 0){
            fitparameters.dt_err[0] = errParab; 
            fitparameters.dt_err[1] = errParab;
        }
        else {
            fitparameters.dt_err[0] = errLow; 
            fitparameters.dt_err[1] = errUp;
        }
        gMinuit->mnerrs(2, errLow, errUp, errParab, gcc); // aoPs error
        fitparameters.aoPs_err[0] = errLow; fitparameters.aoPs_err[1] = errUp;
        gMinuit->mnerrs(3, errLow, errUp, errParab, gcc); //c error
        fitparameters.constant_err[0] = errLow; fitparameters.constant_err[1] = errUp;
        gMinuit->mnerrs(4, errLow, errUp, errParab, gcc); //aIoni error
        fitparameters.aIoni_err[0] = errLow; fitparameters.aIoni_err[1] = errUp;
    }
    else{
        fitparameters.t0_err[0] = ep[0]; fitparameters.t0_err[1] = ep[0]; //t0 error
        fitparameters.dt_err[0] = ep[1]; fitparameters.dt_err[1] = ep[1]; //dt error
        fitparameters.aoPs_err[0] = ep[3]; fitparameters.aoPs_err[1] = ep[3]; // aoPs error
        fitparameters.constant_err[0] = ep[4]; fitparameters.constant_err[1] = ep[4]; //c error
        fitparameters.aIoni_err[0] = ep[2]; fitparameters.aIoni_err[1] = ep[2]; //aIoni error
    }


    // =============  Save correlation matrix  =============
    std::cout << "Number of parameters " << npari << std::endl;

    std::vector<double> cov(npari*npari);

    gMinuit->mnemat(&cov[0], npari);

    for (int i = 0; i < npari; ++i) {
        for (int j = 0; j < npari; ++j) {
            double denom = std::sqrt(cov[i*npari + i] * cov[j*npari + j]);
            fitparameters.corr[i*npari + j] = (denom > 0.0) ? cov[i*npari + j] / denom : 0.0;
        }
    }
    

    shiftHistogram(hPDF1, hIoni, p[0], p[4]);
    shiftHistogram(hPDF2, hoPs, p[0]+p[1], p[2]);
    hFit->Reset();
    hFit->Add(hIoni);
    hFit->Add(hoPs);


    fitparameters.hIoni = hIoni;
    fitparameters.hoPs = hoPs;
    fitparameters.hFit = hFit;

    std::cout << "Info: Final fit values:" << std::endl;
    std::cout << " - aIoni: " << fitparameters.aIoni << ", aoPs: " << fitparameters.aoPs << ", dt: " << fitparameters.dt << std::endl;
    std::cout << " - aIoni_err_low: " << fitparameters.aIoni_err[0] << ", aoPs_err_low: " << fitparameters.aoPs_err[0] << ", dt_err_low: " << fitparameters.dt_err[0] << std::endl;
    std::cout << " - aIoni_err_up: " << fitparameters.aIoni_err[1] << ", aoPs_err_up: " << fitparameters.aoPs_err[1] << ", dt_err_up: " << fitparameters.dt_err[1] << std::endl;
    
    std::cout << "Info: Quality of fit:\n" <<
    " - istat: " << istat << " edm: " << edm << std::endl;

    return fitparameters;
}

void PositroniumFitter::LoadPDF(){

    int nEvents;
    double DCR, errDCR;

    std::string fname;

    //2.2MeV Gamma pdf
    std::cout << "\nInfo: Loading Gamma PDF" << std::endl;
    // if(datatype == "data") fname = "/sps/juno/mlecocq/oPs/util/GammadataPDF_"+std::to_string(nbin)+"bins.root";
    if(datatype == "data") fname = "/sps/juno/mlecocq/oPs/util/Ge68dataPDF_"+std::to_string(nbin)+"bins.root";
    // if(datatype == "MC") fname = "/sps/juno/mlecocq/oPs/util/GammaMCPDF_"+std::to_string(nbin)+"bins.root";
    if(datatype == "MC") fname = "/sps/juno/mlecocq/oPs/util/Ge68MCPDF_"+std::to_string(nbin)+"bins.root";

    TFile* f = TFile::Open(fname.c_str(), "READ");

    TTree* t = (TTree*)f->Get("tree");
    t->SetBranchStatus("nEvents", 1);
    t->SetBranchAddress("nEvents", &nEvents);

    t->GetEntry(0);

    TH1F* h1 = (TH1F*)f->Get("hTime");
    // TH1F* h1 = (TH1F*)f->Get("htempTot");

    std::cout << "xmin " << xmin << " xmax " << xmax << std::endl;

    for (int i = 1; i <= hPDF1->GetNbinsX(); i++){
        int x = hPDF1->GetBinCenter(i);
        int source_bin = h1->FindBin(x);

        double content = h1->GetBinContent(source_bin);
        double error = h1->GetBinContent(source_bin);

        hPDF1->SetBinContent(i, content);
        hPDF1->SetBinError(i, error);
    }

    hPDF1->Sumw2();
    hPDF1->Scale(1.0/hPDF1->Integral());
    std::cout << "PDF build with " << nEvents << " events"<< std::endl;
    std::cout << "PDF1 integral " << hPDF1->Integral() << std::endl;


    // 1.022 MeV Gamma pdf
    std::cout << "\nInfo: Loading Ge68 PDF" << std::endl;
    if(datatype == "data") fname = "/sps/juno/mlecocq/oPs/util/Ge68dataPDF_"+std::to_string(nbin)+"bins.root";
    if(datatype == "MC") fname = "/sps/juno/mlecocq/oPs/util/Ge68MCPDF_"+std::to_string(nbin)+"bins.root";
    // if(datatype == "MC") fname = "/sps/juno/mlecocq/oPs/util/GammaMCPDF.root";
    f = TFile::Open(fname.c_str(), "READ");
    
    t = (TTree*)f->Get("tree");
    t->SetBranchStatus("nEvents", 1);
    t->SetBranchAddress("nEvents", &nEvents);

    t->GetEntry(0);


    h1 = (TH1F*)f->Get("hTime");
    int pdfbins = h1->GetNbinsX();
    // h1 = (TH1F*)f->Get("htempTot");

    double norm = pdfbins/nbin;
    std::cout << "Info: pdf bins " << pdfbins << " data bins " << nbin << " norm " << norm <<std::endl;
    std::cout << "Info: xmin " << xmin << " xmax " << xmax << std::endl;
    
    for (int i = 1; i <= hPDF2->GetNbinsX(); i++){
        int x = hPDF2->GetBinCenter(i);
        int source_bin = h1->FindBin(x);
        
        double content = h1->GetBinContent(source_bin);
        double error = h1->GetBinError(source_bin);
        
        hPDF2->SetBinContent(i, content);
        hPDF2->SetBinError(i, error);
    }
    double PDFEntries = hPDF2->Integral();

    hPDF2->Sumw2();
    hPDF2->Scale(1.0/hPDF2->Integral());

    double Energy = 1.022; //MeV
    double MeanTime = h1->GetMean();

    Coeff = (PDFEntries/nEvents) / Energy;
    Coeff*=norm; // Normalise coeff according to the binning (1.0 if same binning)
    std::cout << "PDF Integral: " << PDFEntries << " Events: " << nEvents << " ratio: " << PDFEntries/nEvents << " Coeff: " << Coeff << std::endl;
    std::cout << "PDF integral " << hPDF2->Integral() << std::endl;

    if(debug){
        TCanvas* c = new TCanvas();
        hPDF1->Draw();
        hPDF2->SetLineColor(kMagenta+2);
        hPDF2->Draw("sames");
        c->SaveAs(Form("%sPDFs.root", datatype.c_str()));
    }
    f->Close();
}





void PositroniumFitter::shiftHistogram(TH1F* h, TH1F* h_shift, double shift, double intensity){
    h_shift->Reset(); // Clear output histogram

    int nbins = h->GetNbinsX();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    for (int i = 1; i <= nbins; ++i) {
        double x = h_shift->GetXaxis()->GetBinCenter(i);
       
        if(x < xmin || x > xmax){
            continue;
        } 
        
        double shiftedX = x - shift; 

        double content = 0.0;
        if (shiftedX >= xmin && shiftedX <= xmax) {
            content = h->GetBinContent(h->FindBin(shiftedX)); 
            // content = h->Interpolate(shiftedX); 
        }

        h_shift->SetBinContent(i, intensity * content);
    }

}



void PositroniumFitter::FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
    instance->loop+=1;

    int evt = instance->loop;
    double chi2n = 0.0;
    
    double t0 = par[0];
    double dt = par[1];
    double a2 = par[2]; //aoPs
    double c = par[3];
    // double a1 = instance->hdata->Integral() - a2 - c*instance->nbin; //aIoni
    double a1 = par[4];

    instance->shiftHistogram(instance->hPDF1, instance->hIoni, t0, a1);
    instance->shiftHistogram(instance->hPDF2, instance->hoPs, t0 + dt, a2);

    instance->hFit->Reset();

    int nbins = instance->hdata->GetNbinsX();
    int effective_bins = 0;
    for (int i = 1; i <= nbins; i++){
        
        double t = instance->hdata->GetBinCenter(i);

        // if (t < 200 || t > 600) continue;
        double val = instance->hIoni->GetBinContent(i) + instance->hoPs->GetBinContent(i) + c;
        // double val = instance->hIoni->Interpolate(t) + instance->hoPs->Interpolate(t) + c;

        if(val<=0 ) val = 0.01;
        if(instance->hdata->GetBinContent(i) > 0){    
            chi2n += 2*(val - instance->hdata->GetBinContent(i) + instance->hdata->GetBinContent(i)*log(instance->hdata->GetBinContent(i)/val));
        }
        else{
            chi2n += 2*val; 
        }
        int ii =  instance->hFit->FindBin(t);
        instance->hFit->SetBinContent(ii, val);
    }
    if(!instance->prescan){
        double penalty2 = std::pow((a2 - 1.022*instance->Coeff)/0.4, 2);
        chi2n += penalty2;
        
        // double penalty = std::pow((a1 - (instance->hdata->Integral() - a2 - c*instance->nbin))/0.6, 2);
        // double penalty = std::pow((instance->hIoni->GetMaximum() - instance->hdata->GetMaximum())/0.2, 2);
        // chi2n += penalty;
    }



    // if(instance->loop % 1 == 0){

    //     TCanvas* c = new TCanvas();
    //     instance->hdata->SetLineColor(kBlack);
    //     instance->hdata->SetMarkerStyle(20);
    //     instance->hdata->SetMarkerColor(kBlack);
    //     instance->hdata->Draw("E1X0");
    //     instance->hdata->GetXaxis()->SetRangeUser(150, 600);
        
    //     instance->hFit->SetMarkerColor(kOrange+7);
    //     instance->hFit->Draw("sames");
    //     instance->hIoni->SetMarkerColor(kAzure-4);
    //     instance->hIoni->Draw("sames");
    //     instance->hoPs->SetMarkerColor(kMagenta+2);
    //     instance->hoPs->Draw("sames");
    //     c->SaveAs(Form("png/Fit_%i.png", instance->loop));
    // }

    f = chi2n;
}

double PositroniumFitter::fPDF(double *x, double *par){
    double xx = x[0];
    instance->shiftHistogram(instance->hPDF1, instance->hIoni, par[0], par[2]); //t0 and aIoni
    instance->shiftHistogram(instance->hPDF2, instance->hoPs, par[1], par[3]); //t0+dt and aoPs

    double c = par[4];
    int binIoni = instance->hIoni->FindBin(xx);
    int binoPs = instance->hoPs->FindBin(xx);

    double yIoni = instance->hIoni->GetBinContent(binIoni);
    double yoPs = instance->hoPs->GetBinContent(binoPs);

    return yIoni + yoPs + c;
}
