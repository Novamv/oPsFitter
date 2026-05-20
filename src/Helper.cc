#include "Helper.h"


// ─────────────────────────────────────────────────────────────────────────────
// Constructor / Destructor
// ─────────────────────────────────────────────────────────────────────────────
Helper::Helper(bool debug_, int min_, int max_, int nbin_, std::string datatype_)
    : debug(debug_), hamaOnly(false), xmin(min_), xmax(max_), nbin(nbin_), datatype(datatype_),
      Sim("Sim"), Calib("Calib"), Reco("Reco"), IBD("event")
{
    LoadPMTInfo();
}

Helper::~Helper()
{}

// ─────────────────────────────────────────────────────────────────────────────
// Progress Bar
// ─────────────────────────────────────────────────────────────────────────────
void Helper::ProgressBar(Long64_t current, Long64_t total, int width){
    float progress = float(current) / total;
    int pos = width * progress;

    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

// ─────────────────────────────────────────────────────────────────────────────
// PMT Information (Type, ID, Circle...)
// ─────────────────────────────────────────────────────────────────────────────
void Helper::LoadPMTInfo(){
    TChain PMT("PMTInfo");
    PMT.Add("/sps/juno/mlecocq/PMTInfoLatest.root");

    Long64_t pmtid;
    char type[256];

    PMT.SetBranchStatus("*", 0);
    PMT.SetBranchStatus("CopyNo", 1);
    PMT.SetBranchAddress("CopyNo", &pmtid);
    PMT.SetBranchStatus("PMTType", 1);
    PMT.SetBranchAddress("PMTType", type);

    for (Long64_t i = 0; i < PMT.GetEntries(); ++i){
        PMT.GetEntry(i);
        PmtType[pmtid] = std::string(type);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PDF Handler (Generate, Load)
// ─────────────────────────────────────────────────────────────────────────────
double Helper::shift(TH1F* src, TH1F* dst, double TargetTime){
    
    const int n      = src->GetNbinsX();
    const double min = src->GetXaxis()->GetXmin();
    const double max = src->GetXaxis()->GetXmax();

    // ── Fit sigmoid to the rising edge only ────────────────────────
    // Search for rise region: from first non-zero bin to the peak
    double peak_x     = src->GetBinCenter(src->GetMaximumBin());
    double peak_val   = src->GetMaximum();

    // Fit window: a few ns before the rise to the peak
    double fit_lo = peak_x - 20.0;
    double fit_hi = peak_x + 5.0;

    // f(x) = A / (1 + exp(-k*(x - x0)))
    // par[0] = A    (amplitude ~ peak value)
    // par[1] = k    (steepness of rise)
    // par[2] = x0   (inflection point — our alignment reference)
    TF1 fsig("_fsig", "[0]/(1+exp(-[1]*(x-[2])))", fit_lo, fit_hi);
    fsig.SetParameter(0, peak_val);
    fsig.SetParameter(1, 1.0);      // positive: rising left to right
    fsig.SetParameter(2, peak_x - 5.0);
    fsig.SetParLimits(1, 0.1, 10.0);  // enforce rising slope

    src->Fit(&fsig, "RQ0");

    double x_inflection = fsig.GetParameter(2);

    // Fallback: if fit failed or gave nonsense, use peak bin
    if(x_inflection < min || x_inflection > peak_x){
        x_inflection = peak_x;
    }

    // ── Shift so inflection point lands on TargetTime ──────────────
    const double delta = x_inflection - TargetTime;

    dst->Reset("ICE");
    for (int i = 1; i <= n; ++i){
        const double t  = dst->GetBinCenter(i);
        const double xs = t + delta;

        double new_content = 0.0;
        if(xs >= min && xs <= max) new_content = src->GetBinContent(src->FindBin(xs));
        dst->SetBinContent(i, new_content);
    }

    return x_inflection;
}

bool Helper::CreateGe68PDF(const std::string& fname){
    
    std::cout << "Info: Generating Ge68 PDF with data type " << datatype << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    double ChargeTotLpmt, ChargeTotWP;
    float  RECX, RECY, RECZ;
    ULong64_t TimeStamp;
    std::vector<std::string>* TriggerName = nullptr;
    std::vector<double>* HitTime = nullptr;
    std::vector<int>* PmtID = nullptr;

    TChain Ge68Calib("Calib");
    // TChain Ge68Reco("Reco");
    TChain Ge68Reco("Reco");

    Ge68Calib.Add(fname.c_str());
    Ge68Reco.Add(fname.c_str());
    Ge68Calib.GetListOfFiles()->Print();

    Ge68Calib.SetBranchStatus("*", 0);
    Ge68Calib.SetBranchStatus("ChargeTotLPMT", 1);
    Ge68Calib.SetBranchStatus("ChargeTotWP", 1);
    Ge68Calib.SetBranchStatus("TimeStamp", 1);
    Ge68Calib.SetBranchStatus("TriggerName", 1);
    Ge68Calib.SetBranchStatus("HitTimeCalibTOF", 1);
    Ge68Calib.SetBranchStatus("PmtIDCalib", 1);
    
    Ge68Calib.SetBranchAddress("ChargeTotLPMT", &ChargeTotLpmt);
    Ge68Calib.SetBranchAddress("ChargeTotWP", &ChargeTotWP);
    Ge68Calib.SetBranchAddress("TimeStamp", &TimeStamp);
    Ge68Calib.SetBranchAddress("TriggerName", &TriggerName);
    Ge68Calib.SetBranchAddress("HitTimeCalibTOF", &HitTime);
    Ge68Calib.SetBranchAddress("PmtIDCalib", &PmtID);



    Ge68Reco.SetBranchStatus("*", 0);
    Ge68Reco.SetBranchStatus("Recx", 1);
    Ge68Reco.SetBranchStatus("Recy", 1);
    Ge68Reco.SetBranchStatus("Recz", 1);

    Ge68Reco.SetBranchAddress("Recx", &RECX);
    Ge68Reco.SetBranchAddress("Recy", &RECY);
    Ge68Reco.SetBranchAddress("Recz", &RECZ);


    // =============  Define some variables  =============
    
    TH1F* h = new TH1F("hTime", "hTime", nbin, xmin, xmax);

    int nEntries = Ge68Calib.GetEntries();
    std::cout << "Total entries " << nEntries << std::endl;
    int nEvents = 0;
    int nHits = 0;

    TH1F* hRef = (TH1F*)h->Clone("hRef");
    TH1F* htemp = (TH1F*)h->Clone("htemp");
    TH1F* htempTot = (TH1F*)h->Clone("htempTot");
    TH1F* htemp_shift = (TH1F*)h->Clone("htemp_shift");
    TH1F* htemp_shiftTot = (TH1F*)h->Clone("htemp_shiftTot");

    TH1F* hCharge = new TH1F("hCharge", "hCharge", 100, 0, 3000);

    float lo_bound = 1000;
    float hi_bound = 2500;
    if(hamaOnly) {
        lo_bound = 350; 
        hi_bound = 500;
    }
    TH1F* hExpected = new TH1F("hExpected", "hExpected", 150, lo_bound, hi_bound);

    TCanvas* c1 = new TCanvas("cShift1", "Shift example", 800, 500);
    TCanvas* c2 = new TCanvas("cShift2", "Accumulated shift", 800, 500);

    // =============  Get charge spectrum  =============
    
    for(size_t i = 0; i < nEntries; i++) {
        Ge68Calib.GetEntry(i);
        hCharge->Fill(ChargeTotLpmt);
    }

    float mean_q = hCharge->GetBinCenter(hCharge->GetMaximumBin());
    float sig_q = hCharge->GetStdDev();
    float lo = mean_q - 1.5*sig_q;
    float hi = mean_q + 1.5*sig_q;

    TF1* fgaus = new TF1("fgaus", "gaus", lo, hi);
    hCharge->Fit("fgaus", "R");
    mean_q = fgaus->GetParameter(1);
    sig_q = fgaus->GetParameter(2);

    TCanvas* cCharge = new TCanvas("cCharge", "cCharge", 800, 500);
    hCharge->Draw();
    fgaus->SetLineColor(kRed+1);
    fgaus->Draw("same");
    cCharge->SaveAs("Ge68data_charge.root");

    // =============  Generate PDF  =============

    std::cout << "  - Hama Only: " << hamaOnly << " suffix: " << suffix << std::endl;
    
    uint64_t tLastMuon;
    int n = -1;
    for(size_t i = 0; i < Ge68Calib.GetEntries(); i++){

        if (i % 1000 == 0 || i == nEntries - 1) {
            ProgressBar(i + 1, nEntries);
        }

        Ge68Calib.GetEntry(i);
        if(!(ChargeTotLpmt >= mean_q-3*sig_q && ChargeTotLpmt <= mean_q+3*sig_q)) continue;

        Ge68Reco.GetEntry(i);
        
        double dt = 0.0;
        float rho = 0.0;
        bool good_event = false;
        if(datatype == "data"){
            if(ChargeTotLpmt >= 30000 || ChargeTotWP >= 700) {
                tLastMuon = TimeStamp;
            }
            double dt = TimeStamp - tLastMuon;
            float rho = std::sqrt(RECX*RECX + RECY*RECY);
            good_event = dt*1e-6 > 2 && rho < 750 && std::abs(RECZ) < 500 && ChargeTotLpmt !=0;
        } else {
            good_event = rho < 750 && std::abs(RECZ) < 500 && ChargeTotLpmt !=0;
        }


        if(good_event){
            n++;
            nEvents++;
            htemp->Reset("ICE");
            for (size_t hit = 0; hit < HitTime->size(); hit++){

                if(!hamaOnly && PmtID->at(hit) > 17612) continue;
                if(hamaOnly && (PmtID->at(hit) < 17612 && PmtType[PmtID->at(hit)] != "Hamamatsu")) continue;

                htemp->Fill(HitTime->at(hit));
            }
        
            double inflection = shift(htemp, htemp_shift); //shift histogram
            h->Add(htemp_shift);
            htempTot->Add(htemp);
            htemp_shiftTot->Add(htemp_shift);

            hExpected->Fill(htemp->Integral()/1.022);
            
            if(n == 0) hRef->Add(htemp_shift);
            
            if(n == 5){
                int maxy = hRef->GetMaximum()*1.1;
                int miny = maxy/3;
                gStyle->SetOptStat(0);

                c1->cd();
                hRef->SetLineColor(kCyan+2);
                hRef->Draw("HIST");
                hRef->SetTitle("Hit Time correction;Time [ns];Entries");
                htemp->Draw("HIST SAME");
                htemp_shift->SetLineColor(kMagenta+2);
                htemp_shift->Draw("hist same");

                TLegend* leg = new TLegend(0.50, 0.25, 0.85, 0.35);
                leg->AddEntry(hRef, "Reference Histogram", "l");
                leg->AddEntry(htemp, "Original Hit Time", "l");
                leg->AddEntry(htemp_shift, "Shifted Hit Time", "l");
                leg->SetBorderSize(2); 
                leg->Draw();
                
                TPad* zoomPad = new TPad("zoomPad", "Zoomed View", 0.40, 0.50, 0.90, 0.80);     
                gStyle->SetOptStat(0);                                   
                zoomPad->SetFillStyle(0);     // transparent                                                                 
                zoomPad->SetFrameLineWidth(1);                                                                               
                zoomPad->Draw();                                                                                             
                zoomPad->cd();                                                                                               
                TH1F* frame = zoomPad->DrawFrame(150,miny,350,maxy);
                frame->GetXaxis()->SetLabelSize(0.08);                                                                              
                frame->GetXaxis()->SetTitleSize(0.08);                                                                              
                frame->GetYaxis()->SetLabelSize(0.08); 
                frame->GetYaxis()->SetTitleSize(0.08);
                frame->SetTitle("Zoomed view on peak;Time [ns];Entries");                                                                               
                hRef->Draw("SAME");                                                                                            
                htemp->Draw("SAME"); 
                htemp_shift->Draw("SAME");
            }

            TH1F* htemp_clone = (TH1F*)htemp_shift->Clone(Form("htemp_clone_%d", int(i)));
            c2->cd();
            htemp_clone->Draw("PLC SAME");
        }
    }

    TF1* f = new TF1("f", "[0]", 50, 200);
    f->SetParameter(0, 4000);
    h->Fit("f", "R");

    for(int bin = 1; bin <= h->GetNbinsX(); bin++){
        double val = h->GetBinContent(bin);
        double content = val - f->GetParameter(0);
        if (content < 0) content = 0;
        h->SetBinContent(bin, content);
        htempTot->SetBinContent(bin, htempTot->GetBinContent(bin) - f->GetParameter(0));
    }

    
    f = new TF1("f_gauss", "gaus", lo_bound, hi_bound);
    hExpected->Fit(f, "R");

    double oPsExpected = f->GetParameter(1);
    double oPs_err = f->GetParError(1);

    std::cout << "Creating PDF file" << std::endl;
    // h->Scale(1.0/h->Integral());
    TH1F* hPDF = (TH1F*)h->Clone("hPDF");
    hPDF->Scale(1.0/hPDF->Integral());

    TTree* Tree = new TTree("tree", "Ge68 Tree");
    Tree->Branch("nEvents", &nEvents);
    Tree->Branch("nHits", &nHits);
    Tree->Branch("oPsScale", &oPsExpected);
    Tree->Branch("oPsScale_err", &oPs_err);

    Tree->Fill();

    TFile* fout = new TFile(Form("/sps/juno/mlecocq/oPs/util/Ge68%sPDF_%ibins%s.root", datatype.c_str(), nbin, suffix.c_str()), "RECREATE");
    fout->cd();
    Tree->Write();
    hCharge->Write();
    h->Write();
    hPDF->Write();
    htempTot->Write();
    htemp_shiftTot->Write();
    hExpected->Write();
    c1->Write();
    c2->Write();
    fout->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - begin;

    std::cout << "Info: Ge68 PDF generated in " << duration.count() << " mins" << std::endl;

    hGe68 = (TH1F*)hPDF->Clone();
    delete h;
    return true;
}

bool Helper::GeneratePDF(){

    bool pdf = false;
    if(hamaOnly) {
        suffix = "_Hama";
        std::cout << "Info: Building PDF with Hamamatsu only. File Suffix:" << suffix << std::endl;
    }
    else suffix="";

    std::string fGe68Name;
    if(datatype == "MC") fGe68Name = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/Ge68_1000Evts_"+juno_version+"_*_plainReco.root";
    else fGe68Name = "/sps/juno/JunoData/Calib/2025/0824/RUN.9541.JUNODAQ.Calib-ACU-Ge68-Global-Pos-x0y0z0*_plain.root";
    
    pdf = CreateGe68PDF(fGe68Name);
    return pdf;
}

// ─────────────────────────────────────────────────────────────────────────────
// Load Data and Return data hist
// ─────────────────────────────────────────────────────────────────────────────
void Helper::LoadData(bool delay){
    std::string fname = "";

    if(datatype == "MC"){
        Sim.Reset();
        Calib.Reset();
        std::string prefix = "/sps/juno/mlecocq/oPs/MC/plain/"+ juno_version + "/Uniform/positronium_1000Evts_"+energy+"MeV_";
        if (delay) prefix = "/sps/juno/mlecocq/oPs/MC/plain/"+ juno_version + "/Uniform/gamma_1000Evts_";

        //Sim
        std::cout << "\n Info: Sim" << std::endl;
        fname = prefix + juno_version + "*_plainSim.root";
        Sim.Add(fname.c_str());
        Sim.GetListOfFiles()->Print();
        Sim.SetBranchStatus("*", 0);
        Sim.SetBranchStatus("EntryNb", 1);
        Sim.SetBranchStatus("PDGID",   1);
        Sim.SetBranchStatus("InitT",   1);

        Sim.SetBranchAddress("EntryNb", &SimEvtID);
        Sim.SetBranchAddress("PDGID",   &PDGID);
        Sim.SetBranchAddress("InitT",   &InitT);

        std::cout << "Sim Entries " << Sim.GetEntries() << std::endl;
        for (size_t i = 0; i < Sim.GetEntries(); ++i){
            Sim.GetEntry(i);
            double oPs_dt = 0.0;
            for(int id = 0; id < PDGID->size(); ++id){
                if(PDGID->at(id) == 22){
                   oPs_dt = InitT->at(id);
                    break;
                }
            }
            oPsTrueDecay[i] = oPs_dt;
        }
        
        // Calib
        std::cout << "\n Info: Calib" << std::endl;
        fname = prefix + juno_version + "*_plainReco.root";
        Calib.Add(fname.c_str());
        Calib.GetListOfFiles()->Print();
        Calib.SetBranchStatus("*",               0);
        Calib.SetBranchStatus("EntryNb",         1);
        Calib.SetBranchStatus("PmtIDCalib",      1);
        Calib.SetBranchStatus("HitTimeCalibTOF", 1);

        Calib.SetBranchAddress("EntryNb",         &EntryNb);
        Calib.SetBranchAddress("PmtIDCalib",      &PmtIDCalib);
        Calib.SetBranchAddress("HitTimeCalibTOF", &HitTimeTOF);
        // Calib.SetBranchStatus("HitTimeCalibTrueTOF", 1);
        // Calib.SetBranchAddress("HitTimeCalibTrueTOF", &PromptTOF);

        Reco.Add(fname.c_str());
        Reco.SetBranchStatus("*",         0);
        Reco.SetBranchStatus("Recx",      1);
        Reco.SetBranchStatus("Recy",      1);
        Reco.SetBranchStatus("Recz",      1);
        Reco.SetBranchStatus("RecEnergy", 1);


        Reco.SetBranchAddress("RecEnergy", &RecEnergy);
        Reco.SetBranchAddress("Recx",      &Recx);
        Reco.SetBranchAddress("Recy",      &Recy);
        Reco.SetBranchAddress("Recz",      &Recz);

        nEntries = Calib.GetEntries();
        std::cout << "Calib Entries " << nEntries << std::endl;
    }else{
        fname = "/sps/juno/mlecocq/Data/Physics/Processing/IBD/ReProd25C/*.root";
        IBD.Reset();
        IBD.Add(fname.c_str());
        IBD.GetListOfFiles()->Print();
    
        IBD.SetBranchStatus("*",               0);
        IBD.SetBranchStatus("Tag",             1);
        IBD.SetBranchStatus("RecEnergy",       1);
        IBD.SetBranchStatus("Recx",            1);
        IBD.SetBranchStatus("Recy",            1);
        IBD.SetBranchStatus("Recz",            1);
        IBD.SetBranchStatus("HitTimeCalibTOF", 1);
        IBD.SetBranchStatus("PmtIDCalib",      1);
    
        IBD.SetBranchAddress("Tag",             &Tag);
        IBD.SetBranchAddress("RecEnergy",       &RecEnergy);
        IBD.SetBranchAddress("Recx",            &Recx);
        IBD.SetBranchAddress("Recy",            &Recy);
        IBD.SetBranchAddress("Recz",            &Recz);
        IBD.SetBranchAddress("PmtIDCalib",      &PmtIDCalib);
        IBD.SetBranchAddress("HitTimeCalibTOF", &HitTimeTOF);

        nEntries = IBD.GetEntries();

    }
}
