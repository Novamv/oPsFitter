#include "PDFManager.h"

void ProgressBar(Long64_t current, Long64_t total, int width = 50){
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

PDFManager::PDFManager()
    : data_type("data"), TargetTime(250), xmin(0), xmax(1000), nbin(500), juno_version("J24.1.2"), HamaOnly(false)
{
}

PDFManager::~PDFManager(){}

bool PDFManager::CreateGammaPDF(const std::string& fname) {

    std::cout << "Info: Generating Gamma PDF" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<double>* HitTime = nullptr;
    std::vector<int>* PmtID = nullptr;

    TChain* Gamma;
    if(data_type == "MC") Gamma = new TChain("Calib");
    else if(data_type == "data") Gamma = new TChain("ibd");

    Gamma->Add(fname.c_str());
    Gamma->GetListOfFiles()->Print();

    Gamma->SetBranchStatus("*", 0);
    if(data_type == "MC"){
        Gamma->SetBranchStatus("HitTimeCalibTOF", 1);        
        Gamma->SetBranchAddress("HitTimeCalibTOF", &HitTime);
        Gamma->SetBranchStatus("PmtIDCalib", 1);
        Gamma->SetBranchAddress("PmtIDCalib", &PmtID);
    }
    else if(data_type == "data"){
        Gamma->SetBranchStatus("HitTime_d", 1);        
        Gamma->SetBranchAddress("HitTime_d", &HitTime);
    }


    TH1F* h = new TH1F("hTime", "hTime", nbin, xmin, xmax);

    int nEntries = Gamma->GetEntries();
    int nEvents = 0;
    int nHits = 0;

    TH1F* hRef = (TH1F*)h->Clone("hRef");
    TH1F* htemp = (TH1F*)h->Clone("htemp");
    TH1F* htempTot = (TH1F*)h->Clone("htempTot");
    TH1F* htemp_shift = (TH1F*)h->Clone("htemp_shift");
    TH1F* htemp_shiftTot = (TH1F*)h->Clone("htemp_shiftTot");

    TCanvas* c1 = new TCanvas("cShift1", "Shift example", 800, 500);
    TCanvas* c2 = new TCanvas("cShift2", "Accumulated shift", 800, 500);
    
    std::cout << "Filling histogram" << std::endl;

    gStyle->SetOptStat(0);

    for (size_t i = 0; i < nEntries; i++){

        if (i % 1000 == 0 || i == nEntries - 1) {
            ProgressBar(i + 1, nEntries);
        }

        Gamma->GetEntry(i);
        nEvents++;
        
        htemp->Reset();
        for (size_t hit = 0; hit < HitTime->size(); hit++){
            if(data_type=="MC" && PmtID->at(hit) > 17612) continue;
            nHits++;
            htemp->Fill(HitTime->at(hit));
        }
        
        double max = shift(htemp, htemp_shift, TargetTime); //shift histogram and return max value of this event
        h->Add(htemp_shift);
        htempTot->Add(htemp);
        htemp_shiftTot->Add(htemp_shift);

        if(i==0) hRef->Add(htemp_shift);
        
        if(i == 7 || i == 30){ //visualise shifting for some events

            int maxy = hRef->GetMaximum()*1.2;
            int miny = maxy/3;
            gStyle->SetOptStat(0);

            c1->cd();
            hRef->SetLineColor(kCyan+2);
            hRef->Draw("HIST");
            hRef->SetTitle("Hit Time correction;Time [ns];Entries");
            htemp->Draw("HIST SAME");
            htemp_shift->SetLineColor(kMagenta+2);
            htemp_shift->Draw("HIST SAME");
            
            TLegend* leg = new TLegend(0.50, 0.25, 0.85, 0.35);
            leg->AddEntry(hRef, "Reference Histogram", "l");
            leg->AddEntry(htemp, "Original Hit Time", "l");
            leg->AddEntry(htemp_shift, "Shifted Hit Time", "l");
            leg->SetBorderSize(3); 
            leg->Draw();

            TPad* zoomPad = new TPad("zoomPad", "Zoomed View",0.40, 0.50, 0.90, 0.80);     
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
        // TH1F* htemp_clone = (TH1F*)htemp_shift->Clone(Form("htemp_clone_%d", int(i)));
        // c2->cd();
        // htemp_clone->Draw("PLC SAME");
    }

    TF1* f = new TF1("f", "[0]", 0, 200);
    h->Fit("f", "R");

    for(int bin = 1; bin <= h->GetNbinsX(); bin++){
        double val = h->GetBinContent(bin);
        double content = val - f->GetParameter(0);
        if (content < 0) content = 0;
        h->SetBinContent(bin, val - f->GetParameter(0));
    }

    std::cout << "Creating PDF file" << std::endl;
    // h->Scale(1.0/h->Integral());
    TH1F* hPDF = (TH1F*)h->Clone("hPDF");
    hPDF->Scale(1.0/hPDF->Integral("width"));

    TTree* Tree = new TTree("tree", "Gamma Tree");
    Tree->Branch("nEvents", &nEvents);
    Tree->Branch("nHits", &nHits);

    Tree->Fill();

    TFile* fout = new TFile(Form("/sps/juno/mlecocq/oPs/util/Gamma%sPDF_%ibins.root", data_type.c_str(), nbin), "RECREATE");
    fout->cd();
    Tree->Write();
    h->Write();
    hPDF->Write();
    htempTot->Write();
    htemp_shiftTot->Write();
    c1->Write();
    // c2->Write();
    fout->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - begin;

    std::cout << "Info: Gamma PDF generated in " << duration.count() << " mins" << std::endl;

    hGamma = (TH1F*)h->Clone();
    delete h;
    delete c1;
    delete c2;

    return true;
}



bool PDFManager::CreateSimGe68PDF(const std::string& fname) {
    std::cout << "Info: Generating Ge68 PDF" << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();
    std::vector<double>* HitTime = nullptr;
    std::vector<int>* PmtID = nullptr;
    double TotalPE;

    TChain Ge68("Calib");
    Ge68.Add(fname.c_str());
    Ge68.GetListOfFiles()->Print();

    Ge68.SetBranchStatus("*", 0);
    Ge68.SetBranchStatus("ChargeTotLPMT", 1);
    Ge68.SetBranchStatus("HitTimeCalibTOF", 1);
    Ge68.SetBranchStatus("PmtIDCalib", 1);

    Ge68.SetBranchAddress("ChargeTotLPMT", &TotalPE);
    Ge68.SetBranchAddress("HitTimeCalibTOF", &HitTime);
    Ge68.SetBranchAddress("PmtIDCalib", &PmtID);


    TH1F* h = new TH1F("hTime", "hTime", nbin, xmin, xmax);

    int nEntries = Ge68.GetEntries();
    int nEvents = 0;
    int nHits = 0;

    TH1F* hRef = (TH1F*)h->Clone("hRef");
    TH1F* htemp = (TH1F*)h->Clone("htemp");
    TH1F* htempTot = (TH1F*)h->Clone("htempTot");
    TH1F* htemp_shift = (TH1F*)h->Clone("htemp_shift");
    TH1F* htemp_shiftTot = (TH1F*)h->Clone("htemp_shiftTot");

    TH1F* hCharge = new TH1F("hCharge", "hCharge", 100, 0, 3000);


    TCanvas* c1 = new TCanvas("cShift1", "Shift example", 800, 500);
    TCanvas* c2 = new TCanvas("cShift2", "Accumulated shift", 800, 500);
    
    for(size_t i = 0; i < nEntries; i++) {
        Ge68.GetEntry(i);
        hCharge->Fill(TotalPE);
    }

    float mean_q = hCharge->GetMean();
    float sig_q = hCharge->GetStdDev();
    TF1* fgaus = new TF1("fgaus", "gaus", mean_q-sig_q, mean_q+1.5*sig_q);

    mean_q = fgaus->GetParameter(1);
    sig_q = fgaus->GetParameter(2);

    hCharge->Fit("fgaus", "R0");

    TCanvas* cCharge = new TCanvas("cCharge", "cCharge", 800, 500);
    hCharge->Draw();
    fgaus->SetLineColor(kRed+1);
    fgaus->Draw("same");
    cCharge->SaveAs("Ge68MC_charge.root");

    mean_q = fgaus->GetParameter(1);
    sig_q = fgaus->GetParameter(2);

    std::cout << "Filling histogram" << std::endl;

    int nbplot = 0;
    int pass = -1;
    for (size_t i = 0; i < nEntries; i++){

        if (i % 1000 == 0 || i == nEntries - 1) {
            ProgressBar(i + 1, nEntries);
        }

        Ge68.GetEntry(i);
        if(!(TotalPE >= mean_q-sig_q && TotalPE <= mean_q+sig_q)) continue;
        
        pass++;
        
        nEvents++;
        htemp->Reset();
        for (size_t hit = 0; hit < HitTime->size(); hit++){
            if(PmtID->at(hit) > 17612) continue;
            nHits++;
            htemp->Fill(HitTime->at(hit));
        }
        
        double max = shift(htemp, htemp_shift, TargetTime); //shift histogram and return max value of this event
        h->Add(htemp_shift);
        htempTot->Add(htemp);
        htemp_shiftTot->Add(htemp_shift);

        if(pass == 0) hRef->Add(htemp_shift);

        if(nbplot == 0 && pass > 1){

            nbplot++;

            int maxy = hRef->GetMaximum()*1.2;
            int miny = maxy/3;
            gStyle->SetOptStat(0);

            c1->cd();
            hRef->SetLineColor(kCyan+2);
            hRef->Draw("HIST");
            hRef->SetTitle("Hit Time correction;Time [ns];Entries");
            htemp->Draw("HIST SAME");
            htemp_shift->SetLineColor(kMagenta+2);
            htemp_shift->Draw("HIST SAME");

            TLegend* leg = new TLegend(0.50, 0.25, 0.85, 0.35);
            leg->AddEntry(hRef, "Reference Histogram", "l");
            leg->AddEntry(htemp, "Original Hit Time", "l");
            leg->AddEntry(htemp_shift, "Shifted Hit Time", "l");
            leg->SetBorderSize(2); 
            leg->Draw();
            
            TPad* zoomPad = new TPad("zoomPad", "Zoomed View",0.40, 0.50, 0.90, 0.80);     
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
        // TH1F* htemp_clone = (TH1F*)htemp_shift->Clone(Form("htemp_clone_%d", int(i)));
        // c2->cd();
        // htemp_clone->Draw("PLC SAME");
    }

    TF1* f = new TF1("f", "[0]", 0, 200);
    h->Fit("f", "R");

    for(int bin = 1; bin <= h->GetNbinsX(); bin++){
        double val = h->GetBinContent(bin);
        double content = val - f->GetParameter(0);
        if (content < 0) content = 0;
        h->SetBinContent(bin, content);
    }

    std::cout << "Creating PDF file" << std::endl;
    // h->Scale(1.0/h->Integral());
    TH1F* hPDF = (TH1F*)h->Clone("hPDF");
    hPDF->Scale(1.0/hPDF->Integral("width"));

    TTree* Tree = new TTree("tree", "Ge68 Tree");
    Tree->Branch("nEvents", &nEvents);
    Tree->Branch("nHits", &nHits);

    Tree->Fill();

    TFile* fout = new TFile(Form("/sps/juno/mlecocq/oPs/util/Ge68%sPDF_%ibins.root", data_type.c_str(), nbin), "RECREATE");
    fout->cd();
    Tree->Write();
    h->Write();
    hPDF->Write();
    htempTot->Write();
    htemp_shiftTot->Write();
    c1->Write();
    // c2->Write();
    fout->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - begin;

    std::cout << "Info: Ge68 PDF generated in " << duration.count() << " mins" << std::endl;

    hGe68 = (TH1F*)h->Clone();
    delete h;
    delete c1;
    delete c2;

    return true;
}

bool PDFManager::CreateDataGe68PDF(const std::string& fname){
    
    std::cout << "Info: Generating Ge68 PDF" << std::endl;

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

    TCanvas* c1 = new TCanvas("cShift1", "Shift example", 800, 500);
    TCanvas* c2 = new TCanvas("cShift2", "Accumulated shift", 800, 500);

    for(size_t i = 0; i < nEntries; i++) {
        Ge68Calib.GetEntry(i);
        if(ChargeTotLpmt < 1000) continue;
        hCharge->Fill(ChargeTotLpmt);
    }

    float mean_q = hCharge->GetMean();
    float sig_q = hCharge->GetStdDev();
    TF1* fgaus = new TF1("fgaus", "gaus", mean_q-sig_q, mean_q+1.5*sig_q);

    hCharge->Fit("fgaus", "R");

    mean_q = fgaus->GetParameter(1);
    sig_q = fgaus->GetParameter(2);

    TCanvas* cCharge = new TCanvas("cCharge", "cCharge", 800, 500);
    hCharge->Draw();
    fgaus->SetLineColor(kRed+1);
    fgaus->Draw("same");
    cCharge->SaveAs("Ge68data_charge.root");

    std::cout << "Hama Only: " << HamaOnly << " suffix: " << suffix << std::endl;
    
    uint64_t tLastMuon;
    int n = -1;
    for(size_t i = 0; i < Ge68Calib.GetEntries(); i++){

        if (i % 10000 == 0 || i == nEntries - 1) {
            ProgressBar(i + 1, nEntries);
        }

        Ge68Calib.GetEntry(i);
        if(!(ChargeTotLpmt >= mean_q-sig_q && ChargeTotLpmt <= mean_q+sig_q)) continue;

        Ge68Reco.GetEntry(i);

        if(ChargeTotLpmt >= 30000 || ChargeTotWP >= 700) {
            tLastMuon = TimeStamp;
        }

        double dt = TimeStamp - tLastMuon;
        float rho = std::sqrt(RECX*RECX + RECY*RECY);

        if(dt*1e-6 > 2 && rho < 750 && std::abs(RECZ) < 500 && ChargeTotLpmt !=0 ){
            n++;
            nEvents++;
            htemp->Reset();
            for (size_t hit = 0; hit < HitTime->size(); hit++){

                if(PmtID->at(hit) > 17612) continue;
                if(HamaOnly && PmtType[PmtID->at(hit)] != "Hamamatsu") continue;

                htemp->Fill(HitTime->at(hit));
            }
        
            double max = shift(htemp, htemp_shift, TargetTime); //shift histogram and return max value of this event
            h->Add(htemp_shift);
            htempTot->Add(htemp);
            htemp_shiftTot->Add(htemp_shift);

            if(n == 0) hRef->Add(htemp_shift);
        
            if(n == 7){
                
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
    }

    std::cout << "Creating PDF file" << std::endl;
    // h->Scale(1.0/h->Integral());
    TH1F* hPDF = (TH1F*)h->Clone("hPDF");
    hPDF->Scale(1.0/hPDF->Integral("width"));

    TTree* Tree = new TTree("tree", "Ge68 Tree");
    Tree->Branch("nEvents", &nEvents);
    Tree->Branch("nHits", &nHits);

    Tree->Fill();

    TFile* fout = new TFile(Form("/sps/juno/mlecocq/oPs/util/Ge68%sPDF_%ibins%s.root", data_type.c_str(), nbin, suffix.c_str()), "RECREATE");
    fout->cd();
    Tree->Write();
    h->Write();
    hPDF->Write();
    htempTot->Write();
    htemp_shiftTot->Write();
    c1->Write();
    c2->Write();
    fout->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::ratio<60>> duration = end - begin;

    std::cout << "Info: Ge68 PDF generated in " << duration.count() << " mins" << std::endl;

    hGe68 = (TH1F*)h->Clone();
    delete h;

    return true;
}


bool PDFManager::GeneratePDF(){
    std::cout << "Info: Generating PDF with type " << data_type << std::endl;
    bool pdf = false;
    if(HamaOnly) {
        suffix = "_Hama";
        std::cout << "Info: Building PDF with Hamamatsu only. File Suffix:" << suffix << std::endl;
    }
    else suffix="";

    if(data_type == "MC") {
        fGammaName = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/gamma_1000Evts_"+juno_version+"_*_plain.root";
        // fGammaName = "/sps/juno/mlecocq/oPs/MC/plain/gamma_1000Evts_*_plain.root";
        fGe68Name = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/Ge68_1000Evts_"+juno_version+"_*_plainReco.root";
        pdf = /*CreateGammaPDF(fGammaName) &&*/ CreateSimGe68PDF(fGe68Name);
    }
    else if(data_type == "data") {
        fGammaName = "/sps/juno/mlecocq/LiHe/root/ReProd25B/LiHeTot_Times_OMILREC.root";
        fGe68Name = "/sps/juno/JunoData/Calib/2025/0824/RUN.9541.JUNODAQ.Calib-ACU-Ge68-Global-Pos-x0y0z0*_plain.root";
        pdf = /*CreateGammaPDF(fGammaName) &&*/ CreateDataGe68PDF(fGe68Name);
        
    }

    return pdf;
}

bool PDFManager::LoadPDF(){

    std::cout << "Info: Loading PDFs..." << std::endl;
    TFile* fGamma = TFile::Open(Form("/sps/juno/mlecocq/oPs/util/Gamma%sPDF.root", data_type.c_str()));
    TFile* fGe68 = TFile::Open(Form("/sps/juno/mlecocq/oPs/util/Ge68%sPDF.root", data_type.c_str()));

    if(!fGamma || !fGe68) {
        std::cout << "Error: Could not load PDFs" << std::endl;
        return false;
    }

    hGamma = (TH1F*)fGamma->Get("hTime");
    hGe68 = (TH1F*)fGe68->Get("hTime");

    return true;
}

void PDFManager::setRange(double min, double max) {
    xmin = min; 
    xmax = max;
}


double PDFManager::shift(TH1F* h, TH1F* h_shift, double x){

    h_shift->Reset();
    
    TSpectrum s(50);
    int nFound = s.Search(h, 3, "nodraw goff", 0.5); 
    double* xPeaks = s.GetPositionX();

    double bestX = h->GetBinCenter(h->GetMaximumBin()); 
    double bestHeight = -1;

    for (int i = 0; i < nFound; i++) {
        double x = xPeaks[i];
        int bin = h->FindBin(x);
        double height = h->GetBinContent(bin);

        if (height > bestHeight) {
            bestHeight = height;
            bestX = x;
        }
    }

    double delta = x - bestX;

    for(int i = 1; i <= h->GetNbinsX(); i++){
        double t = h->GetBinCenter(i);
        
        double newtime = t - delta;
        double content = 0.0;
        if(newtime >= xmin && newtime <= xmax){
            content = h->GetBinContent(h->FindBin(newtime));
        }
        h_shift->SetBinContent(i, content);
    }

    return bestX;
}

