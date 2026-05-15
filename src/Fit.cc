#include "Fit.h"
#include "PositroniumFitter.h"
#include "PDFManager.h"
#include "cxxopts.hpp"

#include "TLegend.h"
#include "TLatex.h"
#include <cstdlib>
#include <chrono>



int main(int argc, char **argv)
{
    LoadPMTInfo();

    hPrompt = new TH1F("hPrompt", "hPrompt", 200, 0, 20000);
    hDelay = new TH1F("hDelay", "hDelay", 200, 0, 6000);
    hDeltaT = new TH1F("hDeltaT", "hDeltaT", 200, 0, 1000);
    hoPsDeltaT = new TH1F("hoPsDeltaT", "hoPsDeltaT", 100, 0, 25);
    hTrueoPsDeltaT = new TH1F("hTrueoPsDeltaT", "hTrueoPsDeltaT", 100, 0, 25);
    
    suffix = "";

    // =============  Options  =============
    
    cxxopts::Options options("oPs fit algorithm", "oPs fit argument parser");
    options.add_options()
        ("start", "Starting IBD event", cxxopts::value<int>()->default_value("0"))
        ("end", "Final IBD event", cxxopts::value<int>()->default_value("0"))
        ("data-type", "Data type to fit (MC, data)", cxxopts::value<std::string>()->default_value("MC"))
        ("uniform", "Use MC with events distributed uniformly in the detector", cxxopts::value<bool>()->default_value("false"))
        ("juno-version", "Version of JUNOSW to produce MC", cxxopts::value<std::string>()->default_value("J24.1.2"))
        ("generatePDF", "Build MC/data PDF", cxxopts::value<bool>()->default_value("false"))
        ("SimEnergy", "Positron energy in MC (1, 2, 3 MeV)", cxxopts::value<std::string>()->default_value("3"))
        ("gammaFit", "Fit gamma MC instead of positronium", cxxopts::value<bool>()->default_value("false"))
        ("delay", "Fit delay signal", cxxopts::value<bool>()->default_value("false"))
        ("hama", "Select only hamamatsu PMTs", cxxopts::value<bool>()->default_value("false"))
        ("electronFit", "Fit electron MC instead of positronium", cxxopts::value<bool>()->default_value("false"))
        ("xmin", "Minimum x for histograms", cxxopts::value<double>()->default_value("0"))
        ("xmax", "Maximum x for histograms", cxxopts::value<double>()->default_value("1000"))
        ("nbin", "data binning", cxxopts::value<int>()->default_value("500"))
        ("nFit", "Number of fit attempts", cxxopts::value<int>()->default_value("3"))
        ("debug", "Debugging", cxxopts::value<bool>()->default_value("false"))
        ("reprod", "reprod version", cxxopts::value<std::string>()->default_value("ReProd25B"))
        ("h,help", "Print help information");

    auto result = options.parse(argc, argv);
    if(result.count("help")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    bool debug = result["debug"].as<bool>();
    bool genPDF = result["generatePDF"].as<bool>();
    bool delay_fit = result["delay"].as<bool>();
    HamaOnly = result["hama"].as<bool>();
    if(HamaOnly){
        std::cout << "Info: Using Hamamatsu information only!" << std::endl;
        suffix = "_Hama";
    }
    int start = result["start"].as<int>();
    int end = result["end"].as<int>();
    data_type = result["data-type"].as<std::string>();
    uniformMC = result["uniform"].as<bool>();
    juno_version = result["juno-version"].as<std::string>();
    xmin = result["xmin"].as<double>();
    xmax = result["xmax"].as<double>();
    nbin = result["nbin"].as<int>();
    SimEnergy = result["SimEnergy"].as<std::string>();
    SaveGamma = result["gammaFit"].as<bool>();
    SaveElectron = result["electronFit"].as<bool>();
    maxNFit = result["nFit"].as<int>();
    reprod = result["reprod"].as<std::string>();

    nbpIBD = 0;
    nbdIBD = 0;


    LoadPMTInfo();
    // ================================================
    // =============  Set up data and PDF  =============
    // ================================================

    PDFManager pdfManager; 
    pdfManager.setDataType(data_type);
    pdfManager.setPmtInfo(PmtType);
    pdfManager.setHamaOnly(HamaOnly);
    pdfManager.setBinning(nbin);

    if(genPDF) {
        pdfManager.setJunoVersion(juno_version);
        pdfManager.GeneratePDF();
        if(debug) return 1;
    }
    else pdfManager.LoadPDF();

    hdata = new TH1F(Form("h%s", data_type.c_str()), Form("h%s", data_type.c_str()), nbin, xmin, xmax);

    TH1F* hPDF1 = (TH1F*)pdfManager.hGamma->Clone("hGammaPDF");
    TH1F* hPDF2 = (TH1F*)pdfManager.hGe68->Clone("hGe68PDF");

    if(data_type == "MC") {
        LoadMC(SimEnergy);
        nEntries = Calib.GetEntries();
    }
    else if(data_type == "data") {
        LoadData();
        if(reprod == "ReProd25B") nEntries = IBD.GetEntries();
        else nEntries = newIBD.GetEntries();
    }

    if(!debug){
        LoadTree();
    }
    
    if (end != 0) nEntries = end;

    std::cout << "type " << data_type << std::endl;
    std::cout << "reprod " << reprod << std::endl;


    // ================================================
    // =============  Fit on data  =============
    // ================================================
    
    Params fit_results;
    PositroniumFitter Fitter(nullptr, debug, xmin, xmax, nbin, data_type);

    if(data_type == "data"){
        std::cout << "Data Entries: " << nEntries << std::endl;
        auto begin = std::chrono::high_resolution_clock::now();

        for (Long64_t i = start; i < nEntries; i++){
            
            if(reprod=="ReProd25B"){
                IBD.GetEntry(i);
                Recx = Prompt_vtx[0]; Recy = Prompt_vtx[1]; Recz = Prompt_vtx[2];
            }
            else{
                newIBD.GetEntry(i);

                if(!delay_fit && *tag!="Prompt") continue;
                else if(delay_fit && *tag!="Delay") continue;
            }
            
            if(debug){
                std::cout << "\nEvent " << i << "\n"
                << " -- Prompt energy: " << PromptE  <<std::endl;
            }

            if( PromptE < 1.022 || PromptE > 5.1) continue;
            
            // float r = std::sqrt(Prompt_vtx[0]*Prompt_vtx[0] + Prompt_vtx[1]*Prompt_vtx[1] + Prompt_vtx[2]*Prompt_vtx[2]);
            // // if(r > 14000) continue;
            

            Entry = i;
            hdata->Reset("ICES");
            hdata->GetXaxis()->SetRange(xmin, xmax);

            // =============  Fit Prompt or Delay  =============
            
            if(delay_fit && reprod == "ReProd25B"){
                for (size_t j = 0; j < DelayTOF->size(); j++){
                    if(PmtIDCalib->at(j) > 17612) continue;
                    double t = DelayTOF->at(j);
                    if(t >= xmin && t <= xmax){
                        hdata->Fill(t);
                    }
                }
            }
            else{
                for (size_t j = 0; j < PromptTOF->size(); j++){
                    if(PmtIDCalib->at(j) > 17612) continue;
                    if(HamaOnly && PmtType[PmtIDCalib->at(j)] != "Hamamatsu") continue;
                    double t = PromptTOF->at(j);
                    if(t >= xmin && t <= xmax){
                        hdata->Fill(t);
                    }
                }
            }

            RMS = hdata->GetRMS();
            std_dev = hdata->GetStdDev();
            mean = hdata->GetMean();
            skew = hdata->GetSkewness();
            kurto = hdata->GetKurtosis();

            if(std_dev > 165) continue;

            Fitter.setDataHist(hdata);
        
            for(int nFit = 0; nFit < maxNFit; nFit++){
                Fitter.DoFit();
                fit_results = Fitter.FetchParams();
                // Fitter.SetParameters(fit_results);
            }
            if(debug){
                std::string outname = "FitoPs_" + std::to_string(i) + ".root";
                Plot(fit_results, outname);
                if(i > 1) break;
            }
            else{
                energy = PromptE;
                istat = fit_results.istat;
                edm = fit_results.edm;
                chi2 = fit_results.chi2;
                chi2_ndf = fit_results.chi2ndf;
                tfirst = fit_results.t0;
                dtime = fit_results.dt;
                aIoni = fit_results.aIoni;
                aoPs = fit_results.aoPs;
                constant = fit_results.constant;
                istat = fit_results.istat;

                tfirst_err_low = fit_results.t0_err[0];
                dtime_err_low = fit_results.dt_err[0];
                aIoni_err_low = fit_results.aIoni_err[0];
                aoPs_err_low = fit_results.aoPs_err[0];
                constant_err_low = fit_results.constant_err[0];

                tfirst_err_up = fit_results.t0_err[1];
                dtime_err_up = fit_results.dt_err[1];
                aIoni_err_up = fit_results.aIoni_err[1]; 
                aoPs_err_up = fit_results.aoPs_err[1];
                constant_err_up = fit_results.constant_err[1];

                for(int j = 0; j < 4*4; j++){
                    corr[j] = fit_results.corr[j];
                }

                FitTree->Fill();
            }

        }
        
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::ratio<60>> duration = end - begin;
        std::cout << "All data processed in " << duration.count() << " mins" <<std::endl;

        // if(!debug){
        //     outfile->cd();
        //     FitTree->Write();
        //     pdfManager.hGamma->Write();
        //     pdfManager.hGe68->Write();
        //     outfile->Close();
        // }
    }


    // ================================================
    // =============  Fit on MC  =============
    // ================================================
    
    int nFit = 0;
    if(data_type == "MC"){
        std::cout << "Calib Entries " << Calib.GetEntries() << " Sim Entries " << Sim.GetEntries() << std::endl;
        auto begin = std::chrono::high_resolution_clock::now();

        for (Long64_t i = start; i < nEntries; i++){
            
            Calib.GetEntry(i);
            Reco.GetEntry(i);
            if(!SaveGamma) {
                Sim.GetEntry(i);

                Simx = SimX->at(0);
                Simy = SimY->at(0);
                Simz = SimZ->at(0);

                if(PDGID->size() < 2) continue; // no annihilation gamma
                for(int j = 0; j < PDGID->size(); j++){
                    if(PDGID->at(j) == 22) {
                        SimDecayTime = InitT->at(j);
                        break;
                    } 
                }
            }
            else SimDecayTime = 0.0;

            std::cout << "\nInfo: Executing: " << i
            << " oPs decay time " << SimDecayTime << " ns" << std::endl;

            Entry = i;
            charge = ChargeTotLPMT;

            hdata->Reset("ICES");
            hdata->GetXaxis()->SetRange(xmin, xmax);
            for (size_t j = 0; j < PromptTOF->size(); j++){
                int pmtid = PmtIDCalib->at(j);
                if(pmtid > 17612) continue;
                if(HamaOnly && PmtType[pmtid] != "Hamamatsu") continue;

                double t = PromptTOF->at(j);
                if(t >= xmin && t <= xmax) {
                    hdata->Fill(t);
                }
            }
            RMS = hdata->GetRMS();
            std_dev = hdata->GetStdDev();
            mean = hdata->GetMean();
            skew = hdata->GetSkewness();
            kurto = hdata->GetKurtosis();

            Fitter.setDataHist(hdata);
            
            for(int n = 0; n < maxNFit; n++){
                Fitter.DoFit();
                fit_results = Fitter.FetchParams();
                if(fit_results.istat >= 2) break;
                if(n < maxNFit - 1) Fitter.setParameters(fit_results);
            }
            nFit++;

            if(debug){
                std::string outname = "FitoPs_MC_" + std::to_string(i) + ".root";
                Plot(fit_results, outname, i);
                if(nFit > 2) break;
            }
            else{
                chi2 = fit_results.chi2;
                chi2_ndf = fit_results.chi2ndf;
                tfirst = fit_results.t0, dtime = fit_results.dt, aIoni = fit_results.aIoni, aoPs = fit_results.aoPs, constant = fit_results.constant;

                tfirst_err_low = fit_results.t0_err[0];
                dtime_err_low = fit_results.dt_err[0];
                aIoni_err_low = fit_results.aIoni_err[0];
                aoPs_err_low = fit_results.aoPs_err[0];
                constant_err_low = fit_results.constant_err[0];

                tfirst_err_up = fit_results.t0_err[1];
                dtime_err_up = fit_results.dt_err[1];
                aIoni_err_up = fit_results.aIoni_err[1]; 
                aoPs_err_up = fit_results.aoPs_err[1];
                constant_err_up = fit_results.constant_err[1];

                Simdt = SimDecayTime;
                istat = fit_results.istat;
                edm = fit_results.edm;
                chi2_prescan = fit_results.chi2prescan;

                for(int j = 0; j < 4*4; j++){
                    corr[j] = fit_results.corr[j];
                }

                FitTree->Fill();
            }


        }
        
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::ratio<60>> duration = end - begin;
        std::cout << "All data processed in " << duration.count() << " mins" <<std::endl;

    }
    
    if(!debug){
        outfile->cd();
        FitTree->Write();
        // hPDF1->Write();
        // hPDF2->Write();
        outfile->Close();
    }
                                                


    return 0;
    
}








void LoadData(){
    if(reprod == "ReProd25B"){
        std::string fname = "/sps/juno/mlecocq/LiHe/root/ReProd25B/LiHeTot_Times_OMILREC.root";
    
        IBD.Add(fname.c_str());
        IBD.GetListOfFiles()->Print();
    
        IBD.SetBranchStatus("*", 0);
        IBD.SetBranchStatus("E_p", 1);
        IBD.SetBranchStatus("E_d", 1);
        IBD.SetBranchStatus("vtx_p", 1);
        IBD.SetBranchStatus("HitTime_p", 1);
        IBD.SetBranchStatus("HitTime_d", 1);
    
        IBD.SetBranchAddress("E_p", &PromptE);
        IBD.SetBranchAddress("E_d", &DelayE);
        IBD.SetBranchAddress("vtx_p", &Prompt_vtx);
        IBD.SetBranchAddress("HitTime_p", &PromptTOF);
        IBD.SetBranchAddress("HitTime_d", &DelayTOF);
    }
    else{
        std::string fname = "/sps/juno/mlecocq/Data/Physics/Processing/IBD/"+reprod+"/*.root";

        float x,y,z;

        newIBD.Add(fname.c_str());
        newIBD.GetListOfFiles()->Print();

        newIBD.SetBranchStatus("*", 0);
        newIBD.SetBranchStatus("Tag", 1);
        newIBD.SetBranchStatus("RecEnergy", 1);
        newIBD.SetBranchStatus("Recx", 1);
        newIBD.SetBranchStatus("Recy", 1);
        newIBD.SetBranchStatus("Recz", 1);
        newIBD.SetBranchStatus("HitTimeCalibTOF", 1);
        newIBD.SetBranchStatus("PmtIDCalib", 1);

        newIBD.SetBranchAddress("Tag", &tag);
        newIBD.SetBranchAddress("RecEnergy", &PromptE);
        newIBD.SetBranchAddress("Recx", &Recx);
        newIBD.SetBranchAddress("Recy", &Recy);
        newIBD.SetBranchAddress("Recz", &Recz);
        newIBD.SetBranchAddress("HitTimeCalibTOF", &PromptTOF);
        newIBD.SetBranchAddress("PmtIDCalib", &PmtIDCalib);

        Prompt_vtx[0] = x; Prompt_vtx[1] = y; Prompt_vtx[2] = z;
        
    }
}






void LoadMC(std::string energy){

    std::string fname = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/positronium_1000Evts_"+energy+"MeV_"+juno_version+"*_plainReco.root";
    if(uniformMC){
        fname = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/Uniform/positronium_1000Evts_"+energy+"MeV_"+juno_version+"*_plainReco.root";
    }

    if(SaveGamma) {
        fname = "/sps/juno/mlecocq/oPs/MC/plain/gamma_1000Evts_*_plain.root";
        if(uniformMC) fname = "/sps/juno/mlecocq/oPs/MC/plain/J24.1.2/Uniform/gamma_1000Evts_*_Uniform_plain.root";
    }
    else if (SaveElectron) fname = "/sps/juno/mlecocq/Electron/plain/Electron_1000Evts_*_plainReco.root";


    Calib.Add(fname.c_str());
    Calib.SetBranchStatus("*", 0);
    Calib.SetBranchStatus("EntryNb", 1);
    Calib.SetBranchStatus("EntryNb", &CalibEvtID);
    Calib.SetBranchStatus("ChargeTotLPMT", 1);
    Calib.SetBranchStatus("PmtIDCalib", 1);
    Calib.SetBranchAddress("ChargeTotLPMT", &ChargeTotLPMT);
    Calib.SetBranchAddress("PmtIDCalib", &PmtIDCalib);
    Calib.SetBranchStatus("HitTimeCalibTOF", 1);
    Calib.SetBranchAddress("HitTimeCalibTOF", &PromptTOF);
    // Calib.SetBranchStatus("HitTimeCalibTrueTOF", 1);
    // Calib.SetBranchAddress("HitTimeCalibTrueTOF", &PromptTOF);

    Calib.GetListOfFiles()->Print();

    Reco.Add(fname.c_str());
    Reco.SetBranchStatus("*", 0);
    Reco.SetBranchStatus("Recx", 1);
    Reco.SetBranchStatus("Recy", 1);
    Reco.SetBranchStatus("Recz", 1);
    Reco.SetBranchAddress("Recx", &Recx);
    Reco.SetBranchAddress("Recy", &Recy);
    Reco.SetBranchAddress("Recz", &Recz);

    // Load true oPs dt
    if(!SaveGamma){
        fname = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/positronium_1000Evts_"+energy+"MeV_"+juno_version+"_*_plainSim.root";
        if(uniformMC){
            fname = "/sps/juno/mlecocq/oPs/MC/plain/"+juno_version+"/Uniform/positronium_1000Evts_"+energy+"MeV_"+juno_version+"_*_plainSim.root";
        }

        Sim.Add(fname.c_str());
        Sim.SetBranchStatus("*", 0);
        Sim.SetBranchStatus("EntryNb", 1);
        Sim.SetBranchStatus("PDGID", 1);
        Sim.SetBranchStatus("InitT", 1);
        Sim.SetBranchStatus("Vtx", 1);
        Sim.SetBranchStatus("Vty", 1);
        Sim.SetBranchStatus("Vtz", 1);

        Sim.SetBranchAddress("EntryNb", &SimEvtID);
        Sim.SetBranchAddress("PDGID", &PDGID);
        Sim.SetBranchAddress("InitT", &InitT);
        Sim.SetBranchAddress("Vtx", &SimX);
        Sim.SetBranchAddress("Vty", &SimY);
        Sim.SetBranchAddress("Vtz", &SimZ);
    
        for (size_t i = 0; i < Sim.GetEntries(); i++){
            Sim.GetEntry(i);
            TruedT[SimEvtID] = SimDecayTime;
            hTrueoPsDeltaT->Fill(SimDecayTime);
        }
    }

}






void LoadTree(){
    std::cout << "Info: Creating output tree" << std::endl;

    if(data_type == "data") outfile = new TFile(Form("/sps/juno/mlecocq/oPs/root/oPs_%s%s_data.root", reprod.c_str(), suffix.c_str()), "RECREATE");
    else if(data_type == "MC" && (SaveGamma || SaveElectron)) outfile = new TFile(Form("/sps/juno/mlecocq/oPs/root/SimGammaFitResult.root"), "RECREATE");
    else if (data_type == "MC" && !SaveGamma) outfile = new TFile(Form("/sps/juno/mlecocq/oPs/root/SimFitResult_%sMeV%s.root", SimEnergy.c_str(), suffix.c_str()), "RECREATE");

    FitTree = new TTree("Fit", "oPs Fit");
    FitTree->Branch("Entry", &Entry);
    FitTree->Branch("Charge", &charge);
    if(data_type == "data"){
        FitTree->Branch("Energy", &energy);
    }
    FitTree->Branch("RMS", &RMS);
    FitTree->Branch("StdDev", &std_dev);
    FitTree->Branch("Mean", &mean);
    FitTree->Branch("Skew", &skew);
    FitTree->Branch("Kurto", &kurto);
    if(data_type == "MC"){
        FitTree->Branch("Simx", &Simx);
        FitTree->Branch("Simy", &Simy);
        FitTree->Branch("Simz", &Simz);
        FitTree->Branch("Truedt", &SimDecayTime);
    }
    FitTree->Branch("Recx", &Recx);
    FitTree->Branch("Recy", &Recy);
    FitTree->Branch("Recz", &Recz);
    FitTree->Branch("t0", &tfirst);
    FitTree->Branch("dt", &dtime);
    FitTree->Branch("aIoni", &aIoni);
    FitTree->Branch("aoPs", &aoPs);
    FitTree->Branch("constant", &constant);
    FitTree->Branch("t0_err_low", &tfirst_err_low);
    FitTree->Branch("dt_err_low", &dtime_err_low);
    FitTree->Branch("aIoni_err_low", &aIoni_err_low);
    FitTree->Branch("aoPs_err_low", &aoPs_err_low);
    FitTree->Branch("constant_err_low", &constant_err_low);
    FitTree->Branch("t0_err_up", &tfirst_err_up);
    FitTree->Branch("dt_err_up", &dtime_err_up);
    FitTree->Branch("aIoni_err_up", &aIoni_err_up);
    FitTree->Branch("aoPs_err_up", &aoPs_err_up);
    FitTree->Branch("constant_err_up", &constant_err_up);
    FitTree->Branch("corr", &corr, "coor[16]/D");
    FitTree->Branch("tmin", &xmin);
    FitTree->Branch("tmax", &xmax);
    FitTree->Branch("nbin", &nbin);
    FitTree->Branch("chi2", &chi2);
    FitTree->Branch("chi2ndf", &chi2_ndf);
    FitTree->Branch("chi2prescan", &chi2_prescan);
    FitTree->Branch("istat", &istat);
    FitTree->Branch("edm", &edm);

}




void LoadPMTInfo(){
    PMT.Add("/sps/juno/mlecocq/PMTInfoLatest.root");

    Long64_t pmtid;
    char type[256];

    PMT.SetBranchStatus("*", 0);
    PMT.SetBranchStatus("CopyNo", 1);
    PMT.SetBranchAddress("CopyNo", &pmtid);
    PMT.SetBranchStatus("PMTType", 1);
    PMT.SetBranchAddress("PMTType", type);

    for (Long64_t i = 0; i < PMT.GetEntries(); i++){
        PMT.GetEntry(i);
        PmtType[pmtid] = std::string(type);
    }
}





void Plot(Params fitparams, std::string outname, int evt){
    std::cout << "Plot Fit results of evt " << evt << std::endl;

    float tmin = 180; 
    float tmax = 550;
    TLatex latex;
    TCanvas* c1 = new TCanvas("", "", 800, 500);
    TLine* line = new TLine(tmin, fitparams.constant, tmax, fitparams.constant);
    
    hdata->SetTitle("o-Ps Fit;Time [ns];Counts/bins");
    hdata->SetLineColor(kBlack);
    hdata->SetMarkerStyle(20);
    hdata->SetMarkerColor(kBlack);
    hdata->SetMarkerSize(0.5);
    hdata->Draw("E1X0");
    hdata->GetXaxis()->SetRangeUser(tmin, tmax);

    fitparams.hIoni->SetLineColor(kAzure-4); fitparams.hIoni->SetLineStyle(2); fitparams.hIoni->Draw("SAMES HIST");
    fitparams.hoPs->SetLineColor(kViolet-5); fitparams.hoPs->SetLineStyle(2); fitparams.hoPs->Draw("SAMES HIST");
    line->SetLineColor(kSpring+2); line->SetLineStyle(2); line->Draw("SAMES");
    fitparams.hFit->SetLineColor(kRed+1); fitparams.hFit->SetLineStyle(1); fitparams.hFit->Draw("SAMES");

    TLegend* leg = new TLegend(0.5, 0.3, 0.85, 0.50);
    leg->AddEntry(hdata, "Data", "ep");
    leg->AddEntry(fitparams.hIoni, "Ionisation", "l");
    leg->AddEntry(fitparams.hoPs, "o-Ps", "l");
    leg->AddEntry(line, "Background", "l");
    leg->AddEntry(fitparams.hFit, "Total Fit", "l");
    leg->SetBorderSize(2); 
    leg->Draw();

    latex.DrawLatexNDC(0.40, 0.7, Form("Fitted dt = %.2f^{+%.2f}_{-%.2f} ns", fitparams.dt, fitparams.dt_err[1], fitparams.dt_err[0]));
    if(data_type == "MC" && !SaveGamma){
        latex.DrawLatexNDC(0.40, 0.6, Form("True dt = %.2f ns", SimDecayTime));
    }
    else if (data_type == "MC" && SaveGamma){
        latex.DrawLatexNDC(0.40, 0.6, Form("True dt = %.2f ns", 0.0));
    }

    c1->SaveAs(outname.c_str());
}
