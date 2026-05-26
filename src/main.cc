#include <iostream>
#include <vector>

#include "cxxopts.hpp"
#include "Fitter.h"
#include "Helper.h"

#include "TLatex.h"
#include "TLine.h"

double SimDecayTime;
void Plot(TH1F* h, FitParameters fitparams, std::string type, int evt);

int main(int argc, char** argv){

    // ================================================
    // =============  Options  =============
    // ================================================
    
    cxxopts::Options options("oPs fit algorithm", "oPs fit argument parser");
    options.add_options()
        ("begin", "Entry of IBD event", cxxopts::value<int>()->default_value("0"))
        ("end", "Final entry of IBD event", cxxopts::value<int>()->default_value("-1"))
        ("data-type", "Data type to fit [MC, data]", cxxopts::value<std::string>()->default_value("data"))
        ("generatePDF", "Build MC/data PDF", cxxopts::value<bool>()->default_value("false"))
        ("sourcePDF", "Build MC/data PDF", cxxopts::value<std::string>()->default_value("Ge68"))
        ("SimEnergy", "Positron energy in MC [1, 2, 3] MeV", cxxopts::value<std::string>()->default_value("3"))
        ("delay", "Fit delay signal", cxxopts::value<bool>()->default_value("false"))
        ("hamamatsu", "Select only hamamatsu PMTs", cxxopts::value<bool>()->default_value("false"))
        ("xmin", "Minimum x for histograms", cxxopts::value<int>()->default_value("0"))
        ("xmax", "Maximum x for histograms", cxxopts::value<int>()->default_value("1000"))
        ("nbin", "data binning", cxxopts::value<int>()->default_value("1000"))
        ("debug", "Debugging", cxxopts::value<bool>()->default_value("false"))
        ("reprod", "reprod version", cxxopts::value<std::string>()->default_value("ReProd25C"))
        ("h,help", "Print help information");

    auto args = options.parse(argc, argv);
    if(args.count("help")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    bool debug           = args["debug"].as<bool>();
    bool delay           = args["delay"].as<bool>();
    bool buildPDF        = args["generatePDF"].as<bool>();
    bool HamaOnly        = args["hamamatsu"].as<bool>();
    
    int begin            = args["begin"].as<int>();
    int end              = args["end"].as<int>();
    int xmin             = args["xmin"].as<int>();
    int xmax             = args["xmax"].as<int>();
    int nbin             = args["nbin"].as<int>();

    std::string reprod   = args["reprod"].as<std::string>();
    std::string datatype = args["data-type"].as<std::string>();
    std::string source = args["sourcePDF"].as<std::string>();

    std::string suffix = "_All";

    // ================================================
    // =============  Setup  =============
    // ================================================
    
    Helper helper(debug, xmin, xmax, nbin, datatype);
    if(HamaOnly) {
        helper.SetHamaOnly();
        suffix = "_Hama";
    }
    if(buildPDF) helper.GeneratePDF(source);
    
    FitParameters fit_results;
    Fitter PositroniumFitter(debug, HamaOnly, xmin, xmax, nbin, datatype);
    PositroniumFitter.UsePDFSource(source);

    int entry;
    float Recx, Recy, Recz, RecEnergy;
    TFile* outfile;
    TTree* FitTree;
    if(!debug){
        std::string type = (delay) ? "_delay" : "";
        std::string outfname = "/sps/juno/mlecocq/oPs/root/" + datatype + "/subfiles/Fit" + type + "_" + std::to_string(begin) + suffix + ".root";
 
        outfile = new TFile(outfname.c_str(), "RECREATE");
        FitTree = new TTree("fit", "oPs Fit result");

        FitTree->Branch("Entry",       &entry);
        FitTree->Branch("Energy",      &helper.RecEnergy);
        FitTree->Branch("Recx",        &helper.Recx);
        FitTree->Branch("Recy",        &helper.Recy);
        FitTree->Branch("Recz",        &helper.Recz);
        FitTree->Branch("t0",          &fit_results.t0);
        FitTree->Branch("dt",          &fit_results.dt);
        if(datatype == "MC"){
            FitTree->Branch("Truedt",  &SimDecayTime);
        }
        FitTree->Branch("aIoni",       &fit_results.aIoni);
        FitTree->Branch("aoPs",        &fit_results.aoPs);
        FitTree->Branch("constant",    &fit_results.constant);
        FitTree->Branch("chi2",        &fit_results.chi2);
        FitTree->Branch("chi2ndf",     &fit_results.chi2ndf);
        FitTree->Branch("chi2prescan", &fit_results.chi2prescan);
        FitTree->Branch("edm",         &fit_results.edm);
        FitTree->Branch("istat",       &fit_results.istat);
    }

    // ================================================
    // =============  Main  =============
    // ================================================
        
    helper.LoadData(delay); // Load data TChain and variables
    
    TH1F* hdata = new TH1F(Form("h%s", datatype.c_str()), Form("h%s", datatype.c_str()), nbin, xmin, xmax);
    
    int n = -1;
    int Entries = (end != -1 && end <= helper.nEntries) ? end : helper.nEntries;

    std::cout << "\nInfo Processing entries " << begin << " to " << Entries << std::endl;

    for (size_t i = begin; i < Entries; ++i){
        n++;
        entry = i;
        if(datatype == "MC") {
            helper.Sim.GetEntry(i);
            if(helper.InitT->size() <= 1) continue;

            helper.Calib.GetEntry(i);
            helper.Reco.GetEntry(i);
            SimDecayTime = helper.oPsTrueDecay[i];
        }
        if(datatype == "data") {
            helper.IBD.GetEntry(i);
            if(*helper.Tag != "Prompt" && !delay) continue;
            // else if(delay && *helper.Tag != "Delay") continue;
            if(helper.RecEnergy < 1.022 || helper.RecEnergy > 3.5) continue;
        }
        
        hdata->Reset("ICE");
        for (size_t hit = 0; hit < helper.PmtIDCalib->size(); ++hit){
            int pmtid = helper.PmtIDCalib->at(hit);
            if (!HamaOnly && pmtid > 17612) continue;
            if (HamaOnly && (pmtid <= 17612 && helper.PmtType[pmtid] != "Hamamatsu")) continue; // keep hamamatsu and SPMTs;
            
            double t = helper.HitTimeTOF->at(hit);
            if(t >= xmin && t <= xmax) hdata->Fill(t);
        }
        
        PositroniumFitter.setDataHist(hdata);
        PositroniumFitter.DoFit();
        std::cout << " Fit results of evt " << i << " with energy = " << helper.RecEnergy << " MeV" << std::endl;
        fit_results = PositroniumFitter.FetchParams();

        
        if(debug) {
            if(delay) SimDecayTime = 0.0;
            Plot(hdata, fit_results, datatype, entry);
            if ( n >= 0) break;
        }
        else FitTree->Fill();
    }
    
    if(!debug){
        outfile->cd();
        FitTree->Write();
        // helper.hGe68->Write();
        outfile->Close();
    }
    return 0;
}




void Plot(TH1F* h, FitParameters fitparams, std::string type, int evt){
    std::cout << "Plot Fit results of evt " << evt << std::endl;

    std::string outname = "FitoPs_"+std::to_string(evt)+"_"+type+".root";

    float tmin = 180; 
    float tmax = 550;
    TLatex latex;
    TCanvas* c1 = new TCanvas("", "", 800, 500);
    TLine* line = new TLine(tmin, fitparams.constant, tmax, fitparams.constant);
    
    h->SetTitle("o-Ps Fit;Time [ns];Counts / 1.0 ns");
    h->SetLineColor(kBlack);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(kBlack);
    h->SetMarkerSize(0.5);
    h->Draw("E1X0");
    h->GetXaxis()->SetRangeUser(tmin, tmax);

    fitparams.hIoni->SetLineColor(kAzure-4); fitparams.hIoni->SetLineStyle(2); fitparams.hIoni->Draw("SAMES HIST");
    fitparams.hoPs->SetLineColor(kViolet-5); fitparams.hoPs->SetLineStyle(2); fitparams.hoPs->Draw("SAMES HIST");
    line->SetLineColor(kSpring+2); line->SetLineStyle(2); line->Draw("SAMES");
    fitparams.hFit->SetLineColor(kRed+1); fitparams.hFit->SetLineStyle(1); fitparams.hFit->Draw("SAMES");

    TLegend* leg = new TLegend(0.5, 0.3, 0.85, 0.50);
    leg->AddEntry(h, "Data", "ep");
    leg->AddEntry(fitparams.hIoni, "Ionisation", "l");
    leg->AddEntry(fitparams.hoPs, "o-Ps", "l");
    leg->AddEntry(line, "Background", "l");
    leg->AddEntry(fitparams.hFit, "Total Fit", "l");
    leg->SetBorderSize(2); 
    leg->Draw();

    latex.DrawLatexNDC(0.40, 0.7, Form("Fitted dt = %.2f^{+%.2f}_{-%.2f} ns", fitparams.dt, fitparams.dt_err[1], fitparams.dt_err[0]));
    if(type == "MC"){
        latex.DrawLatexNDC(0.40, 0.6, Form("True dt = %.2f ns", SimDecayTime));
    }

    c1->SaveAs(outname.c_str());
}