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
        ("index", "Entry of IBD event", cxxopts::value<int>()->default_value("0"))
        ("data-type", "Data type to fit [MC, data]", cxxopts::value<std::string>()->default_value("data"))
        ("generatePDF", "Build MC/data PDF", cxxopts::value<bool>()->default_value("false"))
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
    
    int index            = args["index"].as<int>();
    int xmin             = args["xmin"].as<int>();
    int xmax             = args["xmax"].as<int>();
    int nbin             = args["nbin"].as<int>();

    std::string reprod   = args["reprod"].as<std::string>();
    std::string datatype = args["data-type"].as<std::string>();

    // ================================================
    // =============  Setup  =============
    // ================================================
    
    Helper helper(debug, xmin, xmax, nbin, datatype);
    if(HamaOnly) helper.SetHamaOnly();
    if(buildPDF) helper.GeneratePDF();
    
    FitParameters fit_results;
    Fitter PositroniumFitter(debug, HamaOnly, xmin, xmax, nbin, datatype);

    int entry;
    TFile* outfile;
    TTree* FitTree;
    if(!debug){
        outfile = new TFile(Form("/sps/juno/mlecocq/oPs/util/oPs_FitResult_%s.root", datatype.c_str()), "RECREATE");
        FitTree = new TTree("fit", "oPs Fit result");

        FitTree->Branch("Entry",       &entry);
        FitTree->Branch("Energy",      &helper.RecEnergy);
        FitTree->Branch("Recx",        &helper.Recx);
        FitTree->Branch("Recy",        &helper.Recy);
        FitTree->Branch("Recz",        &helper.Recz);
        FitTree->Branch("t0",          &fit_results.t0);
        FitTree->Branch("dt",          &fit_results.dt);
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
    for (size_t i = index; i < helper.nEntries; ++i){
        n++;
        entry = i;
        if(datatype == "MC") helper.Calib.GetEntry(i);
        if(datatype == "data") {
            helper.IBD.GetEntry(i);
            if(helper.Tag != "Prompt" && !delay) continue;
            else if(delay && helper.Tag != "Delay") continue;
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
        fit_results = PositroniumFitter.FetchParams();

        SimDecayTime = helper.oPsTrueDecay[i];
        if(delay) SimDecayTime = 0.0;
        if(debug) {
            Plot(hdata, fit_results, datatype, helper.EntryNb);
            if ( n >= 0) break;
        }
        if(!debug) FitTree->Fill();
    }
    
    if(!debug){
        outfile->cd();
        FitTree->Write();
        helper.hGe68->Write();
        outfile->Close();
    }

    delete hdata;

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