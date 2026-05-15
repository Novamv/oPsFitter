#ifndef PDFMANAGER_HH
#define PDFMANAGER_HH

#include <iostream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <regex>

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLegend.h"

namespace fs = std::filesystem;


class PDFManager
{
private:
    
    int nbin;
    bool HamaOnly;
    double xmin, xmax;

    std::string suffix;

    double TargetTime; //ns;
    std::string data_type, fGammaName, fGe68Name, juno_version;
    
    std::map<int, std::string> PmtType;

    double shift(TH1F*, TH1F*, double);
    bool CreateGammaPDF(const std::string&);
    bool CreateSimGe68PDF(const std::string&);
    bool CreateDataGe68PDF(const std::string&);

public:
    TH1F* hGamma, *hGe68;

    PDFManager();
    ~PDFManager();
    
    bool GeneratePDF();
    bool LoadPDF();
    void setPmtInfo(std::map<int, std::string> Types){PmtType = Types;};
    void setHamaOnly(bool hama){HamaOnly = hama;};
    void setDataType(const std::string& type){data_type = type;};
    void setRange(double min, double max);
    void setBinning(int binning){nbin = binning;};
    void setJunoVersion(std::string version){juno_version = version;};
};



#endif
