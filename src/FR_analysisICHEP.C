// C++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TChain.h>
#include <TLatex.h>
#include <TMath.h>
#include "TRandom.h"
#include <TFeldmanCousins.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLorentzVector.h>
//#include "Math/VectorUtil.h"
//#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TFrame.h"
//#include "TGraphAsymmErrors.H"

// other includes
#include "setTDRStyle.C"

using namespace std;

// globals
// Lum
const double LUMI_INT = 58800; // /pb

bool isData=true;
TString varName;
TString Title;
TString fileName;
const int DEBUG_LEVEL = 0; // 0 - ERRORs, 1 - ERRORs & WARNINGs, 2 - ERRORs & WARNINGs & INFO
const double CUT_ELPT = 7.;
const double CUT_MUPT = 5.;
const double CUT_MZLOW = 71.;
const double CUT_MZHIGH = 111.;
const double CUT_MZ1LOW = 40.;
const double CUT_MZ1HIGH = 120.;
const double CUT_MZ2LOW = 12;
const double CUT_MZ2HIGH = 120.;
const double CUT_METHIGH = 100.;
const double CUT_M4LLOW = 100.;
const double CUT_M4LHIGH = 170.;
const double CUT_M4LLOW_FULL = 0.;
const double CUT_M4LHIGH_FULL = 1000.;

const int NLumiSections = 4;

const TString sPlotsStore = "plotsZX/";
const int SORT_EVENTS = false;
const int MODE=1;
/*const TString RootNames[]={"root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/skimZ1L_Data_Run2017-17Nov2017_noDuplicates.root",
    "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M10To50.root",
    "root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/DYJetsToLL_M50.root","root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/TTJets.root","root://cmsio5.rc.ufl.edu//store/user/t2/users/klo/Higgs/DarkZ/NTuples/ZPlusX_Early2017_v1/WZTo3LNu.root","ZZTo4L_13TeV_powheg_pythia8_RunIISpring16MiniAODv2_slimmedZX.root"};
*/
const TString RootNames[] ={"2018_noDuplicates.root", 
			    "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190128_Run2017_ZXCR-Z1LSkim/DYJetsToLL_M10To50.root",
			    "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190128_Run2017_ZXCR-Z1LSkim/DYJetsToLL_M50.root",
			    "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190128_Run2017_ZXCR-Z1LSkim/TTJets.root",
//                          "/eos/user/g/guoj/ZX/TTTo2L2Nu_TuneCP5_13TeV_2018_filter2l.root",
			    "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190128_Run2017_ZXCR-Z1LSkim/WZTo3LNu.root",
			    "/raid/raid7/lucien/Higgs/HZZ4l/NTuple/ZPlusX/ZXCR/20190128_Run2017_ZXCR-Z1LSkim/ZZTo4L_ext1.root"}; 

//ZZ sample missing for now

const TString RootNickNames[]={"Data","DY10","DY50","TT","WZ","ZZ"};
// print out and debugging settings
const int SILENT = true;
int printOutWidth = 12;
int printOutPrecision = 3;

int n[NLumiSections];
int d[NLumiSections];

int n_257969,n_258440,n_258749,n_259891;
int d_257969,d_258440,d_258749,d_259891;

int n_prompt,n_nonpromptconv, n_nonprompttaus,n_fakes,n_others,n_pass, n_nonpromptformphot,n_bdfake;

double l_257969,l_258440,l_258749,l_259891;

double lumiSections[]={433.683,414.77,401.103,297.629};

fstream CR2P2F_4e,CR2P2F_4mu,CR2P2F_2e2mu,CR3P1F_4e,CR3P1F_4mu,CR3P1F_2e2mu;
TH2D* RM2;
double Pr_2e2mu, Pr_2mu2e;
double HistSizemu,HistSizeel,HistSizemu_d,HistSizeel_d;
TH1D *h1D_m4l_2P2F,*h1D_m4l_3P1F,*h1D_m4l_4P;
TH1D *h1D_m4l_2P2F_2e2mu,*h1D_m4l_3P1F_2e2mu,*h1D_m4l_4P_2e2mu,  *h1D_m4l_2P2F_2mu2e,*h1D_m4l_3P1F_2mu2e,*h1D_m4l_4P_2mu2e;
TH1D *h1D_m4l_2P2F_4e,*h1D_m4l_3P1F_4e,*h1D_m4l_4P_4e;
TH1D *h1D_m4l_2P2F_4mu,*h1D_m4l_3P1F_4mu,*h1D_m4l_4P_4mu;
int nFailedLeptons,nFailedLeptonsZ1,nFailedLeptonsZ2;

TString processFileName;

TH1D *Hist_prompt[8],*Hist_fakes[8],*Hist_conv[8],*Hist_BDfakes[8],*CR[4][8];

TH1D  *CR_mass4lREFIT, *CR_Dbkg, *CR_massZ1, *CR_massZ2, *CR_Djet;

TH1D  *CR_shapes[5],*CR_shapes_4e[5], *CR_shapes_4mu[5],*CR_shapes_2e2mu[5];
TH1D  *SR_shapes[7],*SR_shapes_4e[7],*SR_shapes_4mu[7],*SR_shapes_2e2mu[7];
TH1D  *CSR_shapes[7],*CSR_shapes_4e[7],*CSR_shapes_4mu[7],*CSR_shapes_2e2mu[7];

TH1D *SR_mass4Refit[3],*SR_mass4Refit_4e[3],*SR_mass4Refit_4mu[3],*SR_mass4Refit_2e2mu[3];
TH2D *D2plotmass4lDbkg;

double D2plotmass4lDbkg_mass4l_nBin = 30 ;
double D2plotmass4lDbkg_Dbkg_nBin = 30 ;

double D2plotmass4lDbkg_mass4l_min = 0 ;
double D2plotmass4lDbkg_Dbkg_min = 0 ;

double D2plotmass4lDbkg_mass4l_max = 1000 ;
double D2plotmass4lDbkg_Dbkg_max = 1 ;




enum {mass4lREFIT,Dbkg, massZ1, massZ2, Djet,mass4lErr,mass4lErrREFIT};

double CR_mass4lREFIT_min = 70;
double CR_mass4lREFIT_max = 870;
double CR_mass4lREFIT_nBins =400;

double CR_Dbkg_min = 0;
double CR_Dbkg_max = 1;
double CR_Dbkg_nBins = 20;

double CR_massZ1_min = 40;
double CR_massZ1_max = 120;
double CR_massZ1_nBins = 40;

double CR_massZ2_min =12;
double CR_massZ2_max =120;
double CR_massZ2_nBins =54;

double CR_Djet_min =0;
double CR_Djet_max =1;
double CR_Djet_nBins =20;

double CR_mass4lErr_min =0;
double CR_mass4lErr_max =0.2;
double CR_mass4lErr_nBins =80;


double CR_mass4lErrREFIT_min =0;
double CR_mass4lErrREFIT_max =0.2;
double CR_mass4lErrREFIT_nBins =80;







double CR_min[] =  {CR_mass4lREFIT_min,CR_Dbkg_min, CR_massZ1_min, CR_massZ2_min,  CR_Djet_min, CR_mass4lErr_min,CR_mass4lErrREFIT_min};
double CR_max[] =  {CR_mass4lREFIT_max,CR_Dbkg_max, CR_massZ1_max, CR_massZ2_max,  CR_Djet_max, CR_mass4lErr_max,CR_mass4lErrREFIT_max};
double CR_nBins[] =  {CR_mass4lREFIT_nBins,CR_Dbkg_nBins, CR_massZ1_nBins, CR_massZ2_nBins,  CR_Djet_nBins, CR_mass4lErr_nBins,CR_mass4lErrREFIT_nBins};

TString Shape_names[] = {"mass4lREFIT","Dbkg","massZ1","massZ2","Djet","mass4lErr","mass4lErrREFIT"};



TString order[]={"EB_n","EB_d","EE_n","EE_d","MB_n","MB_d","ME_n","ME_d"};


enum {EB_n,EB_d,EE_n,EE_d,MB_n,MB_d,ME_n,ME_d};

enum {Prom,Fake,Conv,BDfake};
enum {_2P2F,_2P2F_4e,_2P2F_4mu,_2P2F_2e2mu,_3P1F,_3P1F_4e,_3P1F_4mu,_3P1F_2e2mu};
TString CRSection[]={"2P2F","2P2F4e","2P2F4mu","2P2F2e2mu","3P1F","3P1F4e","3P1F4mu","3P1F2e2mu"};

// declarations

void FR_analysis();
void analyzeZ1L(TFile* fTemplateTree,TString processFileName, double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
void FR_HistCreation(TString processFileName);
int getHistFR(TChain* tree,
              TH1D* &h1D_FRel_EB,   TH1D* &h1D_FRel_EE,   TH1D* &h1D_FRmu_EB,   TH1D* &h1D_FRmu_EE,
              TH1D* &h1D_FRel_EB_d, TH1D* &h1D_FRel_EE_d, TH1D* &h1D_FRmu_EB_d, TH1D* &h1D_FRmu_EE_d,
              double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);

void analysisInit();
void cmsPreliminary(TCanvas* &c);
void setCavasAndStyles(TString canvasName, TCanvas* &c, TString stat = "", double leftMaring = 0.15, double rightMaring = 0.05, double bottomMaring = 0.15, double topMaring = 0.05);
int setHistProperties(TH1D* &hist, Width_t lineWidth, Style_t lineStyle, Color_t lineColor, Style_t fillStyle=0, Color_t fillColor=0, TString xAxisTitle = "skip", TString yAxisTitle = "skip");
int setHistProperties2D(TH2D* &hist, TString xAxisTitle = "skip", TString yAxisTitle = "skip");
int setLegendProperties(TLegend* &leg, TString sHeader = "skip", Style_t fillStyle=0, Color_t fillColor=0);
int normaliseHist2D(TH2D* &h2D, double norm = 1.);
int fillEmptyBinsHist2D(TH2D* &h2D, double floor);
int normaliseHist1D(TH1D* &h1D, double norm = 1.);
int fillEmptyBinsHist1D(TH1D* &h1D, double floor);
void RunNumCalc(bool Tight, long int RunNum);
void estimateZX(TString processFileName,TString DataFileName, double ptElCut = CUT_ELPT, double ptMuCut = CUT_MUPT, double mZ2Cut = CUT_MZ2LOW);
int getEstimatesFromCR(TTree* tree,
                       TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE,
                       TH1D* &h1D_m4l_SR_2P2F, TH1D* &h1D_m4l_SR_2P2F_4e, TH1D* &h1D_m4l_SR_2P2F_4mu,TH1D* &h1D_m4l_SR_2P2F_2e2mu, TH1D* &h1D_m4l_SR_2P2F_2mu2e, TH1D* &h1D_m4l_SR_3P1F,TH1D* &h1D_m4l_SR_3P1F_4e,TH1D* &h1D_m4l_SR_3P1F_4mu,TH1D* &h1D_m4l_SR_3P1F_2e2mu, TH1D* &h1D_m4l_SR_3P1F_2mu2e,TH1D* &h1D_m4l_Add_2P2F, TH1D* &h1D_m4l_Add_2P2F_4e, TH1D* &h1D_m4l_Add_2P2F_4mu,TH1D* &h1D_m4l_Add_2P2F_2e2mu, TH1D* &h1D_m4l_Add_2P2F_2mu2e,
                       double ptElCut, double ptMuCut, double mZ2Cut);

double getFR(int lep_id, double lep_pt, double lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE);
double getFRIso(int lep_id, double lep_pt, double lep_eta);
void PartOrigin(vector<int> *lep_matchedR03_PdgId,vector<int> *lep_matchedR03_MomId,vector<int> *lep_matchedR03_MomMomId, int lep_Hindex[], vector<int> *lep_id, double FillValue,double FillWeight, TH1D *Hist_prompt, TH1D *Hist_fakes,TH1D *Hist_BDfakes, TH1D *Hist_conv, bool isZX);

void Latex(TString capton, TString title [], double value [], double error [] , int N);
void LatexULerr(TString capton,TString title [], double value [], double Uerror[],double Lerror[], int N);
void CRPlot();
void CRMatchPlot();
void PlotNDMatched(TString Type,TString VarP);
void PlotNDProcess(TString Type, TString VarP);

//_______________________________________________________________________________________________________________________________________________

// temporary solution
long long int lNEvents = 1;
void setNEvents(TString processFileName){
    if(processFileName=="Data_2018") lNEvents = 1;
    if(processFileName=="DY10") lNEvents = 37951928.0;
    if(processFileName=="DY50") lNEvents = 18700012.0;
    if(processFileName=="TT") lNEvents = 100907248.0;
    if(processFileName=="WZ") lNEvents = 8721088.0;
    if(processFileName=="ZZ") lNEvents = 87022400.0;
}

//

void FR_analysis(){
    
    
    CR2P2F_4e.open("CR2P2F_4e.txt");
    CR2P2F_4e << "4e \n";
    
    CR2P2F_4mu.open("CR2P2F_4mu.txt");
    CR2P2F_4mu << "4mu \n";
    
    CR2P2F_2e2mu.open("CR2P2F_2e2mu.txt");
    CR2P2F_2e2mu << "2e2mu \n";
    
    CR3P1F_4e.open("CR3P1F_4e.txt");
    CR3P1F_4e << "4e \n";
    
    CR3P1F_4mu.open("CR3P1F_4mu.txt");
    CR3P1F_4mu << "4mu \n";
    
    CR3P1F_2e2mu.open("CR3P1F_2e2mu.txt");
    CR3P1F_2e2mu << "2e2mu \n";
    
    
    
    
    // analyze Z+1L region and compute fake rates
    
    TString Variables[]={"mZ1","mEt","ptl3","nVtx","etal3","m3l","mt"};
    
    n_257969=n_258440=n_258749=n_259891=0;
    d_257969=d_258440=d_258749=d_259891=0;
    
    for(int Nsect=0; Nsect<NLumiSections; Nsect++)
    {
        n[Nsect] = 0;
        d[Nsect] = 0;
    }
    
    n_prompt=0;
    n_nonpromptconv=0;
    n_nonprompttaus=0;
    n_fakes=0;
    n_others=0;
    n_pass=0;
    n_nonpromptformphot=0;
    n_bdfake=0;
    
    // Tree file for storing the histograms
    TString fOption="RECREATE";
    Title="TEST";
    
    
    
    TEfficiency *teff = new TEfficiency();
    
    l_257969=433.683;
    l_258440=414.77;
    l_258749=401.103;
    l_259891=297.629;
    for(int i=2;i<3;i++)
    {   
        varName=Variables[i];
        //Define the histogram range and bin number
        double var_plotHigh = 200.0, var_plotLow = 0; int var_nBins = 40;
        double PtlBins[] = {5,10,20,30,40,50,80};
        
        if (varName=="mZ1"){
            var_plotHigh = 120;var_plotLow = 60;var_nBins = 20;
        }
        if (varName=="mEt"){
            var_plotHigh = 50;var_plotLow = 0;var_nBins = 10;
        }
        if (varName=="nVtx"){
            var_plotHigh = 35;var_plotLow = 0;var_nBins = 10;
        }
        if (varName=="ptl3"){
            var_plotHigh = 80;var_plotLow = 0;var_nBins = 10;
        }
        if (varName=="etal3"){
            var_plotHigh = 3;var_plotLow = -3;var_nBins = 20;
        }
        if (varName=="m3l"){
            var_plotHigh = 600;var_plotLow = 0;var_nBins = 30;
        }
        if (varName=="mt"){
            var_plotHigh = 100;var_plotLow = 0;var_nBins = 10;
        }
        
        TH1D *FR[4];
        
        double var_plotHighCR = 2000.0;double var_plotLowCR = 0.0;double var_nBinsCR = 100;
        TString ID[]={"Prom","Fake","Conv","BDfake"};
        
        
        for (int count = _2P2F; count <= _3P1F_2e2mu; count++)
        {
            
            CR[Prom][count] = new TH1D("h1D_m4l_Prom_"+CRSection[count],"h1D_m4l_Prom_"+CRSection[count],var_nBinsCR, var_plotLowCR, var_plotHighCR);
            CR[Fake][count] = new TH1D("h1D_m4l_Fake_"+CRSection[count],"h1D_m4l_Fake_"+CRSection[count],var_nBinsCR, var_plotLowCR, var_plotHighCR);
            CR[Conv][count] = new TH1D("h1D_m4l_Conv_"+CRSection[count],"h1D_m4l_Conv_"+CRSection[count],var_nBinsCR, var_plotLowCR, var_plotHighCR);
            CR[BDfake][count] = new TH1D("h1D_m4l_BDfake_"+CRSection[count],"h1D_m4l_BDfake_"+CRSection[count],var_nBinsCR, var_plotLowCR, var_plotHighCR);
        }
        for (int count = EB_n; count<ME_d+1;count++)
        {
            if (varName=="ptl3") {
                Hist_conv[count] = new TH1D("Hist_conv"+order[count] ,"Hist_conv"+order[count],6,PtlBins); Hist_conv[count]->Sumw2();
                Hist_prompt[count] = new TH1D("Hist_prompt"+order[count],"Hist_prompt"+order[count],6,PtlBins); Hist_prompt[count]->Sumw2();
                Hist_fakes[count] = new TH1D("Hist_fakes"+order[count],"Hist_fakes"+order[count],6,PtlBins); Hist_fakes[count]->Sumw2();
                Hist_BDfakes[count] = new TH1D("Hist_BDfakes"+order[count],"Hist_BDfakes"+order[count],6,PtlBins); Hist_BDfakes[count]->Sumw2();
                
            }
            
            else{
                Hist_conv[count] = new TH1D("Hist_conv"+order[count] ,"Hist_conv"+order[count],var_nBins, var_plotLow, var_plotHigh); Hist_conv[count]->Sumw2();
                Hist_prompt[count] = new TH1D("Hist_prompt"+order[count],"Hist_prompt"+order[count],var_nBins, var_plotLow, var_plotHigh); Hist_prompt[count]->Sumw2();
                Hist_fakes[count] = new TH1D("Hist_fakes"+order[count],"Hist_fakes"+order[count],var_nBins, var_plotLow, var_plotHigh); Hist_fakes[count]->Sumw2();
                Hist_BDfakes[count] = new TH1D("Hist_BDfakes"+order[count],"Hist_BDfakes"+order[count],var_nBins, var_plotLow, var_plotHigh); Hist_BDfakes[count]->Sumw2();
            }
            
        }
        
        
        
        for(int j=0;j<1;j++)
        {
            processFileName=RootNames[j];
            fileName=RootNickNames[j];
            if (j==0) {
                isData=true;
                Title="Data";
            }
            else {
                Title="MC";
                isData=false;
            }
            cout << "[sample: " << processFileName <<"]" << endl;
            
            
            TFile* fTemplateTree = new TFile("Hist_"+Title+"_"+varName+"_"+RootNickNames[j]+".root", fOption);
            setNEvents(fileName);
	    std::cout<<"Nickname = "<<fileName<<" N_entries = "<<lNEvents<<std::endl;
            analyzeZ1L(fTemplateTree, processFileName);
            
            n_257969=n_258440=n_258749=n_259891=0;
            d_257969=d_258440=d_258749=d_259891=0;
            
            
            
            for(int Nsect=0; Nsect<NLumiSections; Nsect++)
            {
                n[Nsect] = 0;
                d[Nsect] = 0;
            }
            
            
            
            
            
            fTemplateTree->Close();
            TString slimmedZXFileName = "Hist_"+Title+"_"+varName+"_"+RootNickNames[j]+".root";
            if ((varName=="ptl3")&&((Title == "Data")||(fileName == "ZZ"))){
                cout<< Title<<endl;
                
                // estimateZX("Hist_Data_ptl3_WZremoved.root",RootNames[j]);
                
                //}
              // RADI  if ((varName=="ptl3")&&(Title == "Data"))
                //{
             //RADI 
                 estimateZX("Hist_Data_ptl3_Data.root",RootNames[j]);
                
                // estimateZX("Hist_Data_ptl3_WZremoved.root",RootNames[j]);
            }
            
            n_prompt=0;
            n_nonpromptconv=0;
            n_nonprompttaus=0;
            n_fakes=0;
            n_others=0;
            n_pass=0;
            n_nonpromptformphot=0;
            n_bdfake=0;
            
        }
        
        
        
        
		TCanvas *c1;
		setCavasAndStyles("c1",c1,"");
		TLegend* leg = 0;
		
		Width_t lineWidth = 2;
		// double leg_xl = 0.50, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;
		
		double leg_xl = 0.30, leg_xr = 0.70, leg_yb = 0.72, leg_yt = 0.90;
		TH1D* h1D_dummy = new TH1D("dummy", "dummy", 20, 60, 120);
		setHistProperties(h1D_dummy,1,1,kBlack,0,0,"mZ1 [GeV]","Events/3 GeV");
		
		//For data only
		// Fake rates
		// electrons FR
		TColor *color;
		//  color->SetRGB(152,152,253);
		
		enum {EB,EE,MB,ME};
		TString smallOrder[]= {"EB","EE","MB","ME"};
		
		for (int count = EB; count <= ME; count++)
		{
		    FR[Conv] = (TH1D*) Hist_conv[2*count]->Clone();
		    FR[Conv] -> Divide(Hist_conv[2*count+1]);
		    
		    FR[Fake] = (TH1D*) Hist_fakes[2*count]->Clone();
		    FR[Fake] -> Divide(Hist_fakes[2*count+1]);
		    
		    FR[Prom] = (TH1D*) Hist_prompt[2*count]->Clone();
		    FR[Prom] -> Divide(Hist_prompt[2*count+1]);
		    
		    FR[BDfake] = (TH1D*) Hist_BDfakes[2*count]->Clone();
		    FR[BDfake] -> Divide(Hist_BDfakes[2*count+1]);
		    
		    
		    cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
		    Hist_prompt[2*count]->SetFillColor(kBlue-7); Hist_prompt[2*count]->SetLineColor(kBlack); leg->AddEntry(Hist_prompt[2*count], "Prompt","F");
		    Hist_fakes[2*count]->SetFillColor(kRed-7); Hist_fakes[2*count]->SetLineColor(kBlack);  leg->AddEntry(Hist_fakes[2*count], "Fakes","F");
		    Hist_conv[2*count]->SetFillColor(kGreen-7);Hist_conv[2*count]->SetLineColor(kBlack);leg->AddEntry(Hist_conv[2*count], "Conversion","F");
		    Hist_BDfakes[2*count]->SetFillColor(kRed-10); Hist_BDfakes[2*count]->SetLineColor(kBlack);  leg->AddEntry(Hist_BDfakes[2*count], "BD Fakes","F");
		    
		    TH1D *n = (TH1D*) Hist_conv[2*count]->Clone();
		    n->Add(Hist_prompt[2*count]);
		    n->Add(Hist_fakes[2*count]);
		    
		    TH1D *n_conv = (TH1D*) Hist_conv[2*count]->Clone();
		    TH1D *n_fakes = (TH1D*) Hist_fakes[2*count]->Clone();
		    TH1D *n_prompt = (TH1D*) Hist_prompt[2*count]->Clone();
		    TH1D *n_BDfakes = (TH1D*) Hist_BDfakes[2*count]->Clone();
		    
		    
		    
		    THStack *S_n=new THStack("n","n");
		    S_n->Add(Hist_conv[2*count]);
		    S_n->Add(Hist_prompt[2*count]);
		    S_n->Add(Hist_fakes[2*count]);
		    S_n->Add(Hist_BDfakes[2*count]);
		    //   S_n->GetXaxis()->SetTitle(varName+"[GeV]");
		    //S_n->GetYaxis()->SetTitle("Events/BinWidth");
		    
		    
		    
		    S_n->Draw("hist goff");
		    leg->Draw("goff");
		    
		    c1->SaveAs(sPlotsStore+"/ID_"+varName+order[2*count]+".pdf");
		    
		    cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
		    Hist_prompt[2*count+1]->SetFillColor(kBlue-7); Hist_prompt[2*count+1]->SetLineColor(kBlue-7); leg->AddEntry(Hist_prompt[2*count+1], "Prompt","F");
		    Hist_fakes[2*count+1]->SetFillColor(kRed-7); Hist_fakes[2*count+1]->SetLineColor(kRed-7);  leg->AddEntry(Hist_fakes[2*count+1], "Fakes","F");
		    Hist_conv[2*count+1]->SetFillColor(kGreen-7);Hist_conv[2*count+1]->SetLineColor(kGreen-7);leg->AddEntry(Hist_conv[2*count+1], "Conversion","F");
		    Hist_BDfakes[2*count+1]->SetFillColor(kRed-10); Hist_BDfakes[2*count+1]->SetLineColor(kRed-10);  leg->AddEntry(Hist_BDfakes[2*count+1], "BD Fakes","F");
		    
		    TH1D *d = (TH1D*) Hist_conv[2*count+1]->Clone();
		    d->Add(Hist_prompt[2*count+1]);
		    d->Add(Hist_fakes[2*count+1]);
		    d->Add(Hist_BDfakes[2*count+1]);
		    
		    
		    
		    TH1D *d_conv = (TH1D*) Hist_conv[2*count+1]->Clone();
		    TH1D *d_fakes = (TH1D*) Hist_fakes[2*count+1]->Clone();
		    TH1D *d_prompt = (TH1D*) Hist_prompt[2*count+1]->Clone();
		    TH1D *d_BDfakes = (TH1D*) Hist_BDfakes[2*count+1]->Clone();
		    
		    n_conv->Divide(d_conv);
		    n_prompt->Divide(d_prompt);
		    n_fakes->Divide(d_fakes);
		    n_BDfakes->Divide(d_BDfakes);
		    
		    
		    double fakes_perc = 100*(d_fakes->Integral()/(d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()));
		    double BDfakes_perc = 100*(d_BDfakes->Integral()/(d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()));
		    double conv_perc = 100*(d_conv->Integral()/(d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()));
		    double prompt_perc = 100*(d_prompt->Integral()/(d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()));
		    
		    
		    double  upperLimitFake= teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()),d_fakes->Integral() , 0.683, true);
		    double lowerLimitFake = teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()), d_fakes->Integral() , 0.683, false);
		    
		    double  upperLimitBDFake= teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()),d_BDfakes->Integral()  , 0.683, true);
		    double lowerLimitBDFake = teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()), d_BDfakes->Integral() , 0.683, false);
		    
		    
		    double  upperLimitConv= teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()),d_conv->Integral()  , 0.683, true);
		    double lowerLimitConv = teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()),d_conv->Integral() , 0.683, false);
		    
		    double  upperLimitPrompt= teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()),d_prompt->Integral()  , 0.683, true);
		    double lowerLimitPrompt = teff->ClopperPearson((d_conv->Integral() + d_fakes->Integral() + d_BDfakes->Integral()+d_prompt->Integral()),d_prompt->Integral() , 0.683, false);
		    
		    
		    cout<< smallOrder[count]+ " percentage ( Z+1 region ) : "<<endl;
		    cout<< "fakes := " <<fakes_perc<<" + "<<100*upperLimitFake-fakes_perc<<" - "<<fakes_perc - 100*lowerLimitFake<<endl;
		    cout<< "BDfakes := " <<BDfakes_perc<<" + "<<100*upperLimitBDFake-BDfakes_perc<<" - "<< BDfakes_perc - 100*lowerLimitBDFake<<endl;
		    cout<< "conv := " <<conv_perc<<" + "<<100*upperLimitConv-conv_perc<<" - "<< conv_perc - 100*lowerLimitConv<<endl;
		    cout<< "prompt := " <<prompt_perc<<" + "<<100*upperLimitPrompt-prompt_perc<<" - "<< prompt_perc - 100*lowerLimitPrompt<<endl;
		    
		    
		    
		    
		    
		    
		    THStack *S_d=new THStack("d","d");
		    S_d->Add(Hist_conv[2*count+1]);
		    S_d->Add(Hist_prompt[2*count+1]);
		    S_d->Add(Hist_fakes[2*count+1]);
		    S_d->Add(Hist_BDfakes[2*count+1]);
		    // S_d->GetXaxis()->SetTitle(varName+"[GeV]");
		    // S_d->GetYaxis()->SetTitle("Events/BinWidth");
		    
		    
		    S_d->Draw("hist goff");
		    leg->Draw("goff");
		    
		    c1->SaveAs(sPlotsStore+"/ID_"+varName+order[2*count+1]+".pdf");
		    
		    
		    
		    
		    n->Divide(d);
		    n->SetName("FR_"+varName+smallOrder[count]);
		    n->Write();
		    
		    cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
		    n->SetFillColor(kBlack);n->SetLineColor(kBlack); leg->AddEntry(n, "FR combined  "+smallOrder[count],"L");
		    
		    
		    n_fakes->SetLineColor(kRed-7);n_fakes->SetFillColor(kWhite);leg->AddEntry(n_fakes, "FR fakes","F");
		    n_BDfakes->SetLineColor(kRed-10);n_BDfakes->SetFillColor(kWhite);leg->AddEntry(n_BDfakes, "FR BD fakes","F");
		    n_prompt->SetLineColor(kBlue-7);n_prompt->SetFillColor(kWhite); leg->AddEntry(n_prompt, "FR conv","F");
		    n_conv->SetLineColor(kGreen-7);n_conv->SetFillColor(kWhite);leg->AddEntry(n_conv, "FR conv","F");
		    n_prompt->SetMinimum(0);
		    n_prompt->Draw("goff");
		    n_fakes->Draw("goff same");
		    n_BDfakes->Draw("goff same");
		    n->Draw("goff e1 same");
		    n_conv->Draw("goff same");
		    n_prompt->Draw("goff same");
		    leg->Draw("goff");
		    c1->SaveAs(sPlotsStore+"/FR_"+varName+smallOrder[count]+".pdf");
		    
		    
		    TFile* fTemplateTree1 = new TFile("Hist_ID_"+varName+smallOrder[count]+"_.root", fOption);
		    fTemplateTree1->cd();
		    Hist_prompt[2*count+1]->Write();
		    Hist_fakes[2*count+1]->Write();
		    Hist_BDfakes[2*count+1]->Write();
		    Hist_conv[2*count+1]->Write();
		    Hist_prompt[2*count]->Write();
		    Hist_fakes[2*count]->Write();
		    Hist_BDfakes[2*count]->Write();
		    Hist_conv[2*count]->Write();
		    S_n->Write();
		    S_d->Write();
		    FR[Conv]->SetName("FR_conv");
		    FR[Fake]->SetName("FR_fake");
		    FR[BDfake]->SetName("FR_BD_fake");
		    FR[Prom]->SetName("FR_prom");
		    FR[Conv]->Write();
		    FR[Fake]->Write();
		    FR[BDfake]->Write();
		    FR[Prom]->Write();
		    
		    
		    
		    
		    for (int count = _2P2F; count <= _3P1F_2e2mu; count++)
		    {
			CR[Prom][count] -> Write();
			CR[Fake][count] -> Write();
			CR[BDfake][count] -> Write();
			CR[Conv][count] -> Write();
		    }
		    
		    fTemplateTree1->Close();
		}
		
		
		
		const Int_t n = 3;
		
		//Categorization_______________________________________________________________________________________________________________________________________
		
		int numCateg = 4;
		int LatexNumLines = 3;
		double Categor_Perc[numCateg][_3P1F_2e2mu] ;
		double Upper_err_Perc[numCateg][ _3P1F_2e2mu];
		double Lower_err_Perc[numCateg][_3P1F_2e2mu];
		
		for (int count = _2P2F ; count <= _3P1F_2e2mu ; count++ )
		    for (int Identity = Prom; Identity <= BDfake ; Identity++)
		    {
			Categor_Perc[Identity][count]=0;
			Upper_err_Perc[Identity][count]=0;
			Lower_err_Perc[Identity][count]=0;
			
		    }
		
		for (int count = _2P2F ; count <= _3P1F_2e2mu ; count++ )
		{
		    
		    double Total =  CR[Conv][count]->Integral() + CR[Fake][count]->Integral() + CR[BDfake][count]->Integral();
		    
		    
		    
		    for (int Identity = Fake; Identity <= BDfake ; Identity++)
		    {
			// Identiy - 1 : because of enum {Prom, Fake, Conv, BDfake};
			Categor_Perc[Identity][count] = 100.0*CR[Identity][count]->Integral()/Total;
			
			
			double  upperLimit= teff->ClopperPearson(Total, CR[Identity][count]->Integral(), 0.683, true);
			double lowerLimit = teff->ClopperPearson(Total, CR[Identity][count]->Integral(), 0.683, false);
			
			
			Upper_err_Perc[Identity][count] = upperLimit;
			Lower_err_Perc[Identity][count] = lowerLimit;
			
		    }
		    
		    double Fake_err =100.0*(Upper_err_Perc[Fake][count] - Lower_err_Perc[Fake][count]) / 2;
		    double BDfake_err = 100.0*(Upper_err_Perc[BDfake][count] - Lower_err_Perc[BDfake][count]) / 2;
		    double Conv_err = 100.0*(Upper_err_Perc[Conv][count] - Lower_err_Perc[Conv][count]) / 2;
		    
		    double UError[] = {100.0*(Upper_err_Perc[Fake][count])-Categor_Perc[Fake][count] ,100.0*(Upper_err_Perc[BDfake][count])-Categor_Perc[BDfake][count] ,100.0*(Upper_err_Perc[Conv][count]) - Categor_Perc[Conv][count] };
		    double LError[] = {Categor_Perc[Fake][count] - 100.0*(Lower_err_Perc[Fake][count]),Categor_Perc[BDfake][count] - 100.0*(Lower_err_Perc[BDfake][count]),Categor_Perc[Conv][count] - 100.0*(Lower_err_Perc[Conv][count])};
		    
		    double Error [] = {Fake_err,BDfake_err,Conv_err};
		    TString Title_Sec [] = {"Fakes" ,"BD fakes" , "Conversions" };
		    
		    double Perc [n] = {Categor_Perc[Fake][count] , Categor_Perc[BDfake][count], Categor_Perc[Conv][count]};
		    
		    // Latex( CRSection[count]+" categorisation",Title_Sec, Perc,Error, LatexNumLines);
		    
		    LatexULerr( CRSection[count]+" categorisation",Title_Sec, Perc,UError,LError, LatexNumLines);
		}
		
		
		for (int section = EB ; section <= ME; section ++ )
		{
		    
		    double Avr_FR_MC = (Hist_fakes[2*section]->Integral() + Hist_BDfakes[2*section]->Integral() + Hist_conv[2*section]->Integral() + Hist_prompt[2*section]->Integral())/(Hist_fakes[2*section+1]->Integral() + Hist_BDfakes[2*section+1]->Integral() + Hist_conv[2*section+1]->Integral() + Hist_prompt[2*section+1]->Integral());
		    
		    double upperLimit_Avr_FR_MC =teff->ClopperPearson((Hist_fakes[2*section+1]->Integral() + Hist_BDfakes[2*section+1]->Integral() + Hist_conv[2*section+1]->Integral()+ Hist_prompt[2*section+1]->Integral()), (Hist_fakes[2*section]->Integral() + Hist_BDfakes[2*section]->Integral() + Hist_conv[2*section]->Integral()+Hist_prompt[2*section]->Integral()), 0.683, true);
		    
		    double lowerLimit_Avr_FR_MC =teff->ClopperPearson((Hist_fakes[2*section+1]->Integral() + Hist_BDfakes[2*section+1]->Integral() + Hist_conv[2*section+1]->Integral()+ Hist_prompt[2*section+1]->Integral()), (Hist_fakes[2*section]->Integral() + Hist_BDfakes[2*section]->Integral() + Hist_conv[2*section]->Integral()+Hist_prompt[2*section]->Integral()), 0.683, false);
		    
		    
		    
		    double Avr_FR_fake =Hist_fakes[2*section]->Integral()/Hist_fakes[2*section+1]->Integral();
		    double upperFake = teff->ClopperPearson(Hist_fakes[2*section+1]->Integral(),Hist_fakes[2*section]->Integral(),0.683,true);
		    double lowerFake = teff->ClopperPearson(Hist_fakes[2*section+1]->Integral(),Hist_fakes[2*section]->Integral(),0.683,false);
		    
		    
		    double Avr_FR_BDfake =Hist_BDfakes[2*section]->Integral()/Hist_BDfakes[2*section+1]->Integral();
		    double upperBDfake = teff->ClopperPearson(Hist_BDfakes[2*section+1]->Integral(),Hist_BDfakes[2*section]->Integral(),0.683,true);
		    double lowerBDfake = teff->ClopperPearson(Hist_BDfakes[2*section+1]->Integral(),Hist_BDfakes[2*section]->Integral(),0.683,false);
		    
		    
		    double Avr_FR_conv =Hist_conv[2*section]->Integral()/Hist_conv[2*section+1]->Integral();
		    double upperConv = teff->ClopperPearson(Hist_conv[2*section+1]->Integral(),Hist_conv[2*section]->Integral(),0.683,true);
		    double lowerConv = teff->ClopperPearson(Hist_conv[2*section+1]->Integral(),Hist_conv[2*section]->Integral(),0.683,false);
		    
		    double Avr_FR_prompt =Hist_prompt[2*section]->Integral()/Hist_prompt[2*section+1]->Integral();
		    double upperPrompt = teff->ClopperPearson(Hist_prompt[2*section+1]->Integral(),Hist_prompt[2*section]->Integral(),0.683,true);
		    double lowerPrompt = teff->ClopperPearson(Hist_prompt[2*section+1]->Integral(),Hist_prompt[2*section]->Integral(),0.683,false);
		    
		    
		    double Avr_FR_Cat, Avr_FR_Cat_Uperr,Avr_FR_Cat_Lower;
		    std::vector<TString> Titles= {"Average FR","Avr Fake FR","Avr BD fake FR","Avr Conv FR","Avr Prompt FR"};
		    std::vector<double>  AVRFR ={Avr_FR_MC,Avr_FR_fake,Avr_FR_BDfake,Avr_FR_conv,Avr_FR_prompt};
		    std::vector<double> AVRFR_err = {(upperLimit_Avr_FR_MC-lowerLimit_Avr_FR_MC)/2 , (upperFake - lowerFake)/2  , (upperBDfake - lowerBDfake)/2  ,(upperConv - lowerConv)/2,(upperPrompt - lowerPrompt)/2 };
		    std::vector<double> AVRFR_Uerr ={upperLimit_Avr_FR_MC-Avr_FR_MC , upperFake -Avr_FR_fake   , upperBDfake - Avr_FR_BDfake  ,upperConv - Avr_FR_conv ,upperPrompt - Avr_FR_prompt };
		    std::vector<double> AVRFR_Lerr ={Avr_FR_MC - lowerLimit_Avr_FR_MC , Avr_FR_fake - lowerFake   , Avr_FR_BDfake - lowerBDfake , Avr_FR_conv - lowerConv, Avr_FR_prompt - lowerPrompt};
		    
		    for (int count = _2P2F ; count <=  _3P1F_2e2mu ; count++ )
		    {
			
			Avr_FR_Cat =  (Categor_Perc[Fake][count]/100.0)*Avr_FR_fake + (Categor_Perc[BDfake][count]/100.0)*Avr_FR_BDfake + (Categor_Perc[Conv][count]/100.0)*Avr_FR_conv;
			
			Avr_FR_Cat_Uperr = (Categor_Perc[Fake][count]/100.0) * (upperFake - Avr_FR_fake) +(Upper_err_Perc[Fake][count]-Categor_Perc[Fake][count]/100.0) * Avr_FR_fake +
			(Categor_Perc[BDfake][count]/100.0)*(upperBDfake - Avr_FR_BDfake)  + (Upper_err_Perc[BDfake][count] - Categor_Perc[BDfake][count]/100.0) * Avr_FR_BDfake +
			(Categor_Perc[Conv][count]/100.0)*(upperConv - Avr_FR_conv) + (Upper_err_Perc[Conv][count] - Categor_Perc[Conv][count]/100.0) * Avr_FR_conv;
			
			Avr_FR_Cat_Lower = (Categor_Perc[Fake][count]/100.0) * (Avr_FR_fake - lowerFake) + (Categor_Perc[Fake][count]/100.0-Lower_err_Perc[Fake][count] )* Avr_FR_fake +
			(Categor_Perc[BDfake][count]/100.0) * ( Avr_FR_BDfake - lowerBDfake) + (Categor_Perc[BDfake][count]/100.0-Lower_err_Perc[BDfake][count]) * Avr_FR_BDfake +
			(Categor_Perc[Conv][count]/100.0) * (Avr_FR_conv - lowerConv) + (Categor_Perc[Conv][count]/100.0 - Lower_err_Perc[Conv][count]) * Avr_FR_conv;
			AVRFR.push_back(Avr_FR_Cat);
			Titles.push_back(CRSection[count]);
			AVRFR_err.push_back((Avr_FR_Cat_Uperr - Avr_FR_Cat_Lower)/2);
			AVRFR_Uerr.push_back(Avr_FR_Cat_Uperr);
			AVRFR_Lerr.push_back(Avr_FR_Cat_Lower);
			
			if ((count == _2P2F_4mu) || (count == _3P1F_4mu))
			{
			    
			    double Pr = 1  - Upper_err_Perc[Conv][count];
			    Avr_FR_Cat =  (Pr*Categor_Perc[Fake][count]/100.0)*Avr_FR_fake + (Pr*Categor_Perc[BDfake][count]/100.0)*Avr_FR_BDfake + (Categor_Perc[Conv][count]/100.0+Upper_err_Perc[Conv][count])*(upperConv);
			    
			    
			    AVRFR.push_back(Avr_FR_Cat);
			    Titles.push_back(CRSection[count]+ "+ upper err");
			    AVRFR_err.push_back((Avr_FR_Cat_Uperr - Avr_FR_Cat_Lower)/2);
			    AVRFR_Uerr.push_back(Avr_FR_Cat_Uperr);
			    AVRFR_Lerr.push_back(Avr_FR_Cat_Lower);
			    
			}
			
		    }
		    
		    double _AVRFR [20], _AVRFR_err[20], _AVRFR_Lerr[20], _AVRFR_Uerr[20];
		    TString _Titles[20];
		    std::copy(AVRFR.begin(), AVRFR.end(), _AVRFR);
		    std::copy(AVRFR_err.begin(), AVRFR_err.end(), _AVRFR_err);
		    std::copy(AVRFR_Uerr.begin(), AVRFR_Uerr.end(), _AVRFR_Uerr);
		    std::copy(AVRFR_Lerr.begin(), AVRFR_Lerr.end(), _AVRFR_Lerr);
		    std::copy(Titles.begin(), Titles.end(), _Titles);
		    
		    
		    //cout << section<<endl;
		    
		    // Latex("FR "+ smallOrder[section] , _Titles, _AVRFR,_AVRFR_err,12);
		    LatexULerr("FR "+ smallOrder[section] , _Titles, _AVRFR,_AVRFR_Uerr,_AVRFR_Lerr,15);
		    
		}
		
		
		
        /*
         double Avr_FR_MC_2P2F_4mu =  (Perc_2P2F_4mu[0]/100.0)*Avr_FR_fake + (Perc_2P2F_4mu[1]/100.0)*Avr_FR_BDfake + (Perc_2P2F_4mu[2]/100.0)*Avr_FR_conv;
         double Avr_FR_MC_2P2F_4mu_err = Fake_err_2p2f_4mu/100*Avr_FR_fake + Avr_FR_fake_err*(Perc_2P2F_4mu[0]/100.0)+
         BDfake_err_2p2f_4mu/100*Avr_FR_BDfake + Avr_FR_BDfake_err*(Perc_2P2F_4mu[1]/100.0)+
         Conv_err_2p2f_4mu/100*Avr_FR_conv + Avr_FR_conv_err*(Perc_2P2F_4mu[2]/100.0);
         
         
         double Avr_FR_MC_3P1F_4mu =  (Perc_3P1F_4mu[0]/100.0)*Avr_FR_fake + (Perc_3P1F_4mu[1]/100.0)*Avr_FR_BDfake + (Perc_3P1F_4mu[2]/100.0)*Avr_FR_conv;
         double Avr_FR_MC_3P1F_4mu_err = Fake_err_3P1F_4mu/100*Avr_FR_fake + Avr_FR_fake_err*(Perc_3P1F_4mu[0]/100.0)+
         BDfake_err_3P1F_4mu/100*Avr_FR_BDfake + Avr_FR_BDfake_err*(Perc_3P1F_4mu[1]/100.0)+
         Conv_err_3P1F_4mu/100*Avr_FR_conv + Avr_FR_conv_err*(Perc_3P1F_4mu[2]/100.0);
         
         
         
         TString Titles[] = {"Average FR","Avr Fake FR","Avr BD fake FR","Avr Conv FR","2P2F 4mu", "3P1F 4mu"};
         double  AVRFR[] ={Avr_FR_MC,Avr_FR_fake,Avr_FR_BDfake,Avr_FR_conv,Avr_FR_MC_2P2F_4mu,Avr_FR_MC_3P1F_4mu};
         double AVRFR_err[] = {Avr_FR_MC_err,Avr_FR_fake_err,Avr_FR_BDfake_err,Avr_FR_conv_err,Avr_FR_MC_2P2F_4mu_err,Avr_FR_MC_3P1F_4mu_err};
         
         Latex("Muon Barrel FR", Titles, AVRFR,AVRFR_err,6);
         */
        
    }
   
     CRPlot();
     CRMatchPlot();
     PlotNDMatched("n","ptl3");
     PlotNDMatched("d","ptl3");
     PlotNDProcess("n","ptl3");
     PlotNDProcess("d","ptl3");
     PlotNDMatched("n","mZ1");
     PlotNDMatched("d","mZ1");
     PlotNDMatched("n","mEt");
     PlotNDMatched("d","mEt");
     PlotNDProcess("n","mZ1");
     PlotNDProcess("d","mZ1");
     PlotNDProcess("n","mEt");
     PlotNDProcess("d","mEt");
    
    CR2P2F_4e.close();
    CR2P2F_4mu.close();
    CR2P2F_2e2mu.close();
    
    CR3P1F_4e.close();
    CR3P1F_4mu.close();
    CR3P1F_2e2mu.close();
    
    
    
}
void CRMatchPlot()
{
    TCanvas *c1;
    setCavasAndStyles("c1",c1,""); 
    
    //TString histName[] = {"2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","2P2F_2mu2e","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu","3P1F_2mu2e"};
    TString histName[] = {"2P2F","2P2F4e","2P2F4mu","2P2F2e2mu","3P1F","3P1F4e","3P1F4mu","3P1F2e2mu"};
    TString dataName[] = {"2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu"};
    TString Fillcolors[] = {"#AF0096","#91FFFF","#0000ff","#669965","#3CB371"};
    TString Linecolors[] = {"#AF0096","#91FFFF","#0000ff","#669965","#3CB371"};
    TString StackTitle = "Proba";
    TString ProcessNames[] = {"Data","WZ","ZZ","TT","DY50","DY10"};
    TH1D *histPlot,*dataPlot;
    
    int fileSize = 6;
    int plotSize = 8;
    TString variableN = "ptl3";
    Title = "Data";
    TFile *_file1 = new TFile("Hist_"+Title+"_"+variableN+"_"+Title+".root");
    TFile *_file = new TFile("Hist_ID_"+variableN+"EB_.root");
    for (int p = 0; p< plotSize; p++)
    {
        
        TLegend* leg = 0;
        double leg_xl = 0.45, leg_xr = 0.95, leg_yb = 0.55, leg_yt = 0.85 ;
        leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Control Region "+histName[p]);
        
        THStack *A = new THStack(StackTitle,StackTitle);
        
        
        
        
        dataPlot = (TH1D*)_file1->Get("h1D_m4l_"+dataName[p]);
        leg->AddEntry(dataPlot,Title,"LP");
        if (p<=3){
            histPlot = (TH1D*)_file->Get("h1D_m4l_Conv_"+histName[p]);
            histPlot -> SetFillColor(kGreen-7);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            
            leg->AddEntry(histPlot,"Conversions","F");
            A->Add(histPlot);
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_Prom_"+histName[p]);
            histPlot -> SetFillColor(kBlue-7);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            leg->AddEntry(histPlot,"Prompt","F");
            A->Add(histPlot);
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_Fake_"+histName[p]);
            histPlot -> SetFillColor(kRed-7);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            leg->AddEntry(histPlot,"LF Fakes","F");
            A->Add(histPlot);
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_BDfake_"+histName[p]);
            histPlot -> SetFillColor(kRed-10);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            leg->AddEntry(histPlot,"HF fakes","F");
            A->Add(histPlot);
        }
        else{
            
            
            
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_Fake_"+histName[p]);
            histPlot -> SetFillColor(kRed-7);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            
            leg->AddEntry(histPlot,"LF Fakes","F");
            A->Add(histPlot);
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_BDfake_"+histName[p]);
            histPlot -> SetFillColor(kRed-10);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            
            leg->AddEntry(histPlot,"HF fakes","F");
            A->Add(histPlot);
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_Conv_"+histName[p]);
            histPlot -> SetFillColor(kGreen-7);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            
            leg->AddEntry(histPlot,"Conversions","F");
            A->Add(histPlot);
            
            histPlot = (TH1D*)_file->Get("h1D_m4l_Prom_"+histName[p]);
            histPlot -> SetFillColor(kBlue-7);
            histPlot -> SetLineColor(kBlack);
            histPlot -> SetFillStyle(1001);
            
            
            histPlot->SetLineStyle(0);
            histPlot->SetMarkerStyle(20);
            histPlot->GetXaxis()->SetLabelFont(42);
            histPlot->GetXaxis()->SetLabelOffset(0.007);
            histPlot->GetXaxis()->SetLabelSize(0.05);
            histPlot->GetXaxis()->SetTitleSize(0.06);
            histPlot->GetXaxis()->SetTitleOffset(0.9);
            histPlot->GetXaxis()->SetTitleFont(42);
            histPlot->GetYaxis()->SetLabelFont(42);
            histPlot->GetYaxis()->SetLabelOffset(0.007);
            histPlot->GetYaxis()->SetLabelSize(0.05);
            histPlot->GetYaxis()->SetTitleSize(0.06);
            histPlot->GetYaxis()->SetTitleOffset(1.25);
            histPlot->GetYaxis()->SetTitleFont(42);
            histPlot->GetZaxis()->SetLabelFont(42);
            histPlot->GetZaxis()->SetLabelOffset(0.007);
            histPlot->GetZaxis()->SetLabelSize(0.05);
            histPlot->GetZaxis()->SetTitleSize(0.06);
            histPlot->GetZaxis()->SetTitleFont(42);
            
            leg->AddEntry(histPlot,"Prompt","F");
            A->Add(histPlot);
            
            
            
        }
        
        
        
        
        
        
        // setup environment & canvas
        
    //    TCanvas *c1 = new TCanvas("c1", "myPlots",0,67,600,600);
        c1 = new TCanvas("c1", "myPlots",0,67,600,600);
	gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        c1->Range(-102.5,-10.38415,847.5,69.4939);
        c1->SetFillColor(0);
        c1->SetBorderMode(0);
        c1->SetBorderSize(2);
        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetLeftMargin(0.15);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.05);
        c1->SetBottomMargin(0.13);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        
        //setHistProperties(dataPlot,1,1,kBlack,0,0,"m_{4l}[GeV]","Events/20 GeV");
        //dataPlot->SetMaximum(1.5*dataPlot->GetMaximum());
        dataPlot->SetFillColor(kBlack);
        dataPlot->SetLineColor(kBlack);
        dataPlot->SetFillStyle(0);
        dataPlot->SetMarkerStyle(20);
        dataPlot->GetXaxis()->SetTitle("m_{4#font[12]{l}} (GeV)");
        dataPlot->GetXaxis()->SetRange(3,40);
        dataPlot->GetXaxis()->SetLabelFont(42);
        dataPlot->GetXaxis()->SetLabelSize(0.05);
        dataPlot->GetXaxis()->SetTitleSize(0.06);
        dataPlot->GetXaxis()->SetTitleOffset(0.9);
        dataPlot->GetXaxis()->SetTitleFont(42);
        dataPlot->GetYaxis()->SetTitle("Events / 20 GeV");
        dataPlot->GetYaxis()->SetLabelFont(42);
        dataPlot->GetYaxis()->SetLabelSize(0.05);
        dataPlot->GetYaxis()->SetTitleSize(0.06);
        dataPlot->GetYaxis()->SetTitleOffset(1.25);
        dataPlot->GetYaxis()->SetTitleFont(42);
        dataPlot->GetZaxis()->SetLabelFont(42);
        dataPlot->GetZaxis()->SetLabelSize(0.035);
        dataPlot->GetZaxis()->SetTitleSize(0.035);
        dataPlot->GetZaxis()->SetTitleFont(42);
        
        dataPlot->Draw("e1 goff");
        dataPlot->GetXaxis()->SetRangeUser(50,800);
        A->Draw("hist same goff");
        dataPlot->Draw("same e1 goff");
        if (p>3)
        {
            TFile *_file = new TFile("plotsZX/estimatesZXData.root");
            histPlot = (TH1D*)_file->Get("h1D_m4l_Add_"+dataName[p-4]);
            leg->AddEntry(histPlot,"2P2F contr smooth","F");
            histPlot -> SetFillColor(kWhite);
            histPlot -> SetLineColor(kMagenta);
            histPlot -> SetLineWidth(2);
            histPlot -> SetFillStyle(4000);
            histPlot -> Smooth();
            histPlot ->Draw("hist same goff");
        }
        leg->Draw("goff");
        
        TLatex *   tex = new TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)");
        tex->SetNDC();
        tex->SetTextAlign(31);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.15,0.96,"CMS");
        tex->SetNDC();
        tex->SetTextFont(61);
        tex->SetTextSize(0.0375);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.23,0.96,"Preliminary");
        tex->SetNDC();
        tex->SetTextFont(52);
        tex->SetTextSize(0.0285);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        
        c1->SaveAs(histName[p]+"_match_auto.pdf");
        c1->SaveAs(histName[p]+"_match_auto.C");
        c1->SaveAs(histName[p]+"_match_auto.root");
        
        
    }
    
    
}
void CRPlot()
{
    TCanvas *c1;
    setCavasAndStyles("c1",c1,"");   
    
    TString histName[] = {"2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","2P2F_2mu2e","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu","3P1F_2mu2e"};
    TString Fillcolors[] = {"#cc0099","#99ccff","#996666","#669966","#669966"};
    TString Linecolors[] = {"#990066","#000099","#5f3f3f","#669966","#003300"};
    TString StackTitle = "Proba";
    TString ProcessNames[] = {"Data","WZ","ZZ","TT","DY10","DY50"};
    TString LabelNames[] = {"Data","#WZ#","#Z\\gamma^*,ZZ#","#t\\bar{t}+jets#","#Z+jets#","DY10"};
    
    TH1D *histPlot,*dataPlot;
    
    int fileSize = 6;
    int plotSize = 10;
    TString variableN = "ptl3";
    for (int p = 0; p< plotSize; p++)
    {
        std::cout<<"Ulazi u plotting"<<std::endl;     
        TLegend *leg = new TLegend(0.45,0.5,1.05,0.85,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        TLegendEntry *entry;
        if (p>4) entry=leg->AddEntry("NULL","Control Region "+ histName[p] ,"h");
        else entry=leg->AddEntry("NULL","Control Region "+ histName[p] ,"h");
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(62);
        entry=leg->AddEntry("NULL","Data","LP");
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","Z + jets","F");
        int ci;
        ci = TColor::GetColor("#669966");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        ci = TColor::GetColor("#003300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","t#bar{t} + jets","F");
        
        ci = TColor::GetColor("#996666");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        ci = TColor::GetColor("#5f3f3f");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","WZ","F");
        
        ci = TColor::GetColor("#cc0099");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        ci = TColor::GetColor("#990066");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","Z#gamma*,ZZ","F");
        
        ci = TColor::GetColor("#99ccff");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        if (p>4){  ci = TColor::GetColor("#000099");
            entry->SetLineColor(ci);
            entry->SetLineStyle(1);
            entry->SetLineWidth(1);
            entry->SetMarkerColor(1);
            entry->SetMarkerStyle(21);
            entry->SetMarkerSize(1);
            entry->SetTextFont(42);
            entry=leg->AddEntry("NULL","2P2F contribution","F");
            entry->SetFillStyle(4000);
            entry->SetLineColor(6);
            entry->SetLineStyle(1);
            entry->SetLineWidth(2);
            entry->SetMarkerColor(1);
            entry->SetMarkerStyle(21);
            entry->SetMarkerSize(1);
            entry->SetTextFont(42);}
        
        
        
        
        THStack *A = new THStack(StackTitle,StackTitle);
        
        TH1F *Proba_stack_5_stack_1 = new TH1F("Proba_stack_5_stack_1","Proba",100,0,2000);
        Proba_stack_5_stack_1->SetMinimum(-5.608576);
        Proba_stack_5_stack_1->SetMaximum(51.10072);
        Proba_stack_5_stack_1->SetDirectory(0);
        Proba_stack_5_stack_1->SetStats(0);
        Proba_stack_5_stack_1->SetLineColor(TColor::GetColor("#000099"));
        Proba_stack_5_stack_1->SetLineStyle(0);
        Proba_stack_5_stack_1->SetMarkerStyle(20);
        Proba_stack_5_stack_1->GetXaxis()->SetLabelFont(42);
        Proba_stack_5_stack_1->GetXaxis()->SetLabelOffset(0.007);
        Proba_stack_5_stack_1->GetXaxis()->SetLabelSize(0.05);
        Proba_stack_5_stack_1->GetXaxis()->SetTitleSize(0.06);
        Proba_stack_5_stack_1->GetXaxis()->SetTitleOffset(0.9);
        Proba_stack_5_stack_1->GetXaxis()->SetTitleFont(42);
        Proba_stack_5_stack_1->GetYaxis()->SetLabelFont(42);
        Proba_stack_5_stack_1->GetYaxis()->SetLabelOffset(0.007);
        Proba_stack_5_stack_1->GetYaxis()->SetLabelSize(0.05);
        Proba_stack_5_stack_1->GetYaxis()->SetTitleSize(0.06);
        Proba_stack_5_stack_1->GetYaxis()->SetTitleOffset(1.25);
        Proba_stack_5_stack_1->GetYaxis()->SetTitleFont(42);
        Proba_stack_5_stack_1->GetZaxis()->SetLabelFont(42);
        Proba_stack_5_stack_1->GetZaxis()->SetLabelOffset(0.007);
        Proba_stack_5_stack_1->GetZaxis()->SetLabelSize(0.05);
        Proba_stack_5_stack_1->GetZaxis()->SetTitleSize(0.06);
        Proba_stack_5_stack_1->GetZaxis()->SetTitleFont(42);
        A->SetHistogram(Proba_stack_5_stack_1);
        
        
        for(int i = 0;i < fileSize; i++)
        {
            
            if(i == 0) Title = "Data";
            else Title = "MC";
            TFile *_file = new TFile("Hist_"+Title+"_"+variableN+"_"+ProcessNames[i]+".root");
            
            if (i==0)
            {
                dataPlot = (TH1D*)_file->Get("h1D_m4l_"+histName[p]);
                // if (p>3) leg->AddEntry(dataPlot,ProcessNames[i]+" smooth","LP");
                //else
                //leg->AddEntry(dataPlot,ProcessNames[i],"LP");
                cout<<"Data " + histName[p]+ " := "  <<dataPlot->Integral()<<endl;
                
            }
            else
            {
                histPlot = (TH1D*)_file->Get("h1D_m4l_"+histName[p]);
                histPlot -> SetFillColor(TColor::GetColor(Fillcolors[i-1]));
                histPlot -> SetLineColor(TColor::GetColor(Linecolors[i-1]));
                histPlot -> SetFillStyle(1001);
                
                histPlot->SetLineStyle(0);
                histPlot->SetMarkerStyle(20);
                histPlot->GetXaxis()->SetLabelFont(42);
                histPlot->GetXaxis()->SetLabelOffset(0.007);
                histPlot->GetXaxis()->SetLabelSize(0.05);
                histPlot->GetXaxis()->SetTitleSize(0.06);
                histPlot->GetXaxis()->SetTitleOffset(0.9);
                histPlot->GetXaxis()->SetTitleFont(42);
                histPlot->GetYaxis()->SetLabelFont(42);
                histPlot->GetYaxis()->SetLabelOffset(0.007);
                histPlot->GetYaxis()->SetLabelSize(0.05);
                histPlot->GetYaxis()->SetTitleSize(0.06);
                histPlot->GetYaxis()->SetTitleOffset(1.25);
                histPlot->GetYaxis()->SetTitleFont(42);
                histPlot->GetZaxis()->SetLabelFont(42);
                histPlot->GetZaxis()->SetLabelOffset(0.007);
                histPlot->GetZaxis()->SetLabelSize(0.05);
                histPlot->GetZaxis()->SetTitleSize(0.06);
                histPlot->GetZaxis()->SetTitleFont(42);
                
                cout<<histName[p] + " := "<<  histPlot -> Integral()<<endl;
                // leg->AddEntry(histPlot,ProcessNames[i],"F");
                A->Add(histPlot);
            }
            
            
        }
        
        
        // setup environment & canvas
        
       // TCanvas *c1 = new TCanvas("c1", "myPlots",0,67,600,600);
        c1 = new TCanvas("c1", "myPlots",0,67,600,600);
        gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        c1->Range(-102.5,-10.38415,847.5,69.4939);
        c1->SetFillColor(0);
        c1->SetBorderMode(0);
        c1->SetBorderSize(2);
        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetLeftMargin(0.15);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.05);
        c1->SetBottomMargin(0.13);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        
        // TCanvas *c1;
        // setCavasAndStyles("c1",c1,"");
        //setHistProperties(dataPlot,1,1,TColor::GetColor("#000099"),0,0,"m_{4l}[GeV]","Events/20 GeV");
        dataPlot->SetFillColor(kBlack);
        dataPlot->SetLineColor(kBlack);
        dataPlot->SetFillStyle(0);
        dataPlot->SetMarkerStyle(20);
        dataPlot->GetXaxis()->SetTitle("m_{4#font[12]{l}} (GeV)");
        dataPlot->GetXaxis()->SetRange(3,40);
        dataPlot->GetXaxis()->SetLabelFont(42);
        dataPlot->GetXaxis()->SetLabelSize(0.05);
        dataPlot->GetXaxis()->SetTitleSize(0.06);
        dataPlot->GetXaxis()->SetTitleOffset(0.9);
        dataPlot->GetXaxis()->SetTitleFont(42);
        dataPlot->GetYaxis()->SetTitle("Events / 20 GeV");
        dataPlot->GetYaxis()->SetLabelFont(42);
        dataPlot->GetYaxis()->SetLabelSize(0.05);
        dataPlot->GetYaxis()->SetTitleSize(0.06);
        dataPlot->GetYaxis()->SetTitleOffset(1.25);
        dataPlot->GetYaxis()->SetTitleFont(42);
        dataPlot->GetZaxis()->SetLabelFont(42);
        dataPlot->GetZaxis()->SetLabelSize(0.035);
        dataPlot->GetZaxis()->SetTitleSize(0.035);
        dataPlot->GetZaxis()->SetTitleFont(42);
        
        
        
        //dataPlot->SetMaximum(1.5*dataPlot->GetMaximum());
        //if (p>3) dataPlot->SetMaximum(2*dataPlot->GetMaximum());
        // dataPlot->SetMarkerStyle(20);
        //dataPlot->SetMarkerSize(1.2);
        //gStyle->SetOptFit(1);
        //gStyle->SetOptStat(0);
        //gStyle->SetOptTitle(0);
        
        
        // if (p>3) dataPlot->Smooth();
        
        dataPlot->Draw("e1 goff");
        dataPlot->GetXaxis()->SetRangeUser(50,800);
        
        
        
        
        
        //Ovo ukljuci za MC fajlove
        A->Draw("hist same goff");
        dataPlot->Draw("same e1 goff");
        if (p>4)
        {
            TFile *_file = new TFile("plotsZX/estimatesZXData.root");
            histPlot = (TH1D*)_file->Get("h1D_m4l_Add_"+histName[p-5]);
            // leg->AddEntry(histPlot,"2P2F contr smooth","F");
            
            cout<< "2P2F contr:= "<<histPlot->Integral()<<endl;;
            
            histPlot -> SetFillColor(kWhite);
            histPlot -> SetLineColor(kMagenta);
            histPlot -> SetLineWidth(2);
            histPlot -> SetFillStyle(4000);
            histPlot -> Smooth();
            //Ukljuci za MC fajlove
            histPlot ->Draw("hist same goff");
        }
        leg->Draw("goff");
        
        TLatex *   tex = new TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)");
        tex->SetNDC();
        tex->SetTextAlign(31);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.15,0.96,"CMS");
        tex->SetNDC();
        tex->SetTextFont(61);
        tex->SetTextSize(0.0375);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.23,0.96,"Preliminary");
        tex->SetNDC();
        tex->SetTextFont(52);
        tex->SetTextSize(0.0285);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        
        c1->SaveAs(histName[p]+"_auto.pdf");
        c1->SaveAs("m"+histName[p]+"_auto.C");
        c1->SaveAs(histName[p]+"_auto.root");
        
        
    }
    
    
}
void PlotNDProcess(TString Type, TString VarP)
{
    
    
    //TString histName[] = {"2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu"};
    TString histName[] = {"el_EB","el_EE","mu_EB","mu_EE"};
    TString histNameLabel[] = {"electron barrel","electron endcap","muon barrel","muon endcap"};
    TString Fillcolors[] = {"#cc0099","#99ccff","#996666","#669966","#669966"};
    TString Linecolors[] = {"#990066","#000099","#5f3f3f","#669966","#003300"};
    TString StackTitle = "Proba";
    TString ProcessNames[] = {"Data","WZ","ZZ","TT","DY10","DY50"};
    TH1D *histPlot,*dataPlot;
    
    int fileSize = 6;
    int plotSize = 4;
    TString variableN = "ptl3";
    for (int p = 0; p< plotSize; p++)
    {double leg_xl,leg_xr,leg_yb,leg_yt;
        if (VarP == "mZ1")
        {
            
            
            leg_xl = 0.25; leg_xr = 0.75; leg_yb = 0.55; leg_yt = 0.85 ;}
        else {leg_xl = 0.58; leg_xr = 0.98; leg_yb = 0.55; leg_yt = 0.85;}
        TString h = "";
        if (Type == "n") h ="numerator";
        else h = "denominator";
        
        TLegend *leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetLineColor(1);
        leg->SetLineStyle(1);
        leg->SetLineWidth(1);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        TLegendEntry *entry;
        if (p>3) entry=leg->AddEntry("NULL",histNameLabel[p]+" "+h,"h");
        else entry=leg->AddEntry("NULL",histNameLabel[p]+" "+h,"h");
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(62);
        entry=leg->AddEntry("NULL","Data","LP");
        entry->SetLineColor(1);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(20);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","Z + jets","F");
        int ci;
        ci = TColor::GetColor("#669966");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        ci = TColor::GetColor("#003300");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","t#bar{t} + jets","F");
        
        ci = TColor::GetColor("#996666");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        ci = TColor::GetColor("#5f3f3f");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","WZ","F");
        
        ci = TColor::GetColor("#cc0099");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        ci = TColor::GetColor("#990066");
        entry->SetLineColor(ci);
        entry->SetLineStyle(1);
        entry->SetLineWidth(1);
        entry->SetMarkerColor(1);
        entry->SetMarkerStyle(21);
        entry->SetMarkerSize(1);
        entry->SetTextFont(42);
        entry=leg->AddEntry("NULL","Z#gamma*,ZZ","F");
        
        ci = TColor::GetColor("#99ccff");
        entry->SetFillColor(ci);
        entry->SetFillStyle(1001);
        
        if (p>3){  ci = TColor::GetColor("#000099");
            entry->SetLineColor(ci);
            entry->SetLineStyle(1);
            entry->SetLineWidth(1);
            entry->SetMarkerColor(1);
            entry->SetMarkerStyle(21);
            entry->SetMarkerSize(1);
            entry->SetTextFont(42);
            entry=leg->AddEntry("NULL","2P2F contribution","F");
            entry->SetFillStyle(4000);
            entry->SetLineColor(6);
            entry->SetLineStyle(1);
            entry->SetLineWidth(2);
            entry->SetMarkerColor(1);
            entry->SetMarkerStyle(21);
            entry->SetMarkerSize(1);
            entry->SetTextFont(42);}
        
        //if (VarP == "mZ1")
        //{leg_xl = 0.25; leg_xr = 0.75; leg_yb = 0.55; leg_yt = 0.85 ;}
        //else {leg_xl = 0.55; leg_xr = 0.98; leg_yb = 0.55; leg_yt = 0.85;}
        
        //leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Processes "+h+" "+histNameLabel[p]);
        
        THStack *A = new THStack(StackTitle,StackTitle);
        for(int i = 0;i < fileSize; i++)
        {
            
            
            if(i == 0) Title = "Data";
            else Title = "MC";
            TFile *_file = new TFile("Hist_"+Title+"_"+VarP+"_"+ProcessNames[i]+".root");
            
            if (i==0)
            {
                dataPlot = (TH1D*)_file->Get("Data_FR"+histName[p]+"_"+Type);
                dataPlot->SetFillColor(kBlack);
                dataPlot->SetLineColor(kBlack);
                dataPlot->SetFillStyle(0);
                dataPlot->SetMarkerStyle(20);
                dataPlot->GetXaxis()->SetRange(3,40);
                dataPlot->GetXaxis()->SetLabelFont(42);
                dataPlot->GetXaxis()->SetLabelSize(0.05);
                dataPlot->GetXaxis()->SetTitleSize(0.06);
                dataPlot->GetXaxis()->SetTitleOffset(0.9);
                dataPlot->GetXaxis()->SetTitleFont(42);
                dataPlot->GetYaxis()->SetLabelFont(42);
                dataPlot->GetYaxis()->SetLabelSize(0.05);
                dataPlot->GetYaxis()->SetTitleSize(0.06);
                dataPlot->GetYaxis()->SetTitleOffset(1.6);
                dataPlot->GetYaxis()->SetTitleFont(42);
                dataPlot->GetZaxis()->SetLabelFont(40);
                dataPlot->GetZaxis()->SetLabelSize(0.035);
                dataPlot->GetZaxis()->SetTitleSize(0.035);
                dataPlot->GetZaxis()->SetTitleFont(42);
                dataPlot->GetXaxis()->SetTitle(VarP+" [GeV]");
                
                if (VarP=="mZ1") {dataPlot->GetYaxis()->SetTitle("Events/3 GeV");  }
                if (VarP=="mZ1") dataPlot->GetXaxis()->SetTitle("m_{Z_{1}} (GeV)");
                if (VarP=="mEt") dataPlot->GetYaxis()->SetTitle("Events/5 GeV");
                if (VarP=="mEt") dataPlot->GetXaxis()->SetTitle("E_{T_{miss}} (GeV)");
                
                if (VarP=="ptl3")
                {
                    dataPlot->GetXaxis()->SetRangeUser(5,80);
                    dataPlot->GetYaxis()->SetRangeUser(0,1.1*dataPlot->GetMaximum());
                }
                
                if (VarP=="mEt")
                {
                    dataPlot->GetXaxis()->SetRangeUser(0,50);
                    dataPlot->GetYaxis()->SetRangeUser(0,1.5*dataPlot->GetMaximum());
                }
                
                if (VarP=="mZ1")
                {
                    dataPlot->GetXaxis()->SetRangeUser(60,120);
                    dataPlot->GetYaxis()->SetRangeUser(0,1.5*dataPlot->GetMaximum());
                }
                
                
                //  leg->AddEntry(dataPlot,ProcessNames[i],"LP");
            }
            else
            {
                histPlot = (TH1D*)_file->Get("Data_FR"+histName[p]+"_"+Type);
                histPlot -> SetFillColor(TColor::GetColor(Fillcolors[i-1]));
                histPlot -> SetLineColor(TColor::GetColor(Linecolors[i-1]));
                histPlot -> SetFillStyle(1001);
                
                histPlot->SetLineStyle(0);
                histPlot->SetMarkerStyle(20);
                histPlot->GetXaxis()->SetLabelFont(42);
                histPlot->GetXaxis()->SetLabelOffset(0.007);
                histPlot->GetXaxis()->SetLabelSize(0.05);
                histPlot->GetXaxis()->SetTitleSize(0.06);
                histPlot->GetXaxis()->SetTitleOffset(0.9);
                histPlot->GetXaxis()->SetTitleFont(42);
                histPlot->GetYaxis()->SetLabelFont(42);
                histPlot->GetYaxis()->SetLabelOffset(0.007);
                histPlot->GetYaxis()->SetLabelSize(0.05);
                histPlot->GetYaxis()->SetTitleSize(0.06);
                histPlot->GetYaxis()->SetTitleOffset(1.25);
                histPlot->GetYaxis()->SetTitleFont(42);
                histPlot->GetZaxis()->SetLabelFont(42);
                histPlot->GetZaxis()->SetLabelOffset(0.007);
                histPlot->GetZaxis()->SetLabelSize(0.05);
                histPlot->GetZaxis()->SetTitleSize(0.06);
                histPlot->GetZaxis()->SetTitleFont(42);
                //   leg->AddEntry(histPlot,ProcessNames[i],"F");
                A->Add(histPlot);
            }
            
            
        }
        
        // setup environment & canvas
        TCanvas *c1 = new TCanvas("c1", "myPlots",0,75,600,600);
        gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        c1->Range(-102.5,-10.38415,847.5,69.4939);
        c1->SetFillColor(0);
        c1->SetBorderMode(0);
        c1->SetBorderSize(2);
        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetLeftMargin(0.18);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.05);
        c1->SetBottomMargin(0.13);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        
        dataPlot->SetMaximum(1.5*dataPlot->GetMaximum());
        if (p>3) dataPlot->SetMaximum(2*dataPlot->GetMaximum());
        //dataPlot->SetMarkerStyle(20);
        //dataPlot->SetMarkerSize(1.2);
        dataPlot->Draw("e1 goff");
        if (VarP == "mZ1")
            dataPlot->GetXaxis()->SetRangeUser(60,120);
        else if (VarP == "mEt")
            dataPlot->GetXaxis()->SetRangeUser(0,50);
        // dataPlot->GetXaxis()->SetRangeUser(0,1.5*dataPlot->GetMaximum());
        //dataPlot->GetYaxis()->SetRangeUser(0,1.5*dataPlot->GetMaximum());
        A->Draw("hist same goff");
        dataPlot->Draw("same e1 goff");
        leg->Draw("goff");
        
        TLatex *   tex = new TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)");
        tex->SetNDC();
        tex->SetTextAlign(31);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.15,0.96,"CMS");
        tex->SetNDC();
        tex->SetTextFont(61);
        tex->SetTextSize(0.0375);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.23,0.96,"Preliminary");
        tex->SetNDC();
        tex->SetTextFont(52);
        tex->SetTextSize(0.0285);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        
        
        c1->SaveAs(VarP+histName[p]+"_"+Type+"_process_auto.pdf");
        
        
    }
    
    
}

void PlotNDMatched(TString Type,TString VarP)
{
    
    TString histNameLabel[] = {"electron endcap","electron barrel","muon endcap","muon barrel"};
    
    TString histName[] = {"EE","EB","ME","MB"};
    TString dataName[] = {"el_EE","el_EB","mu_EE","mu_EB"};
    TString StackTitle = "Proba";
    TH1D *histPlot,*dataPlot;
    
    int fileSize = 5;
    int plotSize = 4;
    TString variableN = "ptl3";
    for (int p = 0; p< plotSize; p++)
    {
        double leg_xl,leg_xr,leg_yb,leg_yt;
        TString h = "";
        if (Type == "n") h ="numerator";
        else h = "denominator";
        TLegend* leg = 0;
        if (VarP == "mZ1")
        {leg_xl = 0.25; leg_xr = 0.75; leg_yb = 0.55; leg_yt = 0.85 ;}
        else {leg_xl = 0.65; leg_xr = 1.00; leg_yb = 0.55; leg_yt = 0.85;}
        leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Matched "+h+" "+histNameLabel[p]);
        
        TFile *_file1= new TFile("Hist_Data_"+VarP+"_Data.root");
        
        dataPlot = (TH1D*)_file1->Get("Data_FR"+dataName[p]+"_"+Type);
        dataPlot->SetFillColor(kBlack);
        dataPlot->SetLineColor(kBlack);
        dataPlot->SetFillStyle(0);
        dataPlot->SetMarkerStyle(20);
        dataPlot->GetXaxis()->SetRange(3,40);
        dataPlot->GetXaxis()->SetLabelFont(42);
        dataPlot->GetXaxis()->SetLabelSize(0.05);
        dataPlot->GetXaxis()->SetTitleSize(0.06);
        dataPlot->GetXaxis()->SetTitleOffset(0.9);
        dataPlot->GetXaxis()->SetTitleFont(42);
        dataPlot->GetYaxis()->SetLabelFont(42);
        dataPlot->GetYaxis()->SetLabelSize(0.05);
        dataPlot->GetYaxis()->SetTitleSize(0.06);
        dataPlot->GetYaxis()->SetTitleOffset(1.6);
        dataPlot->GetYaxis()->SetTitleFont(42);
        dataPlot->GetZaxis()->SetLabelFont(40);
        dataPlot->GetZaxis()->SetLabelSize(0.035);
        dataPlot->GetZaxis()->SetTitleSize(0.035);
        dataPlot->GetZaxis()->SetTitleFont(42);
        dataPlot->GetXaxis()->SetTitle(VarP+" [GeV]");
        
        if (VarP=="mZ1") dataPlot->GetYaxis()->SetTitle("Events/3 GeV");
        if (VarP=="mZ1") dataPlot->GetXaxis()->SetTitle("m_{Z_{1}} (GeV)");
        if (VarP=="mEt") dataPlot->GetYaxis()->SetTitle("Events/5 GeV");
        if (VarP=="mEt") dataPlot->GetXaxis()->SetTitle("E_{T_{miss}} (GeV)");
        
        
        TString Title = "Data";
        leg->AddEntry(dataPlot,Title,"LP");
        
        TFile *_file = new TFile("Hist_ID_"+VarP+histName[p]+"_.root");
        // THStack *A = (THStack*)_file->Get(Type);
        THStack *A = new THStack(StackTitle,StackTitle);
        
        
        histPlot = (TH1D*)_file->Get("Hist_conv"+histName[p]+"_"+Type);
        histPlot -> SetFillColor(kGreen-7);
        histPlot -> SetLineColor(kBlack);
        histPlot -> SetFillStyle(1001);
        
        histPlot->SetLineStyle(0);
        histPlot->SetMarkerStyle(20);
        histPlot->GetXaxis()->SetLabelFont(42);
        histPlot->GetXaxis()->SetLabelOffset(0.007);
        histPlot->GetXaxis()->SetLabelSize(0.05);
        histPlot->GetXaxis()->SetTitleSize(0.06);
        histPlot->GetXaxis()->SetTitleOffset(0.9);
        histPlot->GetXaxis()->SetTitleFont(42);
        histPlot->GetYaxis()->SetLabelFont(42);
        histPlot->GetYaxis()->SetLabelOffset(0.007);
        histPlot->GetYaxis()->SetLabelSize(0.05);
        histPlot->GetYaxis()->SetTitleSize(0.06);
        histPlot->GetYaxis()->SetTitleOffset(1.25);
        histPlot->GetYaxis()->SetTitleFont(42);
        histPlot->GetZaxis()->SetLabelFont(42);
        histPlot->GetZaxis()->SetLabelOffset(0.007);
        histPlot->GetZaxis()->SetLabelSize(0.05);
        histPlot->GetZaxis()->SetTitleSize(0.06);
        histPlot->GetZaxis()->SetTitleFont(42);
        
        
        leg->AddEntry(histPlot,"Conversions","F");
        A->Add(histPlot);
        
        histPlot = (TH1D*)_file->Get("Hist_prompt"+histName[p]+"_"+Type);
        histPlot -> SetFillColor(kBlue-7);
        histPlot -> SetLineColor(kBlack);
        histPlot -> SetFillStyle(1001);
        histPlot->SetLineStyle(0);
        histPlot->SetMarkerStyle(20);
        histPlot->GetXaxis()->SetLabelFont(42);
        histPlot->GetXaxis()->SetLabelOffset(0.007);
        histPlot->GetXaxis()->SetLabelSize(0.05);
        histPlot->GetXaxis()->SetTitleSize(0.06);
        histPlot->GetXaxis()->SetTitleOffset(0.9);
        histPlot->GetXaxis()->SetTitleFont(42);
        histPlot->GetYaxis()->SetLabelFont(42);
        histPlot->GetYaxis()->SetLabelOffset(0.007);
        histPlot->GetYaxis()->SetLabelSize(0.05);
        histPlot->GetYaxis()->SetTitleSize(0.06);
        histPlot->GetYaxis()->SetTitleOffset(1.25);
        histPlot->GetYaxis()->SetTitleFont(42);
        histPlot->GetZaxis()->SetLabelFont(42);
        histPlot->GetZaxis()->SetLabelOffset(0.007);
        histPlot->GetZaxis()->SetLabelSize(0.05);
        histPlot->GetZaxis()->SetTitleSize(0.06);
        histPlot->GetZaxis()->SetTitleFont(42);
        
        leg->AddEntry(histPlot,"Prompt","F");
        A->Add(histPlot);
        
        histPlot = (TH1D*)_file->Get("Hist_fakes"+histName[p]+"_"+Type);
        histPlot -> SetFillColor(kRed-7);
        histPlot -> SetLineColor(kBlack);
        histPlot -> SetFillStyle(1001);
        histPlot->SetLineStyle(0);
        histPlot->SetMarkerStyle(20);
        histPlot->GetXaxis()->SetLabelFont(42);
        histPlot->GetXaxis()->SetLabelOffset(0.007);
        histPlot->GetXaxis()->SetLabelSize(0.05);
        histPlot->GetXaxis()->SetTitleSize(0.06);
        histPlot->GetXaxis()->SetTitleOffset(0.9);
        histPlot->GetXaxis()->SetTitleFont(42);
        histPlot->GetYaxis()->SetLabelFont(42);
        histPlot->GetYaxis()->SetLabelOffset(0.007);
        histPlot->GetYaxis()->SetLabelSize(0.05);
        histPlot->GetYaxis()->SetTitleSize(0.06);
        histPlot->GetYaxis()->SetTitleOffset(1.25);
        histPlot->GetYaxis()->SetTitleFont(42);
        histPlot->GetZaxis()->SetLabelFont(42);
        histPlot->GetZaxis()->SetLabelOffset(0.007);
        histPlot->GetZaxis()->SetLabelSize(0.05);
        histPlot->GetZaxis()->SetTitleSize(0.06);
        histPlot->GetZaxis()->SetTitleFont(42);
        
        leg->AddEntry(histPlot,"LF Fakes","F");
        A->Add(histPlot);
        
        histPlot = (TH1D*)_file->Get("Hist_BDfakes"+histName[p]+"_"+Type);
        histPlot -> SetFillColor(kRed-10);
        histPlot -> SetLineColor(kBlack);
        histPlot -> SetFillStyle(1001);
        histPlot->SetLineStyle(0);
        histPlot->SetMarkerStyle(20);
        histPlot->GetXaxis()->SetLabelFont(42);
        histPlot->GetXaxis()->SetLabelOffset(0.007);
        histPlot->GetXaxis()->SetLabelSize(0.05);
        histPlot->GetXaxis()->SetTitleSize(0.06);
        histPlot->GetXaxis()->SetTitleOffset(0.9);
        histPlot->GetXaxis()->SetTitleFont(42);
        histPlot->GetYaxis()->SetLabelFont(42);
        histPlot->GetYaxis()->SetLabelOffset(0.007);
        histPlot->GetYaxis()->SetLabelSize(0.05);
        histPlot->GetYaxis()->SetTitleSize(0.06);
        histPlot->GetYaxis()->SetTitleOffset(1.25);
        histPlot->GetYaxis()->SetTitleFont(42);
        histPlot->GetZaxis()->SetLabelFont(42);
        histPlot->GetZaxis()->SetLabelOffset(0.007);
        histPlot->GetZaxis()->SetLabelSize(0.05);
        histPlot->GetZaxis()->SetTitleSize(0.06);
        histPlot->GetZaxis()->SetTitleFont(42);
        
        leg->AddEntry(histPlot,"HF fakes","F");
        A->Add(histPlot);
        
        // setup environment & canvas
        TCanvas *c1 = new TCanvas("c1", "myPlots",0,75,600,600);
        gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        c1->Range(-102.5,-10.38415,847.5,69.4939);
        c1->SetFillColor(0);
        c1->SetBorderMode(0);
        c1->SetBorderSize(2);
        c1->SetTickx(1);
        c1->SetTicky(1);
        c1->SetLeftMargin(0.18);
        c1->SetRightMargin(0.05);
        c1->SetTopMargin(0.05);
        c1->SetBottomMargin(0.13);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        c1->SetFrameFillStyle(0);
        c1->SetFrameBorderMode(0);
        dataPlot->SetMaximum(1.5*dataPlot->GetMaximum());
        if (p>3) dataPlot->SetMaximum(2*dataPlot->GetMaximum());
        
        // dataPlot->SetMaximum(1.5*dataPlot->GetMaximum());
        //dataPlot->SetMarkerStyle(20);
        // dataPlot->SetMarkerSize(1.2);
        if (VarP == "mZ1")
            dataPlot->GetXaxis()->SetRangeUser(60,120);
        else if (VarP == "mEt")
            dataPlot->GetXaxis()->SetRangeUser(0,50);
        //  else
        dataPlot->Draw("e1 goff");
        // if (VarP == "ptl3")
        //   dataPlot->GetXaxis()->SetRangeUser(0,80);
        //  dataPlot->GetYaxis()->SetRangeUser(0,1.5*dataPlot->GetMaximum());
        
        //dataPlot->GetXaxis()->SetRangeUser(0,1.5*dataPlot->GetMaximum());
        A->Draw("hist same goff");
        dataPlot->Draw("same e1 goff");
        leg->Draw("goff");
        
        
        TLatex *   tex = new TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)");
        tex->SetNDC();
        tex->SetTextAlign(31);
        tex->SetTextFont(42);
        tex->SetTextSize(0.03);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.15,0.96,"CMS");
        tex->SetNDC();
        tex->SetTextFont(61);
        tex->SetTextSize(0.0375);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        tex = new TLatex(0.23,0.96,"Preliminary");
        tex->SetNDC();
        tex->SetTextFont(52);
        tex->SetTextSize(0.0285);
        tex->SetLineWidth(2);
        tex->Draw("goff");
        
        
        c1->SaveAs(VarP+"_"+Type+"_"+histName[p]+"_match_auto.pdf");
        
        
    }
    
}


void RunNumCalc(bool Tight, long int RunNum)
{
    if (RunNum<=257969) {
        if (Tight) n_257969++;
        else d_257969++;
    }
    if ((RunNum<=258440)&&(RunNum>257969)) {
        if (Tight) n_258440++;
        else d_258440++;
    }
    if ((RunNum<=258749)&& (RunNum>258440)) {
        if (Tight) n_258749++;
        else d_258749++;
    }
    if ((RunNum<=259891)&& (RunNum>258749)) {
        if (Tight) n_259891++;
        else d_259891++;
    }
    
    
    
    
}

void analyzeZ1L(TFile* fTemplateTree,TString processFileName, double ptElCut, double ptMuCut, double mZ2Cut){
    TChain* zxTree = new TChain("passedEvents");
    if (!isData)
        zxTree = new TChain("Ana/passedEvents");
    
    
    
    zxTree->Add(processFileName);
    
    // common properties
    Width_t lineWidth = 2;
    double leg_xl = 0.50, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;
    
    // define dummy histogram for FRs
    
    //Define the histogram range and bin number
    double var_plotHigh = 2000.0; double var_plotLow = 0.0; double var_nBins = 100;
    
    if (varName=="mZ1"){
        var_plotHigh = 120;var_plotLow = 60;var_nBins = 20;
    }
    if (varName=="mEt"){
        var_plotHigh = 50;var_plotLow = 0;var_nBins = 10;
    }
    if (varName=="nVtx"){
        var_plotHigh = 35;var_plotLow = 0;var_nBins = 10;
    }
    if (varName=="ptl3"){
        var_plotHigh = 80;var_plotLow = 0;var_nBins = 10;
    }
    if (varName=="etal3"){
        var_plotHigh = 3;var_plotLow = -3;var_nBins = 20;
    }
    if (varName=="m3l"){
        var_plotHigh = 600;var_plotLow = 0;var_nBins = 30;
    }
    if (varName=="mt"){
        var_plotHigh = 100;var_plotLow = 0;var_nBins = 10;
    }
    double PtlBins[] = {5,10,20,30,40,50,80};
    double PtlBinsMu[] = {5,7,10,20,30,40,50,80};
    TH1D* h1D_dummy = new TH1D("dummy", "dummy", var_nBins, var_plotLow, var_plotHigh);
    
    
    
    
    TString varAxLabel = varName+" [GeV]";
    if ((varName=="nVtx")||(varName=="etal3")) {
        varAxLabel=varName;
    }
    
    
    double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    TString sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    TString sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummy,1,1,kBlack,0,0,varAxLabel,"Misidentification rate");
    
    // define FR numerator histograms
    TH1D* h1D_FRel_EB   = new TH1D("h1D_FRel_EB","h1D_FRel_EB",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EB->Sumw2();
    TH1D* h1D_FRel_EE   = new TH1D("h1D_FRel_EE","h1D_FRel_EE",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EE->Sumw2();
    TH1D* h1D_FRmu_EB   = new TH1D("h1D_FRmu_EB","h1D_FRmu_EB",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EB->Sumw2();
    TH1D* h1D_FRmu_EE   = new TH1D("h1D_FRmu_EE","h1D_FRmu_EE",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EE->Sumw2();
    
    // define FR denominator histograms
    TH1D* h1D_FRel_EB_d   = new TH1D("h1D_FRel_EB_d","h1D_FRel_EB_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EB_d->Sumw2();
    TH1D* h1D_FRel_EE_d   = new TH1D("h1D_FRel_EE_d","h1D_FRel_EE_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRel_EE_d->Sumw2();
    TH1D* h1D_FRmu_EB_d   = new TH1D("h1D_FRmu_EB_d","h1D_FRmu_EB_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EB_d->Sumw2();
    TH1D* h1D_FRmu_EE_d   = new TH1D("h1D_FRmu_EE_d","h1D_FRmu_EE_d",var_nBins, var_plotLow, var_plotHigh); h1D_FRmu_EE_d->Sumw2();
    
    
    
    if (varName=="ptl3")
    {
        TH1D* h1D_dummy = new TH1D("dummy", "dummy",6, PtlBins);
        
        
        
        
        
        
        h1D_FRel_EB   = new TH1D("h1D_FRel_EB","h1D_FRel_EB",6, PtlBins); h1D_FRel_EB->Sumw2();
        h1D_FRel_EE   = new TH1D("h1D_FRel_EE","h1D_FRel_EE",6, PtlBins); h1D_FRel_EE->Sumw2();
        
        h1D_FRmu_EB   = new TH1D("h1D_FRmu_EB","h1D_FRmu_EB",7, PtlBinsMu); h1D_FRmu_EB->Sumw2();
        h1D_FRmu_EE   = new TH1D("h1D_FRmu_EE","h1D_FRmu_EE",7, PtlBinsMu); h1D_FRmu_EE->Sumw2();
        
        // define FR denominator histograms
        h1D_FRel_EB_d   = new TH1D("h1D_FRel_EB_d","h1D_FRel_EB_d",6, PtlBins); h1D_FRel_EB_d->Sumw2();
        h1D_FRel_EE_d   = new TH1D("h1D_FRel_EE_d","h1D_FRel_EE_d",6, PtlBins); h1D_FRel_EE_d->Sumw2();
        
        
        h1D_FRmu_EB_d   = new TH1D("h1D_FRmu_EB_d","h1D_FRmu_EB_d",7, PtlBinsMu); h1D_FRmu_EB_d->Sumw2();
        h1D_FRmu_EE_d   = new TH1D("h1D_FRmu_EE_d","h1D_FRmu_EE_d",7, PtlBinsMu); h1D_FRmu_EE_d->Sumw2();
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    var_plotHigh = 2000.0; var_plotLow = 0.0; var_nBins = 100;
    
    // define CR histograms
    h1D_m4l_2P2F = new TH1D("h1D_m4l_2P2F","h1D_m4l_2P2F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_2P2F->Sumw2();
    h1D_m4l_3P1F = new TH1D("h1D_m4l_3P1F","h1D_m4l_3P1F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_3P1F->Sumw2();
    h1D_m4l_4P   = new TH1D("h1D_m4l_4P","h1D_m4l_4P",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_4P->Sumw2();
    
    // define CR histograms
    h1D_m4l_2P2F_2e2mu = new TH1D("h1D_m4l_2P2F_2e2mu","h1D_m4l_2P2F_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_2P2F_2e2mu->Sumw2();
    h1D_m4l_3P1F_2e2mu = new TH1D("h1D_m4l_3P1F_2e2mu","h1D_m4l_3P1F_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_3P1F_2e2mu->Sumw2();
    h1D_m4l_4P_2e2mu   = new TH1D("h1D_m4l_4P_2e2mu","h1D_m4l_4P_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_4P_2e2mu->Sumw2();
    
    // define CR histograms
    h1D_m4l_2P2F_2mu2e = new TH1D("h1D_m4l_2P2F_2mu2e","h1D_m4l_2P2F_2mu2e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_2P2F_2mu2e->Sumw2();
    h1D_m4l_3P1F_2mu2e = new TH1D("h1D_m4l_3P1F_2mu2e","h1D_m4l_3P1F_2mu2e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_3P1F_2mu2e->Sumw2();
    h1D_m4l_4P_2mu2e   = new TH1D("h1D_m4l_4P_2mu2e","h1D_m4l_4P_2mu2e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_4P_2mu2e->Sumw2();
    
    // define CR histograms
    h1D_m4l_2P2F_4e = new TH1D("h1D_m4l_2P2F_4e","h1D_m4l_2P2F_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_2P2F_4e->Sumw2();
    h1D_m4l_3P1F_4e = new TH1D("h1D_m4l_3P1F_4e","h1D_m4l_3P1F_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_3P1F_4e->Sumw2();
    h1D_m4l_4P_4e   = new TH1D("h1D_m4l_4P_4e","h1D_m4l_4P_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_4P_4e->Sumw2();
    
    // define CR histograms
    h1D_m4l_2P2F_4mu = new TH1D("h1D_m4l_2P2F_4mu","h1D_m4l_2P2F_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_2P2F_4mu->Sumw2();
    h1D_m4l_3P1F_4mu = new TH1D("h1D_m4l_3P1F_4mu","h1D_m4l_3P1F_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_3P1F_4mu->Sumw2();
    h1D_m4l_4P_4mu   = new TH1D("h1D_m4l_4P_4mu","h1D_m4l_4P_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_4P_4mu->Sumw2();
    
    
    for (int s = mass4lREFIT; s<=Djet; s++)
    {
        CR_shapes[s] = new TH1D("CR_"+Shape_names[s],"CR_"+Shape_names[s],CR_nBins[s], CR_min[s], CR_max[s]); CR_shapes[s]->Sumw2();
        CR_shapes_4e[s] = new TH1D("CR_4e_"+Shape_names[s],"CR_4e"+Shape_names[s],CR_nBins[s], CR_min[s], CR_max[s]); CR_shapes_4e[s]->Sumw2();
        CR_shapes_4mu[s] = new TH1D("CR_4mu"+Shape_names[s],"CR_4mu"+Shape_names[s],CR_nBins[s], CR_min[s], CR_max[s]); CR_shapes_4mu[s]->Sumw2();
        CR_shapes_2e2mu[s] = new TH1D("CR_2e2mu"+Shape_names[s],"CR_2e2mu"+Shape_names[s],CR_nBins[s], CR_min[s], CR_max[s]); CR_shapes_2e2mu[s]->Sumw2();
    }
    
    
    
    
    
    RM2 = new TH2D("a","a",40,0,130,20,0,5);
    
    
    // tree for a set of variables, for selected events
    // TTree* slimmedTreeZX = new TTree("selectedEvents","selectedEvents");
    
    getHistFR(zxTree, h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE, h1D_FRel_EB_d, h1D_FRel_EE_d, h1D_FRmu_EB_d, h1D_FRmu_EE_d,ptElCut, ptMuCut, mZ2Cut);
    
    
    
    
    
    TH1D* h1D_FRel_EE_n=(TH1D*)h1D_FRel_EE->Clone();
    TH1D* h1D_FRmu_EE_n=(TH1D*)h1D_FRmu_EE->Clone();
    TH1D* h1D_FRel_EB_n=(TH1D*)h1D_FRel_EB->Clone();
    TH1D* h1D_FRmu_EB_n=(TH1D*)h1D_FRmu_EB->Clone();
    
    
    
    
    TCanvas *c;
    setCavasAndStyles("c",c,"");
    
    
    h1D_FRel_EE->SetFillColor(kBlue);
    h1D_FRel_EE->Draw("e1 goff");
    
    // divide hists to get the fake rates
    h1D_FRel_EB_n->Divide(h1D_FRel_EB_d);
    h1D_FRel_EE_n->Divide(h1D_FRel_EE_d);
    h1D_FRmu_EB_n->Divide(h1D_FRmu_EB_d);
    h1D_FRmu_EE_n->Divide(h1D_FRmu_EE_d);
    
    h1D_FRel_EB_n->GetYaxis()->SetRangeUser(0.01,0.35); 
    h1D_FRel_EE_n->GetYaxis()->SetRangeUser(0.01,0.35);
    h1D_FRmu_EB_n->GetYaxis()->SetRangeUser(0.04,0.35);
    h1D_FRmu_EE_n->GetYaxis()->SetRangeUser(0.04,0.35);
    //Save histograms in .root file
    
    
    h1D_FRel_EB->SetName("Data_FRel_EB_n"); h1D_FRel_EB->Write();
    h1D_FRel_EE->SetName("Data_FRel_EE_n"); h1D_FRel_EE->Write();
    h1D_FRmu_EB->SetName("Data_FRmu_EB_n"); h1D_FRmu_EB->Write();
    h1D_FRmu_EE->SetName("Data_FRmu_EE_n"); h1D_FRmu_EE->Write();
    
    h1D_FRel_EB_d->SetName("Data_FRel_EB_d"); h1D_FRel_EB_d->Write();
    h1D_FRel_EE_d->SetName("Data_FRel_EE_d"); h1D_FRel_EE_d->Write();
    h1D_FRmu_EB_d->SetName("Data_FRmu_EB_d"); h1D_FRmu_EB_d->Write();
    h1D_FRmu_EE_d->SetName("Data_FRmu_EE_d"); h1D_FRmu_EE_d->Write();
    
    
    h1D_FRel_EB_n->SetName("Data_FRel_EB"); h1D_FRel_EB_n->Write();
    h1D_FRel_EE_n->SetName("Data_FRel_EE"); h1D_FRel_EE_n->Write();
    h1D_FRmu_EB_n->SetName("Data_FRmu_EB"); h1D_FRmu_EB_n->Write();
    h1D_FRmu_EE_n->SetName("Data_FRmu_EE"); h1D_FRmu_EE_n->Write();
    
    // store slimmed tree and histograms in .root file
    
    h1D_m4l_2P2F->SetName("h1D_m4l_2P2F"); h1D_m4l_2P2F->Write();
    h1D_m4l_2P2F_2e2mu->SetName("h1D_m4l_2P2F_2e2mu"); h1D_m4l_2P2F_2e2mu->Write();
    h1D_m4l_2P2F_2mu2e->SetName("h1D_m4l_2P2F_2mu2e"); h1D_m4l_2P2F_2mu2e->Write();
    h1D_m4l_2P2F_4e->SetName("h1D_m4l_2P2F_4e"); h1D_m4l_2P2F_4e->Write();
    h1D_m4l_2P2F_4mu->SetName("h1D_m4l_2P2F_4mu"); h1D_m4l_2P2F_4mu->Write();
    
    h1D_m4l_3P1F->SetName("h1D_m4l_3P1F"); h1D_m4l_3P1F->Write();
    h1D_m4l_3P1F_2e2mu->SetName("h1D_m4l_3P1F_2e2mu"); h1D_m4l_3P1F_2e2mu->Write();
    h1D_m4l_3P1F_2mu2e->SetName("h1D_m4l_3P1F_2mu2e"); h1D_m4l_3P1F_2mu2e->Write();
    h1D_m4l_3P1F_4e->SetName("h1D_m4l_3P1F_4e"); h1D_m4l_3P1F_4e->Write();
    h1D_m4l_3P1F_4mu->SetName("h1D_m4l_3P1F_4mu"); h1D_m4l_3P1F_4mu->Write();
    
    
    
    // plot the histograms
    // setup environment & canvas
    TCanvas *c1;
    setCavasAndStyles("c1",c1,"");
    TLegend* leg = 0;
    
    RM2->GetXaxis()->SetTitle("m_{#font[12]{Z_2}} (GeV)");
    RM2->GetXaxis()->SetRange(3,40);
    RM2->GetXaxis()->SetLabelFont(42);
    RM2->GetXaxis()->SetLabelSize(0.05);
    RM2->GetXaxis()->SetTitleSize(0.06);
    RM2->GetXaxis()->SetTitleOffset(0.9);
    RM2->GetXaxis()->SetTitleFont(42);
    RM2->GetYaxis()->SetTitle("#Delta R");
    RM2->GetYaxis()->SetLabelFont(42);
    RM2->GetYaxis()->SetLabelSize(0.05);
    RM2->GetYaxis()->SetTitleSize(0.06);
    RM2->GetYaxis()->SetTitleOffset(1.25);
    RM2->GetYaxis()->SetTitleFont(42);
    RM2->GetZaxis()->SetLabelFont(42);
    RM2->GetZaxis()->SetLabelSize(0.035);
    RM2->GetZaxis()->SetTitleSize(0.035);
    RM2->GetZaxis()->SetTitleFont(42);
    
    
    RM2->Draw("COL goff");
    TLatex *   tex = new TLatex(0.95,0.96,"2.8 fb^{-1} (13 TeV)");
    tex->SetNDC();
    tex->SetTextAlign(31);
    tex->SetTextFont(42);
    tex->SetTextSize(0.03);
    tex->SetLineWidth(2);
    tex->Draw("goff");
    tex = new TLatex(0.15,0.96,"CMS");
    tex->SetNDC();
    tex->SetTextFont(61);
    tex->SetTextSize(0.0375);
    tex->SetLineWidth(2);
    tex->Draw("goff");
    tex = new TLatex(0.23,0.96,"Preliminary");
    tex->SetNDC();
    tex->SetTextFont(52);
    tex->SetTextSize(0.0285);
    tex->SetLineWidth(2);
    tex->Draw("goff");
    c1->SaveAs(sPlotsStore+"/2DRM2.pdf");
    
    
    setHistProperties(h1D_dummy,1,1,kBlack,0,0,varAxLabel,"Fake Rate");
    
    //For data only
    // Fake rates
    // electrons FR
    
    /* if (varName == "ptl3"){
     
     TFile *f = new TFile("Hist_Data_ptl3_WZremoved.root");
     
     TH1D *EB = (TH1D*) f->Get("Data_FRel_EB"); EB->SetLineColor(kRed+2);
     TH1D *EE = (TH1D*) f->Get("Data_FRel_EE"); EE->SetLineColor(kBlue+2);
     TH1D *MB = (TH1D*) f->Get("Data_FRmu_EB"); MB->SetLineColor(kRed+2);
     TH1D *ME = (TH1D*) f->Get("Data_FRmu_EE"); ME->SetLineColor(kBlue+2);
     
     //h1D_dummy->SetMaximum(4.0*h1D_FRel_EE_n->GetMaximum());
     h1D_dummy->SetMaximum(0.3);
     h1D_dummy->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
     setHistProperties(h1D_FRel_EB_n,lineWidth,1,kRed-7); h1D_FRel_EB_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_EB_n, "electrons, EB","L");
     setHistProperties(h1D_FRel_EE_n,lineWidth,1,kBlue-7); h1D_FRel_EE_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_EE_n, "electrons, EE","L");
     
     setHistProperties(EB,lineWidth,1,kRed+2); EB->Draw("same E1 goff"); leg->AddEntry(EB, "electrons, EB WZ removed","L");
     setHistProperties(EE,lineWidth,1,kBlue+2); EE->Draw("same E1 goff"); leg->AddEntry(EE, "electrons, EE WZ removed","L");
     
     
     
     leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_electrons_"+varName+"_"+fileName+"all.pdf");
     
     
     // muons FR
     //h1D_dummy->SetMaximum(4.0*h1D_FRmu_EE_n->GetMaximum());
     h1D_dummy->SetMaximum(0.3);
     h1D_dummy->Draw(); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
     setHistProperties(h1D_FRmu_EB_n,lineWidth,1,kRed-7); h1D_FRmu_EB_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_EB_n, "muons, MB","L");
     setHistProperties(h1D_FRmu_EE_n,lineWidth,1,kBlue-7); h1D_FRmu_EE_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_EE_n, "muons, ME","L");
     
     setHistProperties(MB,lineWidth,1,kRed+2); MB->Draw("same E1 goff"); leg->AddEntry(MB, "muons, MB WZ removed","L");
     setHistProperties(ME,lineWidth,1,kBlue+2); ME->Draw("same E1 goff"); leg->AddEntry(ME, "muons, ME WZ removed","L");
     
     leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_muons_"+varName+"_"+fileName+"all.pdf");
     
     
     
     
     
     }*/
    
    
    h1D_dummy->SetMaximum(0.35);
    h1D_dummy->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRel_EB_n,lineWidth,1,kRed-7); h1D_FRel_EB_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_EB_n, "electrons, EB","L");
    setHistProperties(h1D_FRel_EE_n,lineWidth,1,kBlue-7); h1D_FRel_EE_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRel_EE_n, "electrons, EE","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_electrons_"+varName+"_"+fileName+".pdf");
    // muons FR
    h1D_dummy->SetMaximum(0.35);
    h1D_dummy->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRmu_EB_n,lineWidth,1,kRed-7); h1D_FRmu_EB_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_EB_n, "muons, MB","L");
    setHistProperties(h1D_FRmu_EE_n,lineWidth,1,kBlue-7); h1D_FRmu_EE_n->Draw("same E1 goff"); leg->AddEntry(h1D_FRmu_EE_n, "muons, ME","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/FR_muons_"+varName+"_"+fileName+".pdf");
    
    
    
    setHistProperties(h1D_dummy,1,1,kBlack,0,0,varAxLabel,"Events/Bin Width");
    
    // Nominator and denominator
    // electrons FR
    //if (isData)
    //{
    if (h1D_FRel_EB->GetMaximum()>h1D_FRel_EE->GetMaximum()) HistSizeel=2*h1D_FRel_EB->GetMaximum();
    else HistSizeel=2*h1D_FRel_EE->GetMaximum();
    if (h1D_FRmu_EB->GetMaximum()>h1D_FRmu_EE->GetMaximum()) HistSizemu=2*h1D_FRmu_EB->GetMaximum();
    else HistSizemu=2*h1D_FRmu_EE->GetMaximum();
    
    if (h1D_FRel_EB_d->GetMaximum()>h1D_FRel_EE_d->GetMaximum()) HistSizeel_d=2*h1D_FRel_EB_d->GetMaximum();
    else HistSizeel=2*h1D_FRel_EE->GetMaximum();
    if (h1D_FRmu_EB_d->GetMaximum()>h1D_FRmu_EE_d->GetMaximum()) HistSizemu_d=2*h1D_FRmu_EB_d->GetMaximum();
    else HistSizemu_d=2*h1D_FRmu_EE_d->GetMaximum();
    // }
    
    
    
    h1D_dummy->SetMaximum(HistSizeel);
    h1D_dummy->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRel_EB,lineWidth,1,kRed-7);h1D_FRel_EB->SetLineColor(kRed-7);h1D_FRel_EB->SetFillColor(0); h1D_FRel_EB->Draw("hist same"); leg->AddEntry(h1D_FRel_EB_n, "electrons, Barrel","L");
    setHistProperties(h1D_FRel_EE,lineWidth,1,kBlue-7);h1D_FRel_EE->SetLineColor(kBlue-7);h1D_FRel_EE->SetFillColor(0); h1D_FRel_EE->Draw("hist same"); leg->AddEntry(h1D_FRel_EE_n, "electrons, Endcap","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/N_electrons_"+varName+"_"+fileName+".pdf");
    // muons FR
    h1D_dummy->SetMaximum(HistSizemu);
    h1D_dummy->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRmu_EB,lineWidth,1,kRed-7); h1D_FRmu_EB->SetLineColor(kRed-7);h1D_FRmu_EB->SetFillColor(0); h1D_FRmu_EB->Draw("hist same"); leg->AddEntry(h1D_FRmu_EB_n, "muons, Barrel","L");
    setHistProperties(h1D_FRmu_EE,lineWidth,1,kBlue-7); h1D_FRmu_EE->SetLineColor(kBlue-7);h1D_FRmu_EE->SetFillColor(0); h1D_FRmu_EE->Draw("hist same"); leg->AddEntry(h1D_FRmu_EE_n, "muons, Endcap","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/N_muons_"+varName+"_"+fileName+".pdf");
    
    h1D_dummy->SetMaximum(HistSizeel_d);
    h1D_dummy->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRel_EB_d,lineWidth,1,kRed-7);h1D_FRel_EB_d->SetLineColor(kRed-7);h1D_FRel_EB_d->SetFillColor(0); h1D_FRel_EB_d->Draw("hist same"); leg->AddEntry(h1D_FRel_EB_d, "electrons, Barrel","L");
    setHistProperties(h1D_FRel_EE_d,lineWidth,1,kBlue-7);h1D_FRel_EE_d->SetLineColor(kBlue-7);h1D_FRel_EE_d->SetFillColor(0); h1D_FRel_EE_d->Draw("hist same"); leg->AddEntry(h1D_FRel_EE_n, "electrons, Endcap","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/D_electrons_"+varName+"_"+fileName+".pdf");
    // muons FR
    h1D_dummy->SetMaximum(HistSizemu_d);
    h1D_dummy->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"pp #rightarrow Z + l_{loose}");
    setHistProperties(h1D_FRmu_EB_d,lineWidth,1,kRed-7); h1D_FRmu_EB_d->SetLineColor(kRed-7);h1D_FRmu_EB_d->SetFillColor(0); h1D_FRmu_EB_d->Draw("hist same"); leg->AddEntry(h1D_FRmu_EB_n, "muons, Barrel","L");
    setHistProperties(h1D_FRmu_EE_d,lineWidth,1,kBlue-7); h1D_FRmu_EE_d->SetLineColor(kBlue-7);h1D_FRmu_EE_d->SetFillColor(0); h1D_FRmu_EE_d->Draw("hist same"); leg->AddEntry(h1D_FRmu_EE_n, "muons, Endcap","L");
    leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/D_muons_"+varName+"_"+fileName+".pdf");
    
    
    
    
    
    
    
    
    
    TEfficiency *teff = new TEfficiency();
    double upperLimit_257969 = teff->ClopperPearson(d_257969, n_257969, 0.683, true);
    double lowerLimit_257969 = teff->ClopperPearson(d_257969, n_257969, 0.683, false);
    
    double upperLimit_258440 = teff->ClopperPearson(d_258440, n_258440, 0.683, true);
    double lowerLimit_258440 = teff->ClopperPearson(d_258440, n_258440, 0.683, false);
    
    double upperLimit_258749 = teff->ClopperPearson(d_258749, n_258749, 0.683, true);
    double lowerLimit_258749 = teff->ClopperPearson(d_258749, n_258749, 0.683, false);
    
    double upperLimit_259891 = teff->ClopperPearson(d_259891, n_259891, 0.683, true);
    double lowerLimit_259891 = teff->ClopperPearson(d_259891, n_259891, 0.683, false);
    
    
    if (isData)
    {
        
        TCanvas *c2 = new TCanvas("c1","FR split into groups with ~400pb^{-1}:"+fileName+"_"+varName,200,10,700,500);
        c2->SetFillColor(10);
        c2->SetGrid();
        c2->GetFrame()->SetFillColor(21);
        c2->GetFrame()->SetBorderSize(12);
        const Int_t n = 4;
        Double_t x[n]  = {l_257969/2.0,l_257969+l_258440/2.0,l_257969+l_258440+l_258749/2.0,l_257969+l_258440+l_258749+l_259891/2.0};
        Double_t y[n]  ={1.0*n_257969/d_257969,1.0*n_258440/d_258440,1.0*n_258749/d_258749,1.0*n_259891/d_259891};
        
        Double_t ex_h[n] = {l_257969/2.0,l_258440/2.0,l_258749/2.0,l_259891/2.0};
        Double_t ex_l[n] = {l_257969/2.0,l_258440/2.0,l_258749/2.0,l_259891/2.0};
        
        
        Double_t ey_h[n] = {upperLimit_257969-1.0*n_257969/d_257969,upperLimit_258440-1.0*n_258440/d_258440,upperLimit_258749-1.0*n_258749/d_258749,upperLimit_259891-1.0*n_259891/d_259891};
        Double_t ey_l[n] = {1.0*n_257969/d_257969-lowerLimit_257969,1.0*n_258440/d_258440-lowerLimit_258440,1.0*n_258749/d_258749-lowerLimit_258749,1.0*n_259891/d_259891-lowerLimit_259891};
        TGraphAsymmErrors *gr = new TGraphAsymmErrors(n,x,y,ex_l,ex_h,ey_l,ey_h);
        
        gr->SetTitle("TGraphErrors Example");
        gr->SetMarkerColor(10);
        gr->SetMarkerStyle(10);
        gr->SetMaximum(0.05);
        gr->SetMinimum(0);
        gr->Draw("AP goff");
        cout<<"1 "<<1.0*n_257969/d_257969<<" 2 "<<1.0*n_258440/d_258440<<" 3 "<<1.0*n_258749/d_258749<<" 4 "<<1.0*n_259891/d_259891<<endl;
        c2->SaveAs(sPlotsStore+"/FR_RunNum_"+varName+"_"+fileName+".pdf");
    }
    
    // define dummy histogram for CRs
    var_plotHigh = 2000.0;  var_plotLow = 0.0;  var_nBins = 100;
    TH1D* h1D_dummyCR = new TH1D("dummyCR", "dummyCR", var_nBins, var_plotLow, var_plotHigh);
    varAxLabel = "m_{4l} [GeV]";
    binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummyCR,1,1,kBlack,0,0,varAxLabel,"Events/"+sBinWidth);
    
    
    // Control Regions
    // CR 2P2F
    TCanvas *c2;
    setCavasAndStyles("c2",c2,"");
    TLegend* leg2 = 0;
    
    if (isData and varName == "ptl3"){
        
        // CR 3P1F
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_3P1F->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 3P+1F");
        setHistProperties(h1D_m4l_3P1F,lineWidth,1,kRed-7); h1D_m4l_3P1F->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_3P1F, "4l","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_3P1F.pdf");
        // CR 4P
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_4P->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Signal Region");
        setHistProperties(h1D_m4l_4P,lineWidth,1,kRed-7); h1D_m4l_4P->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_4P, "4l","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_4P.pdf");
        
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_2P2F->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 2P+2F");
        setHistProperties(h1D_m4l_2P2F,lineWidth,1,kRed-7); h1D_m4l_2P2F->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_2P2F, "4l","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_2P2F.pdf");
        
        
        
        //4e
        
        // CR 3P1F_4e
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_3P1F_4e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 3P+1F");
        setHistProperties(h1D_m4l_3P1F_4e,lineWidth,1,kRed-7); h1D_m4l_3P1F_4e->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_3P1F_4e, "4e","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_3P1F_4e.pdf");
        // CR 4P
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_4P_4e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Signal Region");
        setHistProperties(h1D_m4l_4P_4e,lineWidth,1,kRed-7); h1D_m4l_4P_4e->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_4P_4e, "4e","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_4P_4e.pdf");
        
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_2P2F_4e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 2P+2F");
        setHistProperties(h1D_m4l_2P2F_4e,lineWidth,1,kRed-7); h1D_m4l_2P2F_4e->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_2P2F_4e, "4e","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_2P2F_4e.pdf");
        
        
        
        //4mu
        // CR 3P1F_4mu
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_3P1F_4mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 3P+1F");
        setHistProperties(h1D_m4l_3P1F_4mu,lineWidth,1,kRed-7); h1D_m4l_3P1F_4mu->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_3P1F_4mu, "4mu","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_3P1F_4mu.pdf");
        // CR 4P
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_4P_4mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Signal Region");
        setHistProperties(h1D_m4l_4P_4mu,lineWidth,1,kRed-7); h1D_m4l_4P_4mu->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_4P_4mu, "4mu","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_4P_4mu.pdf");
        
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_2P2F_4mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 2P+2F");
        setHistProperties(h1D_m4l_2P2F_4mu,lineWidth,1,kRed-7); h1D_m4l_2P2F_4mu->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_2P2F_4mu, "4mu","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_2P2F_4mu.pdf");
        
        
        //2e2mu
        // CR 2P2F_2e2mu_2e2mu
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_3P1F_2e2mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 3P+1F");
        setHistProperties(h1D_m4l_3P1F_2e2mu,lineWidth,1,kRed-7); h1D_m4l_3P1F_2e2mu->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_3P1F_2e2mu, "2e2mu","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_3P1F_2e2mu.pdf");
        // CR 4P
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_4P_2e2mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Signal Region");
        setHistProperties(h1D_m4l_4P_2e2mu,lineWidth,1,kRed-7); h1D_m4l_4P_2e2mu->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_4P_2e2mu, "2e2mu","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_4P_2e2mu.pdf");
        
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_2P2F_2e2mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c2); leg2 = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg2,"Region 2P+2F");
        setHistProperties(h1D_m4l_2P2F_2e2mu,lineWidth,1,kRed-7); h1D_m4l_2P2F_2e2mu->Draw("same E1 goff"); leg2->AddEntry(h1D_m4l_2P2F_2e2mu, "2e2mu","L");
        leg2->Draw("goff"); c2->SaveAs(sPlotsStore+"/CR_2P2F_2e2mu.pdf");
        
        
        
        
        
        
        
        
        
    }
    
    
    
    
}



void IDfinder(vector<int> *lep_matchedR03_PdgId,vector<int> *lep_matchedR03_MomId,vector<int> *lep_matchedR03_MomMomId )
{
    
    
    
    
}
int getHistFR(TChain* tree,
              TH1D* &h1D_FRel_EB,   TH1D* &h1D_FRel_EE,   TH1D* &h1D_FRmu_EB,   TH1D* &h1D_FRmu_EE,
              TH1D* &h1D_FRel_EB_d, TH1D* &h1D_FRel_EE_d, TH1D* &h1D_FRmu_EB_d, TH1D* &h1D_FRmu_EE_d,
              double ptElCut, double ptMuCut, double mZ2Cut){
    
    
    // define vars and branches
    float mass3l, mass4l, massZ1, massZ2, met, met_phi,mt;
    float eventWeight, dataMCWeight, crossSection, genWeight;
    long int Run, LumiSect, Event;
    int finalState, nVtx;
    bool passedFullSelection, passedZ1LSelection, passedZ4lZXCRSelection, passedZ4lSelection,passedZXCRSelection, passedZ4lZ1LSelection;
    //    TClonesArray *lep_p4 = new TClonesArray("TLorentzVector", 10);
    //    TClonesArray *met_p4 = new TClonesArray("TLorentzVector", 1);
    vector<float> *lep_pt = 0; TBranch *b_lep_pt = 0;
    vector<float> *lep_eta = 0; TBranch *b_lep_eta = 0;
    vector<float> *lep_phi = 0; TBranch *b_lep_phi = 0;
    vector<float> *lep_mass = 0; TBranch *b_lep_mass = 0;
    int lep_Hindex[4];
    vector<int> *lep_id = 0; TBranch *b_lep_id = 0;
    vector<int> *lep_tightId = 0; TBranch *b_lep_tightId = 0;
    vector<float> *lep_RelIso = 0; TBranch *b_lep_RelIso = 0;
    vector<int> *lep_matchedR03_PdgId = 0; TBranch *b_lep_matchedR03_PdgId = 0;
    vector<int> *lep_matchedR03_MomId = 0; TBranch *b_lep_matchedR03_MomId = 0;
    vector<int> *lep_matchedR03_MomMomId = 0; TBranch *b_lep_matchedR03_MomMomId = 0;
    
    
    
    
    // counters
    int nFailedIDOnlyLeptonsZ2 = 0;
    int nEvt2P2FLeptonsIDonly_lt06 = 0;
    int nEvt2P2FLeptonsIsoID_lt06 = 0;
    int nEvtPassedZ1LSelection = 0;
    int nEvtPassedZXCRSelection = 0;
    int nEvtZ1WithFailedLeptons = 0;
    int nEvt1P3FLeptons = 0;
    int nEvt2P2FLeptons = 0;
    int nEvt3P1FLeptons = 0;
    int nEvt4PLeptons = 0;
    
    int nEvt1P3FLeptons4e = 0;
    int nEvt2P2FLeptons4e = 0;
    int nEvt3P1FLeptons4e = 0;
    int nEvt4PLeptons4e = 0;
    
    int nEvt1P3FLeptons4mu = 0;
    int nEvt2P2FLeptons4mu = 0;
    int nEvt3P1FLeptons4mu = 0;
    int nEvt4PLeptons4mu = 0;
    
    int nEvt1P3FLeptons2e2mu = 0;
    int nEvt2P2FLeptons2e2mu = 0;
    int nEvt3P1FLeptons2e2mu = 0;
    int nEvt4PLeptons2e2mu = 0;
    
    
    
    // get branches
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);
    tree->SetBranchAddress("crossSection",&crossSection);
    tree->SetBranchAddress("eventWeight",&eventWeight);
    tree->SetBranchAddress("genWeight",&genWeight);
    tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
    tree->SetBranchAddress("passedZ4lZ1LSelection",&passedZ4lZ1LSelection);
    tree->SetBranchAddress("passedZ4lZXCRSelection",&passedZ4lZXCRSelection);
    tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
    tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
    tree->SetBranchAddress("nVtx",&nVtx);
    tree->SetBranchAddress("finalState",&finalState);
    tree->SetBranchAddress("mass4l",&mass4l);
    //tree->SetBranchAddress("massZ1",&massZ1);
    tree->SetBranchAddress("massZ2",&massZ2);
    
   // tree->SetBranchAddress("lNEvents",&lNEvents);
    
    //    tree->SetBranchAddress("mass3l",&mass3l);
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    //    tree->GetBranch("lep_p4")->SetAutoDelete(kFALSE);
    //    tree->SetBranchAddress("lep_p4",&lep_p4);
    tree->SetBranchAddress("lep_pt",&lep_pt,&b_lep_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta,&b_lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi,&b_lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass,&b_lep_mass);
    //    tree->SetBranchAddress("met_p4",&met_p4);
    tree->SetBranchAddress("met",&met);
    tree->SetBranchAddress("met_phi",&met_phi);
    tree->SetBranchAddress("lep_id",&lep_id,&b_lep_id);
    tree->SetBranchAddress("lep_tightId",&lep_tightId,&b_lep_tightId);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIso,&b_lep_RelIso);
    tree->SetBranchAddress("lep_matchedR03_PdgId",&lep_matchedR03_PdgId,&b_lep_matchedR03_PdgId);
    tree->SetBranchAddress("lep_matchedR03_MomId",&lep_matchedR03_MomId,&b_lep_matchedR03_MomId);
    tree->SetBranchAddress("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId,&b_lep_matchedR03_MomMomId);
    
    
    
    // fill histograms
    Long64_t nentries = tree->GetEntries();
    //    nentries = 1000;
    cout << "nentries: " << nentries << endl;
    
    
    
    
    if(!SILENT) {
        cout << setw(printOutWidth) << "nEntries: " << nentries << setw(printOutWidth) << endl;
        cout << "Event" << ":" << "Run" << ":" << "LumiSect" << ":" << "mass3l" << ":" << "massZ1" << ":" << "weight" << endl;
    }
    // sort tree in "Event" number, ascending
    int *index = new int[nentries];
    if (SORT_EVENTS) {
        tree->Draw("Event","","goff");
        int nentries_int = static_cast<int>(nentries);
        TMath::Sort(nentries_int, tree->GetV1(), index, false);
    }
    
    double RunFR_n=0,RunFR_d=0;
    int SkippedEvents=0;
    long int prevEvent = -1;
    long int prevRun = -1;
   for(int iEvt=0; iEvt < nentries; iEvt++){
    //for(int iEvt=0; iEvt < 100000; iEvt++){
        if (SORT_EVENTS) {
            tree->GetEntry(index[iEvt]);	//index[iEvt]);}// take sorted entries
        } else {
            tree->GetEntry(iEvt);
        }// take unsroted entries
        /*
         cout<<"Size pdg id =  "<<lep_matchedR03_PdgId->size()<<endl;
         cout<<"Size Mom id = "<<lep_matchedR03_MomId->size()<<endl;
         cout<<"Size MomMom id = "<<lep_matchedR03_MomMomId->size()<<endl;
         */
        //        cout<<"RunNum= "<< Run<<"; Event= "<<Event<<endl;
        /*  if ((prevEvent==Event)&&(prevRun==Run)) {
         SkippedEvents++;
         //     cout<<"Skipping event " << Event <<"; passedFullSelection "<<passedFullSelection<<"; passedZ1LSelection "<<passedZ1LSelection<<"; passedZXCRSelection "<<passedZXCRSelection<<"; passedZ4lSelection "<<passedZ4lSelection<<"; finalState "<<finalState<<"; mass4l "<<mass4l<<endl;
         cout<<"Skipping event " << Event<<"; Run "<<Run <<endl;
         continue;
         } else {
         prevEvent = Event;
         prevRun = Run;
         }*/
        
        // print out
        if(iEvt%50000==0) cout << "   event: " << iEvt << "/" << nentries << endl;
        
        // weight
        float weight = eventWeight*dataMCWeight;
        
        if (isData) weight=1;
        
        //  cout<<lNEvents<<endl;
        
        //double Lumi=1547.18;
        if (!isData){
            if (fileName=="DY10")
                // weight
                weight = genWeight*18610.0*LUMI_INT/lNEvents;
            
            //weight = eventWeight*dataMCWeight*crossSection*LUMI_INT/lNEvents;
            
            if (fileName=="DY50")
                // weight
                weight = genWeight*6104.0*LUMI_INT/lNEvents;
            if (fileName=="TT")
                // weight
                weight = genWeight*815.96*LUMI_INT/lNEvents;
            if (fileName=="WZ")
                // weight
                weight = genWeight*4.430*LUMI_INT/lNEvents;
            if (fileName=="ZZ")
                // weight
                weight = genWeight*1.256*LUMI_INT/lNEvents;
        }
        bool passedCRsel = false;
        bool passedZ1Lsel = false;
        
        if (isData) passedCRsel =passedZXCRSelection;
       // if (!isData) passedCRsel = passedZ4lZXCRSelection;
        if (!isData) passedCRsel =passedZXCRSelection;
        if (isData) passedZ1Lsel =passedZ1LSelection;
     //   if (!isData) passedZ1Lsel = passedZ4lZ1LSelection;
        if (!isData) passedZ1Lsel = passedZ1LSelection; 
        //     cout<<weight<<endl;
        
        if (passedZ1Lsel) {
            
            // if(!isData)   weight*=(1.37/1.23);
            
            nEvtPassedZ1LSelection++;
            //        cout << "   nEvtPassedZ1LSelection: " << nEvtPassedZ1LSelection << "/" << iEvt << endl;
            //        cout << Event << ":" << Run << ":" << LumiSect << ":" << massZ1 << ":" << passedZ1LSelection << ":" << lep_Hindex << ":" << b_lep_id << ":" << lep_p4 << endl;
            //        cout << "iEvt: " << iEvt << ", lep_Hindex[2]: " << lep_Hindex[2] << ", lep_id->size():" << lep_id->size() << ", lep_tightId->size():" << lep_tightId->size() << endl;
            
            // get properties of the probe lepton (3rd)
            int lep_tight = lep_tightId->at(lep_Hindex[2]);
            double lep_iso = lep_RelIso->at(lep_Hindex[2]);
            int idL3 = lep_id->at(lep_Hindex[2]);
            TLorentzVector lep;
            lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[2]),lep_eta->at(lep_Hindex[2]),lep_phi->at(lep_Hindex[2]),lep_mass->at(lep_Hindex[2]));
            float pTL3  = lep.Pt();
            float etaL3 = lep.Eta();
            float phiL3 = lep.Phi();
            mt = sqrt((pTL3 + met )*(pTL3 + met )- (pTL3*sin(phiL3)  + met*sin( met_phi))*(pTL3*sin(phiL3)  + met*sin( met_phi))- (pTL3*cos(phiL3) + met*cos(met_phi))*(pTL3*cos(phiL3) + met*cos(met_phi)));
            
            
            TLorentzVector lep_1, lep_2;
            lep_1.SetPtEtaPhiM(lep_pt->at(lep_Hindex[0]),lep_eta->at(lep_Hindex[0]),lep_phi->at(lep_Hindex[0]),lep_mass->at(lep_Hindex[0]));
            lep_2.SetPtEtaPhiM(lep_pt->at(lep_Hindex[1]),lep_eta->at(lep_Hindex[1]),lep_phi->at(lep_Hindex[1]),lep_mass->at(lep_Hindex[1]));
            
            massZ1 = (lep_1+lep_2).M();
            
            //l3summ=(*lep1+*lep2+*lep3);
            //LMet=(*met_p4);
            // double m3l=l3summ.M();
            double m3l = 0;
            // TLorentzVector metl3=(*mEt+*lep3);
            //double mt=sqrt(metl3.E()*metl3.E()-metl3.X()*metl3.X()-metl3.Y()*metl3.Y());
            double mt = 0;
            
            bool TestVar=false;
            double FillVar=0;
            
            if (varName=="mZ1"){
                TestVar= 1;
                FillVar= massZ1;
            }
            if (varName=="mEt"){
                TestVar= 1;
                FillVar= met;
            }
            if (varName=="nVtx"){
                TestVar= (met<25)&&(fabs(massZ1-91.188)<7);
                FillVar= nVtx;
            }
            if (varName=="ptl3"){
                TestVar= (fabs(massZ1-91.188)<7)&&(met<25) ;
                FillVar= pTL3;
            }
            if (varName=="etal3") {
                TestVar= (met<25)&&fabs(massZ1-91.188)<7;
                FillVar= etaL3;
            }
            
            if (varName=="m3l"){
                TestVar= (met<25)&&fabs(massZ1-91.188)<7;
                FillVar= m3l;
            }
            
            if (varName=="mt"){
                TestVar= 1;
                FillVar= mt;
            }
            
            
            // std::cout<<"Inside Z1L"<<met<<std::endl;
            // fill 1D hists
            if ((abs(idL3) == 11) && (fabs(etaL3) < 1.497)&&TestVar)
                
            {
                //RunNumCalc(false,Run);
                h1D_FRel_EB_d->Fill(FillVar, weight);
                // if (TestVar)
                  PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[EB_d], Hist_fakes[EB_d], Hist_BDfakes[EB_d], Hist_conv[EB_d],false);
                //if (lep_tight && (lep_iso<0.35)&&TestVar)
                if (lep_tight && TestVar)
                    
                {
                    
                         PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[EB_n], Hist_fakes[EB_n],Hist_BDfakes[EB_n],Hist_conv[EB_n],false);
                    RunNumCalc(true,Run);
                    h1D_FRel_EB->Fill(FillVar, weight);
                }
                
            }
            if ((abs(idL3) == 11) && (fabs(etaL3) > 1.497)&&TestVar) {
                RunNumCalc(false,Run);
                h1D_FRel_EE_d->Fill(FillVar, weight);
                
                // if (TestVar)
                  PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[EE_d], Hist_fakes[EE_d],Hist_BDfakes[EE_d], Hist_conv[EE_d],false);
                
             //   if (lep_tight && (lep_iso<0.35)&&TestVar)
               if (lep_tight && TestVar)
		 {
                    
                    
                    PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[EE_n], Hist_fakes[EE_n], Hist_BDfakes[EE_n], Hist_conv[EE_n],false);
                    //RunNumCalc(true,Run);
                    h1D_FRel_EE->Fill(FillVar, weight);
                }
            }
            if ((abs(idL3) == 13) && (fabs(etaL3) < 1.2)&&TestVar) {
                //RunNumCalc(false,Run);
                h1D_FRmu_EB_d->Fill(FillVar, weight);
                // if (TestVar)
                  PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[MB_d], Hist_fakes[MB_d], Hist_BDfakes[MB_d], Hist_conv[MB_d],false);
                if (lep_tight && (lep_iso<0.35)&&TestVar){
                    // RunNumCalc(true,Run);
                     PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[MB_n], Hist_fakes[MB_n],Hist_BDfakes[MB_n], Hist_conv[MB_n],false);
                    h1D_FRmu_EB->Fill(FillVar, weight);
                }
            }
            if ((abs(idL3) == 13) && (fabs(etaL3) > 1.2)&&TestVar) {
                //RunNumCalc(false,Run);
                h1D_FRmu_EE_d->Fill(FillVar, weight);
                
                  PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[ME_d], Hist_fakes[ME_d],Hist_BDfakes[ME_d], Hist_conv[ME_d],false);
                
                if (lep_tight && (lep_iso<0.35)&&TestVar){
                    // RunNumCalc(true,Run);
                     PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, FillVar,weight,Hist_prompt[ME_n], Hist_fakes[ME_n],Hist_BDfakes[ME_n], Hist_conv[ME_n],false);
                    h1D_FRmu_EE->Fill(FillVar, weight);
                }
            }
        }
        
        
        
        else  if (passedCRsel && varName == "ptl3") {
            nEvtPassedZXCRSelection++;
            
            // if (!isData)   weight = genWeight*crossSection*LUMI_INT/lNEvents;
            //        cout << "   nEvtPassedZ1LSelection: " << nEvtPassedZ1LSelection << "/" << iEvt << endl;
            //        cout << Event << ":" << Run << ":" << LumiSect << ":" << massZ1 << ":" << passedZ1LSelection << ":" << lep_Hindex << ":" << b_lep_id << ":" << lep_p4 << endl;
            //        cout << "iEvt: " << iEvt << ", lep_Hindex[2]: " << lep_Hindex[2] << ", lep_id->size():" << lep_id->size() << ", lep_tightId->size():" << lep_tightId->size() << endl;
            
            // get properties of the non-Z1 leptons (3rd, 4th)
            int lep_tight[4];
            double lep_iso[4];
            int idL[4];
            for(unsigned int k = 0; k <= 3; k++) {
                lep_tight[k] = lep_tightId->at(lep_Hindex[k]);
                lep_iso[k]= lep_RelIso->at(lep_Hindex[k]);
                idL[k] = lep_id->at(lep_Hindex[k]);
            }
            
            
            TLorentzVector lep;
            lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[2]),lep_eta->at(lep_Hindex[2]),lep_phi->at(lep_Hindex[2]),lep_mass->at(lep_Hindex[2]));
            float pTL3  = lep.Pt();
            float etaL3 = lep.Eta();
            float phiL3 = lep.Phi();
            TLorentzVector lep2;
            lep2.SetPtEtaPhiM(lep_pt->at(lep_Hindex[3]),lep_eta->at(lep_Hindex[3]),lep_phi->at(lep_Hindex[3]),lep_mass->at(lep_Hindex[3]));
            float pTL4  = lep2.Pt();
            float etaL4 = lep2.Eta();
            float phiL4 = lep2.Phi();
            float DR = sqrt((phiL3-phiL4)*(phiL3-phiL4) + (etaL3-etaL4)*(etaL3-etaL4));
            
            TLorentzVector lep_1, lep_2;
            lep_1.SetPtEtaPhiM(lep_pt->at(lep_Hindex[0]),lep_eta->at(lep_Hindex[0]),lep_phi->at(lep_Hindex[0]),lep_mass->at(lep_Hindex[0]));
            lep_2.SetPtEtaPhiM(lep_pt->at(lep_Hindex[1]),lep_eta->at(lep_Hindex[1]),lep_phi->at(lep_Hindex[1]),lep_mass->at(lep_Hindex[1]));
            
            massZ1 = (lep_1+lep_2).M();
            
            RM2->Fill(massZ2,DR);
            
            
            
            // count the failed leptons
            //nFailedLeptonsZ1 = !(lep_tight[0] && ((abs(idL[0])==11 && lep_iso[0]<0.35) || (abs(idL[0])==13 && lep_iso[0]<0.35))) +
           // !(lep_tight[1] && ((abs(idL[1])==11 && lep_iso[1]<0.35) || (abs(idL[1])==13 && lep_iso[1]<0.35)));
           // nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11 && lep_iso[2]<0.35) || (abs(idL[2])==13 && lep_iso[2]<0.35))) +
           // !(lep_tight[3] && ((abs(idL[3])==11 && lep_iso[3]<0.35) || (abs(idL[3])==13 && lep_iso[3]<0.35)));

            nFailedLeptonsZ1 = !(lep_tight[0] && ((abs(idL[0])==11) || (abs(idL[0])==13 && lep_iso[0]<0.35))) +
            !(lep_tight[1] && ((abs(idL[1])==11) || (abs(idL[1])==13 && lep_iso[1]<0.35)));
            nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11) || (abs(idL[2])==13 && lep_iso[2]<0.35))) +
            !(lep_tight[3] && ((abs(idL[3])==11) || (abs(idL[3])==13 && lep_iso[3]<0.35)));           


	     nFailedLeptons   = nFailedLeptonsZ1 + nFailedLeptonsZ2;
            
            
            nFailedIDOnlyLeptonsZ2 = !(lep_tight[2]) + !(lep_tight[3]);
            
            
            
            
            
            
            // fill 1D hists
            if (nFailedLeptons == 0) {
                nEvt4PLeptons++;
                h1D_m4l_4P->Fill(mass4l, weight);
                
                
                
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44){
                    
                    h1D_m4l_4P_4e->Fill(mass4l, weight);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                {
                    
                    h1D_m4l_4P_4mu->Fill(mass4l, weight);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                {
                    
                    h1D_m4l_4P_2e2mu->Fill(mass4l, weight);
                    
                }
                
                
                
            }
            if (nFailedLeptons == 1) {
                
                //  if(!isData) weight*=(1.37/1.23);
                
                nEvt3P1FLeptons++;
                h1D_m4l_3P1F->Fill(mass4l, weight);
                 PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_3P1F], CR[Fake][_3P1F],CR[BDfake][_3P1F],CR[Conv][_3P1F],false);
                
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44)
                {
                    CR3P1F_4e<<Run;
                    CR3P1F_4e<<" : ";
                    CR3P1F_4e<<LumiSect;
                    CR3P1F_4e<<" : ";
                    CR3P1F_4e<<Event<<endl;
                    
                       PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_3P1F_4e], CR[Fake][_3P1F_4e],CR[BDfake][_3P1F_4e],CR[Conv][_3P1F_4e],false);
                    
                    h1D_m4l_3P1F_4e->Fill(mass4l, weight);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                {
                    CR3P1F_4mu<<Run;
                    CR3P1F_4mu<<" : ";
                    CR3P1F_4mu<<LumiSect;
                    CR3P1F_4mu<<" : ";
                    CR3P1F_4mu<<Event<<endl;
                    
                     PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_3P1F_4mu], CR[Fake][_3P1F_4mu], CR[BDfake][_3P1F_4mu],CR[Conv][_3P1F_4mu],false);
                    h1D_m4l_3P1F_4mu->Fill(mass4l, weight);
                    
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                {
                    CR3P1F_2e2mu<<Run;
                    CR3P1F_2e2mu<<" : ";
                    CR3P1F_2e2mu<<LumiSect;
                    CR3P1F_2e2mu<<" : ";
                    CR3P1F_2e2mu<<Event<<endl;
                    
                         PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_3P1F_2e2mu], CR[Fake][_3P1F_2e2mu],CR[BDfake][_3P1F_2e2mu],CR[Conv][_3P1F_2e2mu],false);
                    if (abs(idL[2])+abs(idL[3])==22)  h1D_m4l_3P1F_2mu2e->Fill(mass4l, weight);
                    if (abs(idL[2])+abs(idL[3])==26)  h1D_m4l_3P1F_2e2mu->Fill(mass4l, weight);
                    
                }
                
                
            }
            if (nFailedLeptons == 2){
                
                // if(!isData) weight*=(1.53/1.23) ;
                
                nEvt2P2FLeptons++;
                h1D_m4l_2P2F->Fill(mass4l, weight);
                  PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_2P2F], CR[Fake][_2P2F], CR[BDfake][_2P2F],CR[Conv][_2P2F],false);
                
                
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44)
                {//cout<<"OK 4e 2P2F"<<endl;
                    // TString print =Run+" : "+LumiSect+" : "+Event;
                    CR2P2F_4e<<Run;
                    CR2P2F_4e<<" : ";
                    CR2P2F_4e<<LumiSect;
                    CR2P2F_4e<<" : ";
                    CR2P2F_4e<<Event<<endl;
                    
                    
                      PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_2P2F_4e], CR[Fake][_2P2F_4e], CR[BDfake][_2P2F_4e],CR[Conv][_2P2F_4e],false);
                    
                    
                    h1D_m4l_2P2F_4e->Fill(mass4l, weight);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                {
                    CR2P2F_4mu<<Run;
                    CR2P2F_4mu<<" : ";
                    CR2P2F_4mu<<LumiSect;
                    CR2P2F_4mu<<" : ";
                    CR2P2F_4mu<<Event<<endl;
                    
                       PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_2P2F_4mu], CR[Fake][_2P2F_4mu],CR[BDfake][_2P2F_4mu],CR[Conv][_2P2F_4mu],false);
                    h1D_m4l_2P2F_4mu->Fill(mass4l, weight);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                {
                    
                    CR2P2F_2e2mu<<Run;
                    CR2P2F_2e2mu<<" : ";
                    CR2P2F_2e2mu<<LumiSect;
                    CR2P2F_2e2mu<<" : ";
                    CR2P2F_2e2mu<<Event<<endl;
                    
                      PartOrigin(lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId,lep_Hindex,lep_id, mass4l,weight,CR[Prom][_2P2F_2e2mu], CR[Fake][_2P2F_2e2mu],CR[BDfake][_2P2F_2e2mu],CR[Conv][_2P2F_2e2mu],false);
                    
                     h1D_m4l_2P2F_2e2mu->Fill(mass4l, weight);
                    
                    if (abs(idL[2])+abs(idL[3])==22)  h1D_m4l_2P2F_2mu2e->Fill(mass4l, weight);
                    if (abs(idL[2])+abs(idL[3])==26)  h1D_m4l_2P2F_2e2mu->Fill(mass4l, weight);
                    
                }
                
                //  cout<<"Tezina = "<<weight<<" Mass4l = "<<mass4l<<endl;
                
                if (DR < 2*0.3) {
                    nEvt2P2FLeptonsIsoID_lt06 += 2;
                    nEvt2P2FLeptonsIDonly_lt06 += nFailedIDOnlyLeptonsZ2;
                }
                
            }
            if (nFailedLeptons == 3){
                nEvt1P3FLeptons++;
                //                h1D_m4l_2P2F->Fill(mass4l, weight);
            }
            if (nFailedLeptonsZ1) nEvtZ1WithFailedLeptons++;
            
        }
        if(!SILENT) {
            cout << fixed;
            cout.precision(printOutPrecision);
            cout << Event << ":" << Run << ":" << LumiSect << ":" << mass3l << ":" << massZ1 << ":" << weight << endl;
        }
    }
    
    //    if(!SILENT) {
    cout  << setw(printOutWidth) << "selected events in Z+1L region: " << nEvtPassedZ1LSelection << endl;
    cout  << setw(printOutWidth) << "selected events in ZX CRs:      " << nEvtPassedZXCRSelection << endl;
    cout  << setw(printOutWidth) << "selected events in 1P3F region: " << nEvt1P3FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 2P2F region: " << nEvt2P2FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 3P1F region: " << nEvt3P1FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 4P region:   " << nEvt4PLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in ZX CRs with failing Z1 leptons: " << nEvtZ1WithFailedLeptons << endl;
    cout  << setw(printOutWidth) << "2P2F that fail iso region:   " << 100.0*nEvt2P2FLeptonsIDonly_lt06/nEvt2P2FLeptonsIsoID_lt06 << endl;
    //    }
    
    
    
    cout<<fileName<<endl;
    cout<<"Electrons"<<endl;
    cout<<"Nprompt=  "<<n_prompt<<endl;
    cout<<"Nnonprompttau=  "<<n_nonprompttaus<<endl;
    cout<<"Nnonpromptconv=  "<<n_nonpromptconv<<endl;
    cout<<"Nnonpromptgamma= "<<n_nonpromptformphot<<endl;
    cout<<"NfakesBD= "<<n_bdfake<<endl;
    cout<<"Nfakes= "<<n_fakes<<endl;
    cout<<"Nothers=  "<<n_others<<endl;
    cout<<"Ntot= "<<n_pass<<endl;
    
    
    cout<<"Total skipped events = "<<SkippedEvents<<endl;
    
    
    return 0;
}


//_______________________________________________________________________________________________________________________________________________
double getFR(int lep_id, double lep_pt, double lep_eta, TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE){
    // if (lep_pt>80) cout << "lep_id: " << lep_id << ", lep_pt: " << lep_pt << ", lep_eta: " << lep_eta << endl;
    if ((fabs(lep_id) == 11) && (fabs(lep_eta) < 1.497)) return h1D_FRel_EB->GetBinContent(h1D_FRel_EB->FindBin(lep_pt));
    if ((fabs(lep_id) == 11) && (fabs(lep_eta) > 1.497)) return h1D_FRel_EE->GetBinContent(h1D_FRel_EE->FindBin(lep_pt));
    if ((fabs(lep_id) == 13) && (fabs(lep_eta) < 1.2)) return h1D_FRmu_EB->GetBinContent(h1D_FRmu_EB->FindBin(lep_pt));
    if ((fabs(lep_id) == 13) && (fabs(lep_eta) > 1.2)) return h1D_FRmu_EE->GetBinContent(h1D_FRmu_EE->FindBin(lep_pt));
    // cout << "bla, bla, bla " <<endl;
    return 0;
}

double getFRIso(int lep_id, double lep_pt, double lep_eta){
    //  cout << "lep_id: " << lep_id << ", lep_pt: " << lep_pt << ", lep_eta: " << lep_eta << endl;
    if ((fabs(lep_id) == 11) && (lep_pt < 10) )return 0.31;
    if ((fabs(lep_id) == 11) && ((lep_pt >= 10) && (lep_pt < 30))) return 0.332;
    if ((fabs(lep_id) == 11) && ((lep_pt >= 30) && (lep_pt < 60)))return 0.5;
    if ((fabs(lep_id) == 11) && ((lep_pt >= 60) )) return 0.64;
    
    if ((fabs(lep_id) == 13) && (lep_pt < 10)) return 0.31;
    if ((fabs(lep_id) == 13) && ((lep_pt >= 10) && (lep_pt < 30))) return 0.275;
    if ((fabs(lep_id) == 13) && ((lep_pt >= 30) && (lep_pt < 60))) return 0.35;
    if ((fabs(lep_id) == 13) && ((lep_pt >= 60))) return 0.65;
    // cout << "bla, bla, bla " <<endl;
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int getEstimatesFromCR(TTree* tree,
                       TH1D* h1D_FRel_EB,   TH1D* h1D_FRel_EE,   TH1D* h1D_FRmu_EB,   TH1D* h1D_FRmu_EE,
                       TH1D* &h1D_m4l_SR_2P2F, TH1D* &h1D_m4l_SR_2P2F_4e, TH1D* &h1D_m4l_SR_2P2F_4mu,TH1D* &h1D_m4l_SR_2P2F_2e2mu, TH1D* &h1D_m4l_SR_2P2F_2mu2e, TH1D* &h1D_m4l_SR_3P1F,TH1D* &h1D_m4l_SR_3P1F_4e,TH1D* &h1D_m4l_SR_3P1F_4mu,TH1D* &h1D_m4l_SR_3P1F_2e2mu, TH1D* &h1D_m4l_SR_3P1F_2mu2e,TH1D* &h1D_m4l_Add_2P2F, TH1D* &h1D_m4l_Add_2P2F_4e, TH1D* &h1D_m4l_Add_2P2F_4mu,TH1D* &h1D_m4l_Add_2P2F_2e2mu, TH1D* &h1D_m4l_Add_2P2F_2mu2e,
                       double ptElCut, double ptMuCut, double mZ2Cut){
    
    
    
    
    // define vars and branches
    float mass4l,massZ2,mass4lREFIT,massZ1,mass4lErr;
    Double_t Dbkg,Djet;
    float eventWeight, dataMCWeight, crossSection, genWeight;
    //    long int Run, LumiSect, Event;
    ULong64_t Run, LumiSect, Event;
    int finalState, nVtx;
    bool passedFullSelection, passedZ1LSelection, passedZ4lZXCRSelection, passedZXCRSelection, passedZ4lSelection;
    //    TClonesArray *lep_p4 = new TClonesArray("TLorentzVector", 10);
    vector<float> *lep_pt = 0; TBranch *b_lep_pt = 0;
    vector<float> *lep_eta = 0; TBranch *b_lep_eta = 0;
    vector<float> *lep_phi = 0; TBranch *b_lep_phi = 0;
    vector<float> *lep_mass = 0; TBranch *b_lep_mass = 0;
    int lep_Hindex[4];
    vector<int> *lep_id = 0; TBranch *b_lep_id = 0;
    vector<int> *lep_tightId = 0; TBranch *b_lep_tightId = 0;
    vector<float> *lep_RelIso = 0; TBranch *b_lep_RelIso = 0;
    //    vector<int> *lep_matchedR03_PdgId = 0; TBranch *b_lep_matchedR03_PdgId = 0;
    //    vector<int> *lep_matchedR03_MomId = 0; TBranch *b_lep_matchedR03_MomId = 0;
    //    vector<int> *lep_matchedR03_MomMomId = 0; TBranch *b_lep_matchedR03_MomMomId = 0;
    
    
    
    // counters
    int nEvtPassedZXCRSelection = 0;
    int nEvt2P2FLeptons = 0;
    int nEvt3P1FLeptons = 0;
    int nFailedLeptonsZ2 = 0;
    
    enum {tot,tot_2P2F,tot_3P1F};
    
    // get branches
    tree->SetBranchAddress("Run",&Run);
    tree->SetBranchAddress("LumiSect",&LumiSect);
    tree->SetBranchAddress("Event",&Event);
    tree->SetBranchAddress("crossSection",&crossSection);
    tree->SetBranchAddress("eventWeight",&eventWeight);
    tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
    tree->SetBranchAddress("genWeight",&genWeight);
    tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
    tree->SetBranchAddress("passedZ1LSelection",&passedZ1LSelection);
    tree->SetBranchAddress("passedZ4lZXCRSelection",&passedZ4lZXCRSelection);
    tree->SetBranchAddress("passedZXCRSelection",&passedZXCRSelection);
    tree->SetBranchAddress("passedZ4lSelection",&passedZ4lSelection);
    tree->SetBranchAddress("nVtx",&nVtx);
    tree->SetBranchAddress("finalState",&finalState);
    tree->SetBranchAddress("mass4l",&mass4l);
    tree->SetBranchAddress("mass4lREFIT",&mass4lREFIT);
    tree->SetBranchAddress("mass4lErr",&mass4lErr);
    //Dok Lucian ne popravi flag tree->SetBranchAddress("massZ1",&massZ1);
    
  //  tree->SetBranchAddress("lNEvents",&lNEvents);
    
    tree->SetBranchAddress("massZ2",&massZ2);
    tree->SetBranchAddress("D_bkg_kin",&Dbkg);
    tree->SetBranchAddress("Djet_VAJHU",&Djet);
    
    //    tree->SetBranchAddress("mass3l",&mass3l);
    tree->SetBranchAddress("lep_Hindex",&lep_Hindex);
    //    tree->GetBranch("lep_p4")->SetAutoDelete(kFALSE);
    //    tree->SetBranchAddress("lep_p4",&lep_p4);
    tree->SetBranchAddress("lep_pt",&lep_pt,&b_lep_pt);
    tree->SetBranchAddress("lep_eta",&lep_eta,&b_lep_eta);
    tree->SetBranchAddress("lep_phi",&lep_phi,&b_lep_phi);
    tree->SetBranchAddress("lep_mass",&lep_mass,&b_lep_mass);
    tree->SetBranchAddress("lep_id",&lep_id,&b_lep_id);
    tree->SetBranchAddress("lep_tightId",&lep_tightId,&b_lep_tightId);
    tree->SetBranchAddress("lep_RelIsoNoFSR",&lep_RelIso,&b_lep_RelIso);
    // fill histograms
    Long64_t nentries = tree->GetEntries();
    cout << "nentries: " << nentries << endl;
    
    // sort tree in "Event" number, ascending
    int *index = new int[nentries];
    if (SORT_EVENTS) {
        tree->Draw("Event","","goff");
        int nentries_int = static_cast<int>(nentries);
        TMath::Sort(nentries_int, tree->GetV1(), index, false);
    }
     for(int iEvt=0; iEvt < nentries; iEvt++){
        if (SORT_EVENTS) {
            tree->GetEntry(index[iEvt]);	//index[iEvt]);}// take sorted entries
        } else {
            tree->GetEntry(iEvt);
        }// take unsroted entries
        
        // weight
        float weight = eventWeight*dataMCWeight;
        
        if (isData) weight=1;
        
        if (!isData){
            // weight
            weight = crossSection*LUMI_INT/lNEvents;
        }
        
        //std::cout<<iEvt<<std::endl;
        
        if (isData) weight=1;
        //double Lumi=1547.18;
        if (!isData){
            if (fileName=="DY10")
                // weight
               // weight = genWeight*crossSection*LUMI_INT/lNEvents;
                weight = genWeight*18610.0*LUMI_INT/lNEvents;
            //weight = eventWeight*dataMCWeight*crossSection*LUMI_INT/lNEvents;
            if (fileName=="DY50")
                // weight
                weight = genWeight*6104.0*LUMI_INT/lNEvents;
            if (fileName=="TT")
                // weight
                weight = genWeight*815.96*LUMI_INT/lNEvents;
            if (fileName=="WZ")
                // weight
                weight = genWeight*4.430*LUMI_INT/lNEvents;
            
            if (fileName=="ZZ")
                // weight
                weight = genWeight*1.256*LUMI_INT/lNEvents;
        }
        
        bool passedCRsel = false;
        
        if (isData) passedCRsel =passedZXCRSelection;
        //if (!isData) passedCRsel = passedZ4lZXCRSelection;
        if (!isData) passedCRsel =passedZXCRSelection;       
 
        if (passedCRsel) {
            nEvtPassedZXCRSelection++;
            //        cout << "   nEvtPassedZ1LSelection: " << nEvtPassedZ1LSelection << "/" << iEvt << endl;
            //        cout << Event << ":" << Run << ":" << LumiSect << ":" << massZ1 << ":" << passedZ1LSelection << ":" << lep_Hindex << ":" << b_lep_id << ":" << lep_p4 << endl;
            //        cout << "iEvt: " << iEvt << ", lep_Hindex[2]: " << lep_Hindex[2] << ", lep_id->size():" << lep_id->size() << ", lep_tightId->size():" << lep_tightId->size() << endl;
            
            
            // get properties of the non-Z1 leptons (3rd, 4th)
            int lep_tight[4];
            float lep_iso[4];
            int idL[4];
            float pTL[4];
            float etaL[4];
            float phiL[4];
            float massL[4];
            TLorentzVector dilep;
            for(unsigned int k = 0; k <= 3; k++) {
                lep_tight[k] = lep_tightId->at(lep_Hindex[k]);
                lep_iso[k]= lep_RelIso->at(lep_Hindex[k]);
                idL[k] = lep_id->at(lep_Hindex[k]);
                //                TLorentzVector *lep = (TLorentzVector*) lep_p4->At(lep_Hindex[k]);
                TLorentzVector lep;
                lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[k]),lep_eta->at(lep_Hindex[k]),lep_phi->at(lep_Hindex[k]),lep_mass->at(lep_Hindex[k]));
                pTL[k]  = lep.Pt();
                etaL[k] = lep.Eta();
                
            }
            //Fix dok Lucian ne popravi flag
            TLorentzVector lep_1, lep_2;
            lep_1.SetPtEtaPhiM(lep_pt->at(lep_Hindex[0]),lep_eta->at(lep_Hindex[0]),lep_phi->at(lep_Hindex[0]),lep_mass->at(lep_Hindex[0]));
            lep_2.SetPtEtaPhiM(lep_pt->at(lep_Hindex[1]),lep_eta->at(lep_Hindex[1]),lep_phi->at(lep_Hindex[1]),lep_mass->at(lep_Hindex[1]));
            
            massZ1 = (lep_1+lep_2).M();
            
            TLorentzVector lep;
            lep.SetPtEtaPhiM(lep_pt->at(lep_Hindex[2]),lep_eta->at(lep_Hindex[2]),lep_phi->at(lep_Hindex[2]),lep_mass->at(lep_Hindex[2]));
            float pTL3  = lep.Pt();
            float etaL3 = lep.Eta();
            float phiL3 = lep.Phi();
            TLorentzVector lep2;
            lep2.SetPtEtaPhiM(lep_pt->at(lep_Hindex[3]),lep_eta->at(lep_Hindex[3]),lep_phi->at(lep_Hindex[3]),lep_mass->at(lep_Hindex[3]));
            float pTL4  = lep2.Pt();
            float etaL4 = lep2.Eta();
            float phiL4 = lep2.Phi();
            float DR = sqrt((phiL3-phiL4)*(phiL3-phiL4) + (etaL3-etaL4)*(etaL3-etaL4));
            
            double fr3 = getFR(idL[2], pTL[2], etaL[2], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            double fr4 = getFR(idL[3], pTL[3], etaL[3], h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE);
            double frShape;
            
            frShape= (fr3/(1-fr3)) * (fr4/(1-fr4));
            float Fill_value;
            
            //cout<<"ulaz"<<endl;
            
            for (int i = 0; i<=6;i++)
            {
                
                if (Shape_names[i] == "mass4lREFIT") Fill_value = mass4lREFIT;
                if (Shape_names[i] == "Dbkg") Fill_value = Dbkg;
                if (Shape_names[i] == "massZ1") Fill_value = massZ1;
                if (Shape_names[i] == "massZ2") Fill_value = massZ2;
                if (Shape_names[i] == "Djet") Fill_value = Djet;
                if (Shape_names[i] == "mass4lErr")
                {Fill_value = 1.0*mass4lErr/mass4l; //cout<<Fill_value<<"   " <<weight*frShape<<endl; cout<<SR_shapes[i]->Integral()<<endl;
                }
                
                if (Shape_names[i] == "mass4lErrREFIT")
                {Fill_value = 1.0*mass4lErrREFIT/mass4lREFIT; //cout<<Fill_value<<"   " <<weight*frShape<<endl; cout<<SR_shapes[i]->Integral()<<endl;
                }
                
                
                // cout<<i<<"  "<<Shape_names[i]<<" "<< Fill_value<<" "<<weight* frShape<<endl;
                
                SR_shapes[i] ->Fill(Fill_value,weight*frShape);
                
                
                
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44)
                {
                    SR_shapes_4e[i]->Fill(Fill_value, weight* frShape);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                {
                    SR_shapes_4mu[i]->Fill(Fill_value, weight* frShape);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                {
                    
                    SR_shapes_2e2mu[i]->Fill(Fill_value, weight* frShape);
                }
                
                
                if ((mass4l > 118) && (mass4l < 130))
                {
                    
                    CSR_shapes[i] ->Fill(Fill_value,weight*frShape);
                    
                    
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44)
                    {
                        CSR_shapes_4e[i]->Fill(Fill_value, weight* frShape);
                    }
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                    {
                        CSR_shapes_4mu[i]->Fill(Fill_value, weight* frShape);
                    }
                    if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                    {
                        
                        CSR_shapes_2e2mu[i]->Fill(Fill_value, weight* frShape);
                    }
                    
                    
                    
                    
                }
                
            }
            
            
            
            D2plotmass4lDbkg->Fill(mass4l,Dbkg,weight* frShape);
            // cout<<"izlaz"<<endl;
            
            
            
            
            // count the failed leptons
            //            nFailedLeptonsZ1 = !(lep_tight[0] && ((abs(idL[0])==11 && lep_iso[0]<0.5) || (abs(idL[0])==13 && lep_iso[0]<0.4))) +
            //                                   !(lep_tight[1] && ((abs(idL[1])==11 && lep_iso[1]<0.5) || (abs(idL[1])==13 && lep_iso[1]<0.4)));
            nFailedLeptonsZ2 = !(lep_tight[2] && ((abs(idL[2])==11) || (abs(idL[2])==13 && lep_iso[2]<0.35))) +
            !(lep_tight[3] && ((abs(idL[3])==11) || (abs(idL[3])==13 && lep_iso[3]<0.35)));
            //            nFailedLeptons   = nFailedLeptonsZ1 + nFailedLeptonsZ2;
            
            // prepare the estimates
            if (nFailedLeptonsZ2 == 1) {
                nEvt3P1FLeptons++;
                
                // if(!isData) weight*= (1.37/1.23);
                
                // double fr3=0.06;
                //double fr4=0.06;
                
                double fr = (!(lep_tight[2] && ((abs(idL[2])==11) || (abs(idL[2])==13 && lep_iso[2]<0.35))))*(fr3/(1-fr3)) +
                (!(lep_tight[3] && ((abs(idL[3])==11) || (abs(idL[3])==13 && lep_iso[3]<0.35))))*(fr4/(1-fr4));
                double frX;
                if (!(lep_tight[2] && ((abs(idL[2])==11) || (abs(idL[2])==13 && lep_iso[2]<0.35)))) frX = (fr3/(1-fr3));
                if (!(lep_tight[3] && ((abs(idL[3])==11) || (abs(idL[3])==13 && lep_iso[3]<0.35)))) frX = (fr4/(1-fr4));
                // if (DR < 0.6) fr = 0;
                
                h1D_m4l_SR_3P1F->Fill(mass4l, weight * fr);
                SR_mass4Refit[tot_3P1F]->Fill(mass4lREFIT, weight * fr);
                
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44)
                {
                    h1D_m4l_SR_3P1F_4e->Fill(mass4l, weight* fr);
                    SR_mass4Refit_4e[tot_3P1F]->Fill(mass4lREFIT, weight * fr);
                    
                    
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                {
                    h1D_m4l_SR_3P1F_4mu->Fill(mass4l, weight* fr);
                    SR_mass4Refit_4mu[tot_3P1F]->Fill(mass4lREFIT, weight * fr);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                {
                    
                    //h1D_m4l_SR_3P1F_2e2mu->Fill(mass4l, weight* fr);
                    SR_mass4Refit_2e2mu[tot_3P1F]->Fill(mass4lREFIT, weight * fr);
                    
                    if ((abs(idL[2])+abs(idL[3]))== 26) h1D_m4l_SR_3P1F_2e2mu->Fill(mass4l, weight* fr);
                    if ((abs(idL[2])+abs(idL[3]))== 22) h1D_m4l_SR_3P1F_2mu2e->Fill(mass4l, weight* fr);
                    
                }
                //}
                //                cout<<"mass 4l = "<< mass4l<< "  ; weight = "<<weight<<"  ; fr3 = "<<fr3<<endl;
            }
            if (nFailedLeptonsZ2 == 2){
                
                //  if(!isData) weight*=(1.53/1.23);
                
                nEvt2P2FLeptons++;
                //  double fr3=0.06;
                //double fr4=0.06;
                double fr;
                
                fr= (fr3/(1-fr3)) * (fr4/(1-fr4));
                
                double fr2 =(fr3/(1-fr3)) + (fr4/(1-fr4));
                
                
                h1D_m4l_Add_2P2F->Fill(mass4l, weight * fr2);
                
                /* if (DR < 0.6)
                 {
                 double Riso = 0.3;
                 
                 double fr3iso,fr4iso,fr3isoCor,fr4isoCor;
                 fr3iso = getFRIso(idL[2], pTL[2], etaL[2]);
                 fr4iso = getFRIso(idL[3], pTL[3], etaL[3]);
                 
                 double num_Cor_iso =fr3*fr4/(fr3iso*fr4iso)*(fr3iso*fr4iso*(DR/(2*Riso))*(DR/(2*Riso)) + sqrt( fr3iso*fr4iso )*(1 - (DR/(2*Riso))));
                 double denom_Cor_iso =(1-fr3)*(1-fr4)/((1-fr3iso)*(1-fr4iso))*((1-fr3iso)*(1-fr4iso)*(DR/(2*Riso))*(DR/(2*Riso)) + sqrt( (1-fr3iso)*(1-fr4iso) )*(1 - (DR/(2*Riso))));
                 
                 fr = num_Cor_iso/denom_Cor_iso;
                 cout <<"FR = "<< fr<<endl;
                 }*/
                
                /*else {
                 double Riso = 0.3;
                 double  num = fr3*fr4*(DR/(2*Riso))*(DR/(2*Riso)) + sqrt( fr3*fr4 )*(1 - (DR/(2*Riso)));
                 double  denom = (1-fr3)*(1-fr4)*(DR/(2*Riso))* (DR/(2*Riso)) +  sqrt( (1-fr3)*(1-fr4) )*( 1 - (DR/(2*Riso)));
                 fr = num/denom;
                 //}
                 // */
                h1D_m4l_SR_2P2F->Fill(mass4l, weight * fr);
                SR_mass4Refit[tot_2P2F]->Fill(mass4lREFIT, weight * fr);
                
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==44)
                {
                    h1D_m4l_SR_2P2F_4e->Fill(mass4l, weight* fr);
                    h1D_m4l_Add_2P2F_4e->Fill(mass4l, weight * fr2);
                    
                    SR_mass4Refit_4e[tot_2P2F]->Fill(mass4lREFIT, weight * fr);
                    
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==52)
                {
                    h1D_m4l_SR_2P2F_4mu->Fill(mass4l, weight* fr);
                    h1D_m4l_Add_2P2F_4mu->Fill(mass4l, weight * fr2);
                    
                    SR_mass4Refit_4mu[tot_2P2F]->Fill(mass4lREFIT, weight * fr);
                }
                if ((abs(idL[0])+abs(idL[1])+abs(idL[2])+abs(idL[3]))==48)
                {
                    //h1D_m4l_SR_2P2F_2e2mu->Fill(mass4l, weight* fr);
                    //h1D_m4l_Add_2P2F_2e2mu->Fill(mass4l, weight * fr2);
                    
                    
                    if ((abs(idL[2])+abs(idL[3]))== 26)
                    {
                        h1D_m4l_SR_2P2F_2e2mu->Fill(mass4l, weight* fr);
                        h1D_m4l_Add_2P2F_2e2mu->Fill(mass4l, weight * fr2);
                        
                    }
                    if ((abs(idL[2])+abs(idL[3]))== 22)
                    {
                        h1D_m4l_SR_2P2F_2mu2e->Fill(mass4l, weight* fr);
                        h1D_m4l_Add_2P2F_2mu2e->Fill(mass4l, weight * fr2);
                    }
                    
                    SR_mass4Refit_2e2mu[tot_2P2F]->Fill(mass4lREFIT, weight * fr);
                }
                // }
                //cout<<"mass 4l = "<< mass4l<< "  ; weight = "<<weight<<"  ; fr3 = "<<fr3<<endl;
                
            }
        }
    }
    
    
    
    cout  << setw(printOutWidth) << "selected events in ZX CRs:      " << nEvtPassedZXCRSelection << endl;
    cout  << setw(printOutWidth) << "selected events in 2P2F region: " << nEvt2P2FLeptons << endl;
    cout  << setw(printOutWidth) << "selected events in 3P1F region: " << nEvt3P1FLeptons << endl;
    
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
void estimateZX(TString slimmedZXFileName,TString DataFileName, double ptElCut, double ptMuCut, double mZ2Cut){
    
    // common properties
    Width_t lineWidth = 2;
    double leg_xl = 0.50, leg_xr = 0.90, leg_yb = 0.72, leg_yt = 0.90;
    
    // get the FR histograms and slimmed ZX tree
    TFile *f1 = TFile::Open(slimmedZXFileName);
    TTree* zxTree1 = (TTree*) f1->Get("passedEvents");
    if (!isData)
        zxTree1 = new TChain("Ana/passedEvents");
    TH1D* h1D_FRel_EB = (TH1D*) f1->Get("Data_FRel_EB");
    TH1D* h1D_FRel_EE = (TH1D*) f1->Get("Data_FRel_EE");
    TH1D* h1D_FRmu_EB = (TH1D*) f1->Get("Data_FRmu_EB");
    TH1D* h1D_FRmu_EE = (TH1D*) f1->Get("Data_FRmu_EE");
    
    TFile *f = TFile::Open(DataFileName);
    TTree* zxTree = (TTree*) f->Get("passedEvents");
    
    
    // define dummy histogram for CRs
    double var_plotHigh = 2000.0; double var_plotLow = 0.0; double var_nBins = 100;// povecan broj bin-ova sa 100->200
    TH1D* h1D_dummyCR = new TH1D("dummyCR", "dummyCR", var_nBins, var_plotLow, var_plotHigh);
    TString varAxLabel = "m_{4l} [GeV]";
    double binWidth = ((int) (100*(var_plotHigh - var_plotLow)/var_nBins))/100.;
    TString sUnit = (varAxLabel.Contains(" [GeV]"))?" GeV":" ";
    TString sBinWidth = TString::Format("%.1f",binWidth) + sUnit;
    setHistProperties(h1D_dummyCR,1,1,kBlack,0,0,varAxLabel,"Events/"+sBinWidth);
    
    // define CR histograms
    TH1D* h1D_m4l_SR_2P2F = new TH1D("h1D_m4l_SR_2P2F","h1D_m4l_SR_2P2F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F->Sumw2();
    TH1D* h1D_m4l_SR_3P1F = new TH1D("h1D_m4l_SR_3P1F","h1D_m4l_SR_3P1F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F->Sumw2();
    TH1D* h1D_m4l_Add_2P2F = new TH1D("h1D_m4l_Add_2P2F","h1D_m4l_Add_2P2F",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_Add_2P2F->Sumw2();
    
    // define CR histograms_4e
    TH1D* h1D_m4l_SR_2P2F_4e = new TH1D("h1D_m4l_SR_2P2F_4e","h1D_m4l_SR_2P2F_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F_4e->Sumw2();
    TH1D* h1D_m4l_SR_3P1F_4e= new TH1D("h1D_m4l_SR_3P1F_4e","h1D_m4l_SR_3P1F_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F_4e->Sumw2();
    TH1D* h1D_m4l_Add_2P2F_4e = new TH1D("h1D_m4l_Add_2P2F_4e","h1D_m4l_Add_2P2F_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_Add_2P2F_4e->Sumw2();
    
    // define CR histograms_4mu
    TH1D* h1D_m4l_SR_2P2F_4mu = new TH1D("h1D_m4l_SR_2P2F_4mu","h1D_m4l_SR_2P2F_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F_4mu->Sumw2();
    TH1D* h1D_m4l_SR_3P1F_4mu = new TH1D("h1D_m4l_SR_3P1F_4mu","h1D_m4l_SR_3P1F_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F_4mu->Sumw2();
    TH1D* h1D_m4l_Add_2P2F_4mu = new TH1D("h1D_m4l_Add_2P2F_4mu","h1D_m4l_Add_2P2F_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_Add_2P2F_4mu->Sumw2();
    
    
    // define CR histograms_2e2mu
    TH1D* h1D_m4l_SR_2P2F_2e2mu = new TH1D("h1D_m4l_SR_2P2F_2e2mu","h1D_m4l_SR_2P2F_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F_2e2mu->Sumw2();
    TH1D* h1D_m4l_SR_3P1F_2e2mu = new TH1D("h1D_m4l_SR_3P1F_2e2mu","h1D_m4l_SR_3P1F_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F_2e2mu->Sumw2();
    TH1D* h1D_m4l_Add_2P2F_2e2mu = new TH1D("h1D_m4l_Add_2P2F_2e2mu","h1D_m4l_Add_2P2F_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_Add_2P2F_2e2mu->Sumw2();
    
    // define CR histograms_2mu2e
    TH1D* h1D_m4l_SR_2P2F_2mu2e = new TH1D("h1D_m4l_SR_2P2F_2mu2e","h1D_m4l_SR_2P2F_2mu2e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_2P2F_2mu2e->Sumw2();
    TH1D* h1D_m4l_SR_3P1F_2mu2e = new TH1D("h1D_m4l_SR_3P1F_2mu2e","h1D_m4l_SR_3P1F_2mu2e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_3P1F_2mu2e->Sumw2();
    TH1D* h1D_m4l_Add_2P2F_2mu2e = new TH1D("h1D_m4l_Add_2P2F_2e2mu","h1D_m4l_Add_2P2F_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_Add_2P2F_2mu2e->Sumw2();
    
    enum {tot,tot_2P2F,tot_3P1F};
    TString nm[] = {"tot","2P2F","3P1F"};
    
    for (int i = tot; i<=tot_3P1F ; i++)
    {
        SR_mass4Refit[i] = new  TH1D("m4lREFIT_"+nm[i],"m4lREFIT_"+nm[i],var_nBins,var_plotLow,var_plotHigh); SR_mass4Refit[i]->Sumw2();
        SR_mass4Refit_4e[i] = new  TH1D("m4lREFIT_4e_"+nm[i],"m4lREFIT_4e_"+nm[i],var_nBins,var_plotLow,var_plotHigh);  SR_mass4Refit_4e[i]->Sumw2();
        SR_mass4Refit_4mu[i] = new  TH1D("m4lREFIT_4mu_"+nm[i],"m4lREFIT_4mu_"+nm[i],var_nBins,var_plotLow,var_plotHigh); SR_mass4Refit_4mu[i]->Sumw2();
        SR_mass4Refit_2e2mu[i] = new  TH1D("m4lREFIT_2e2mu_"+nm[i],"m4lREFIT_2e2mu_"+nm[i],var_nBins,var_plotLow,var_plotHigh); SR_mass4Refit_2e2mu[i]->Sumw2();
        
    }
    
    
    for (int i = mass4lREFIT; i<=mass4lErrREFIT; i++)
    {
        SR_shapes[i] = new TH1D("SR_"+Shape_names[i],"SR_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); SR_shapes[i]->Sumw2();
        SR_shapes_4e[i] = new TH1D("SR_4e_"+Shape_names[i],"SR_4e_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); SR_shapes_4e[i]->Sumw2();
        SR_shapes_4mu[i] = new TH1D("SR_4mu_"+Shape_names[i],"SR_4mu_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); SR_shapes_4mu[i]->Sumw2();
        SR_shapes_2e2mu[i] = new TH1D("SR_2e2mu_"+Shape_names[i],"SR_2e2mu_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); SR_shapes_2e2mu[i]->Sumw2();
        
    }
    
    for (int i = mass4lREFIT; i<=mass4lErrREFIT; i++)
    {
        CSR_shapes[i] = new TH1D("CSR_"+Shape_names[i],"CSR_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); CSR_shapes[i]->Sumw2();
        CSR_shapes_4e[i] = new TH1D("CSR_4e_"+Shape_names[i],"CSR_4e_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); CSR_shapes_4e[i]->Sumw2();
        CSR_shapes_4mu[i] = new TH1D("CSR_4mu_"+Shape_names[i],"CSR_4mu_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); CSR_shapes_4mu[i]->Sumw2();
        CSR_shapes_2e2mu[i] = new TH1D("CSR_2e2mu_"+Shape_names[i],"CSR_2e2mu_"+Shape_names[i],CR_nBins[i], CR_min[i], CR_max[i]); CSR_shapes_2e2mu[i]->Sumw2();
        
    }
    
    
    
    D2plotmass4lDbkg = new TH2D("h2D_mass4l_Dbkg","h2D_mass4l_Dbkg",D2plotmass4lDbkg_mass4l_nBin, D2plotmass4lDbkg_mass4l_min, D2plotmass4lDbkg_mass4l_max,D2plotmass4lDbkg_Dbkg_nBin, D2plotmass4lDbkg_Dbkg_min, D2plotmass4lDbkg_Dbkg_max);
    
    D2plotmass4lDbkg->Sumw2();
    
    // get histograms from tree AnaZX/passedEvents
    getEstimatesFromCR(zxTree, h1D_FRel_EB, h1D_FRel_EE, h1D_FRmu_EB, h1D_FRmu_EE, h1D_m4l_SR_2P2F,h1D_m4l_SR_2P2F_4e,h1D_m4l_SR_2P2F_4mu,h1D_m4l_SR_2P2F_2e2mu, h1D_m4l_SR_2P2F_2mu2e, h1D_m4l_SR_3P1F,h1D_m4l_SR_3P1F_4e,h1D_m4l_SR_3P1F_4mu,h1D_m4l_SR_3P1F_2e2mu, h1D_m4l_SR_3P1F_2mu2e, h1D_m4l_Add_2P2F,h1D_m4l_Add_2P2F_4e,h1D_m4l_Add_2P2F_4mu,h1D_m4l_Add_2P2F_2e2mu, h1D_m4l_Add_2P2F_2mu2e, ptElCut, ptMuCut, mZ2Cut);
    
    
    /* cout<<"----MASS4L REFIT---------------------------"<<endl;
     
     cout<<"h1D_m4l_SR_2P2F = "<<SR_mass4Refit[tot_2P2F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_4e = "<<SR_mass4Refit_4e[tot_2P2F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_4mu = "<<SR_mass4Refit_4mu[tot_2P2F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_2e2mu = "<<SR_mass4Refit_2e2mu[tot_2P2F]->Integral()<<endl;
     cout<<"--------------------------------"<<endl;
     cout<<"h1D_m4l_SR_3P1F = "<<SR_mass4Refit[tot_3P1F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_4e = "<<SR_mass4Refit_4e[tot_3P1F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_4mu = "<<SR_mass4Refit_4mu[tot_3P1F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_2e2mu = "<<SR_mass4Refit_2e2mu[tot_3P1F]->Integral()<<endl;
     cout<<"--------------Mass4l-----------------"<<endl;
     
     cout<<"h1D_m4l_SR_2P2F = "<<h1D_m4l_SR_2P2F->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_4e = "<<h1D_m4l_SR_2P2F_4e->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_4mu = "<<h1D_m4l_SR_2P2F_4mu->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_2e2mu = "<<h1D_m4l_SR_2P2F_2e2mu->Integral()<<endl;
     cout<<"--------------------------------"<<endl;
     cout<<"h1D_m4l_SR_3P1F = "<<h1D_m4l_SR_3P1F->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_4e = "<<h1D_m4l_SR_3P1F_4e->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_4mu = "<<h1D_m4l_SR_3P1F_4mu->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_2e2mu = "<<h1D_m4l_SR_3P1F_2e2mu->Integral()<<endl;
     cout<<"--------------------------------"<<endl;*/
    
    
    
    
    /* if (isData){
     TFile *ZZ = new TFile ("plotsZX/estimatesZXZZ.root","OPEN");
     
     TH1D *ZZ_4e = (TH1D*) ZZ->Get("h1D_m4l_SR_3P1F_4e");
     TH1D *ZZ_4mu = (TH1D*) ZZ->Get("h1D_m4l_SR_3P1F_4mu");
     TH1D *ZZ_2e2mu = (TH1D*) ZZ->Get("h1D_m4l_SR_3P1F_2e2mu");
     TH1D *ZZ_2mu2e = (TH1D*) ZZ->Get("h1D_m4l_SR_3P1F_2mu2e");
     TH1D *ZZ_tot = (TH1D*) ZZ->Get("h1D_m4l_SR_3P1F");
     
     
     h1D_m4l_SR_3P1F->Add(ZZ_tot,-1);
     h1D_m4l_SR_3P1F_4e->Add(ZZ_4e,-1);
     h1D_m4l_SR_3P1F_4mu->Add(ZZ_4mu,-1);
     h1D_m4l_SR_3P1F_2e2mu->Add(ZZ_2e2mu,-1);
     h1D_m4l_SR_3P1F_2mu2e->Add(ZZ_2mu2e,-1);
     
     
     TH1D *RZZ_4e = (TH1D*) ZZ->Get("m4lREFIT_4e_3P1F");
     TH1D *RZZ_4mu = (TH1D*) ZZ->Get("m4lREFIT_4mu_3P1F");
     TH1D *RZZ_2e2mu = (TH1D*) ZZ->Get("m4lREFIT_2e2mu_3P1F");
     TH1D *RZZ_tot = (TH1D*) ZZ->Get("m4lREFIT_3P1F");
     
     
     SR_mass4Refit[tot_3P1F]->Add(RZZ_tot,-1);
     SR_mass4Refit_4e[tot_3P1F]->Add(RZZ_4e,-1);
     SR_mass4Refit_4mu[tot_3P1F]->Add(RZZ_4mu,-1);
     SR_mass4Refit_2e2mu[tot_3P1F]->Add(RZZ_2e2mu,-1);
     
     }*/
    
    TH1D* h1D_m4l_SR_tot = new TH1D("h1D_m4l_SR_tot","h1D_m4l_SR_tot",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_tot->Sumw2();
    TH1D* h1D_m4l_SR_tot_4e = new TH1D("h1D_m4l_SR_tot_4e","h1D_m4l_SR_tot_4e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_tot_4e->Sumw2();
    TH1D* h1D_m4l_SR_tot_4mu = new TH1D("h1D_m4l_SR_tot_4mu","h1D_m4l_SR_tot_4mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_tot_4mu->Sumw2();
    TH1D* h1D_m4l_SR_tot_2e2mu = new TH1D("h1D_m4l_SR_tot_2e2mu","h1D_m4l_SR_tot_2e2mu",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_tot_2e2mu->Sumw2();
    TH1D* h1D_m4l_SR_tot_2mu2e = new TH1D("h1D_m4l_SR_tot_2mu2e","h1D_m4l_SR_tot_2mu2e",var_nBins, var_plotLow, var_plotHigh); h1D_m4l_SR_tot_2mu2e->Sumw2();
    
    for (int c = 1; c<=var_nBins; c++)
    {
        
         if (h1D_m4l_SR_3P1F_4e->GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_4e->GetBinContent(c)) {
         
         h1D_m4l_SR_tot_4e->SetBinContent(c,h1D_m4l_SR_2P2F_4e->GetBinContent(c));
         h1D_m4l_SR_tot_4e->SetBinError(c,h1D_m4l_SR_2P2F_4e->GetBinError(c));
         
         }
         
         else{
         h1D_m4l_SR_tot_4e->SetBinContent(c,h1D_m4l_SR_3P1F_4e->GetBinContent(c)-h1D_m4l_SR_2P2F_4e->GetBinContent(c));
         h1D_m4l_SR_tot_4e->SetBinError(c,sqrt(h1D_m4l_SR_2P2F_4e->GetBinError(c)*h1D_m4l_SR_2P2F_4e->GetBinError(c)+h1D_m4l_SR_3P1F_4e->GetBinError(c)*h1D_m4l_SR_3P1F_4e->GetBinError(c)));
         
         }
        
        
        if (h1D_m4l_SR_3P1F_4mu->GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_4mu->GetBinContent(c)) {
            
            h1D_m4l_SR_tot_4mu->SetBinContent(c,h1D_m4l_SR_2P2F_4mu->GetBinContent(c));
            h1D_m4l_SR_tot_4mu->SetBinError(c,h1D_m4l_SR_2P2F_4mu->GetBinError(c));
            
        }
        
        else{
            h1D_m4l_SR_tot_4mu->SetBinContent(c,h1D_m4l_SR_3P1F_4mu->GetBinContent(c)-h1D_m4l_SR_2P2F_4mu->GetBinContent(c));
            h1D_m4l_SR_tot_4mu->SetBinError(c,sqrt(h1D_m4l_SR_2P2F_4mu->GetBinError(c)*h1D_m4l_SR_2P2F_4mu->GetBinError(c)+h1D_m4l_SR_3P1F_4mu->GetBinError(c)*h1D_m4l_SR_3P1F_4mu->GetBinError(c)));
            
        }
        
        
        
        if (h1D_m4l_SR_3P1F_2e2mu->GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_2e2mu->GetBinContent(c)) {
            
            h1D_m4l_SR_tot_2e2mu->SetBinContent(c,h1D_m4l_SR_2P2F_2e2mu->GetBinContent(c));
            h1D_m4l_SR_tot_2e2mu->SetBinError(c,h1D_m4l_SR_2P2F_2e2mu->GetBinError(c));
            
        }
        
        else{
            h1D_m4l_SR_tot_2e2mu->SetBinContent(c,h1D_m4l_SR_3P1F_2e2mu->GetBinContent(c)-h1D_m4l_SR_2P2F_2e2mu->GetBinContent(c));
            h1D_m4l_SR_tot_2e2mu->SetBinError(c,sqrt(h1D_m4l_SR_2P2F_2e2mu->GetBinError(c)*h1D_m4l_SR_2P2F_2e2mu->GetBinError(c)+h1D_m4l_SR_3P1F_2e2mu->GetBinError(c)*h1D_m4l_SR_3P1F_2e2mu->GetBinError(c)));
            
        }
        
        
           if (h1D_m4l_SR_3P1F_2mu2e->GetBinContent(c) <= 2*h1D_m4l_SR_2P2F_2mu2e->GetBinContent(c)) {
         
         h1D_m4l_SR_tot_2mu2e->SetBinContent(c,h1D_m4l_SR_2P2F_2mu2e->GetBinContent(c));
         h1D_m4l_SR_tot_2mu2e->SetBinError(c,h1D_m4l_SR_2P2F_2mu2e->GetBinError(c));
         
         }
         
         else{
         h1D_m4l_SR_tot_2mu2e->SetBinContent(c,h1D_m4l_SR_3P1F_2mu2e->GetBinContent(c)-h1D_m4l_SR_2P2F_2mu2e->GetBinContent(c));
         h1D_m4l_SR_tot_2mu2e->SetBinError(c,sqrt(h1D_m4l_SR_2P2F_2mu2e->GetBinError(c)*h1D_m4l_SR_2P2F_2mu2e->GetBinError(c)+h1D_m4l_SR_3P1F_2mu2e->GetBinError(c)*h1D_m4l_SR_3P1F_2mu2e->GetBinError(c)));
         
         }
        
        
        
    }
    
    cout<<" Provera " << "4e tot 3p1f = "<<SR_mass4Refit_4e[tot_3P1F]->Integral()<< " 4e tot 2p2f = "<<SR_mass4Refit_4e[tot_2P2F]->Integral()<<endl;
    
    double a2 = 0;double a3 = 0;
    
    for (int c = 1; c<=var_nBins; c++)
    {
        if (SR_mass4Refit_4e[tot_3P1F]->GetBinContent(c) <= 2*SR_mass4Refit_4e[tot_2P2F]->GetBinContent(c))
        {
            
            SR_mass4Refit_4e[tot]->SetBinContent(c,SR_mass4Refit_4e[tot_2P2F]->GetBinContent(c));
            SR_mass4Refit_4e[tot]->SetBinError(c,SR_mass4Refit_4e[tot_2P2F]->GetBinError(c));
            
        }
        
        else{
            SR_mass4Refit_4e[tot]->SetBinContent(c,SR_mass4Refit_4e[tot_3P1F]->GetBinContent(c)-SR_mass4Refit_4e[tot_2P2F]->GetBinContent(c));
            SR_mass4Refit_4e[tot]->SetBinError(c,sqrt(SR_mass4Refit_4e[tot_3P1F]->GetBinError(c)*SR_mass4Refit_4e[tot_3P1F]->GetBinError(c)+SR_mass4Refit_4e[tot_2P2F]->GetBinError(c)*SR_mass4Refit_4e[tot_2P2F]->GetBinError(c)));
            
        }
        
        
        
        
    }
    
    
    cout<< SR_mass4Refit_4e[tot]->Integral()<<endl;
    
    
    
    
    
    
    
    // total contribution (3P1F - 2P2F) 4e
   // TH1D* h1D_m4l_SR_tot_4e = (TH1D*) h1D_m4l_SR_3P1F_4e->Clone("h1D_m4l_SR_tot_4e");
    //h1D_m4l_SR_tot_4e->Add(h1D_m4l_SR_2P2F_4e, -1);
    
    
    // total contribution (3P1F - 2P2F) 4mu
    /* TH1D* h1D_m4l_SR_tot_4mu = (TH1D*) h1D_m4l_SR_3P1F_4mu->Clone("h1D_m4l_SR_tot_4mu");
     
     h1D_m4l_SR_tot_4mu->Add(h1D_m4l_SR_2P2F_4mu, -1);
     
     SR_mass4Refit_4mu[tot] = (TH1D*) SR_mass4Refit_4mu[tot_2P2F]->Clone("m4lREFIT_4mu_tot");
     SR_mass4Refit_4mu[tot]->SetTitle("m4lREFIT_4mu_tot");
     
     // total contribution (3P1F - 2P2F) 2e2mu
     TH1D* h1D_m4l_SR_tot_2e2mu = (TH1D*) h1D_m4l_SR_3P1F_2e2mu->Clone("h1D_m4l_SR_tot_2e2mu");
     
     h1D_m4l_SR_tot_2e2mu->Add(h1D_m4l_SR_2P2F_2e2mu, -1);*/
    
    
    // total contribution (3P1F - 2P2F) 2mu2e
    //TH1D* h1D_m4l_SR_tot_2mu2e = (TH1D*) h1D_m4l_SR_3P1F_2mu2e->Clone("h1D_m4l_SR_tot_2mu2e");
    
   //h1D_m4l_SR_tot_2mu2e->Add(h1D_m4l_SR_2P2F_2mu2e, -1);
    
    //Total
    
    h1D_m4l_SR_tot->Add(h1D_m4l_SR_tot_4e);
    h1D_m4l_SR_tot->Add(h1D_m4l_SR_tot_4mu);
    h1D_m4l_SR_tot->Add(h1D_m4l_SR_tot_2e2mu);
    h1D_m4l_SR_tot->Add(h1D_m4l_SR_tot_2mu2e);
    
    SR_mass4Refit_2e2mu[tot] = (TH1D*) SR_mass4Refit_2e2mu[tot_2P2F]->Clone("m4lREFIT_2e2mu_tot");
    SR_mass4Refit_2e2mu[tot]->SetTitle("m4lREFIT_2e2mu_tot");
    
    // total contribution (3P1F - 2P2F)
    /*  TH1D* h1D_m4l_SR_tot = (TH1D*) h1D_m4l_SR_tot_4e->Clone("h1D_m4l_SR_tot");
     h1D_m4l_SR_tot->Add(h1D_m4l_SR_tot_4mu, 1);
     h1D_m4l_SR_tot->Add(h1D_m4l_SR_tot_2e2mu, 1);*/
    
    
    SR_mass4Refit[tot] = (TH1D*)  SR_mass4Refit_4e[tot] -> Clone("m4lREFIT_tot");
    SR_mass4Refit[tot]->Add( SR_mass4Refit_4mu[tot],1);
    SR_mass4Refit[tot]->Add( SR_mass4Refit_2e2mu[tot],1);
    SR_mass4Refit[tot]->SetTitle("m4lREFIT_tot");
    if (isData){
        
        // plot the histograms
        // setup environment & canvas
        TCanvas *c1;
        setCavasAndStyles("c1",c1,"");
        TLegend* leg = 0;
        
        TH1D *Shape, *Shape2,*Shape_t;
        
        
        // Estimates from Control Regions
        // SR from 2P2F
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_2P2F->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 2P+2F");
        setHistProperties(h1D_m4l_SR_2P2F,lineWidth,1,kRed-7); h1D_m4l_SR_2P2F->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_2P2F, "Estimated from data","L");
        Shape2 = (TH1D*) h1D_m4l_SR_2P2F->Clone("Shape2");
        Shape2->Smooth();
        setHistProperties(Shape2,lineWidth,1,kGreen-7);//leg->AddEntry(Shape2, "Smooth","L");
        //Shape2->Draw("hist same");
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_2P2F.pdf");
        
        
        
        
        
        
        
        // SR form 3P1F
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_3P1F->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 3P+1F");
        setHistProperties(h1D_m4l_SR_3P1F,lineWidth,1,kRed-7); h1D_m4l_SR_3P1F->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_3P1F, "Data","L");
        Shape = (TH1D*) h1D_m4l_SR_3P1F->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same goff");
        Shape2->Add(Shape2);//leg->AddEntry(Shape2, "Smooth 2*2P2F","L");
        //Shape2->Draw("hist same");
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_3P1F.pdf");
        
        
        
        // Estimates from Control Regions
        // SR from 2P2F 4e
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_2P2F_4e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 2P+2F 4e");
        setHistProperties(h1D_m4l_SR_2P2F_4e,lineWidth,1,kRed-7); h1D_m4l_SR_2P2F_4e->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_2P2F_4e, "Estimated from data","L");
        Shape2 = (TH1D*) h1D_m4l_SR_2P2F_4e->Clone("Shape2");
        Shape2->Smooth();
        setHistProperties(Shape2,lineWidth,1,kGreen-7);//leg->AddEntry(Shape2, "Smooth","L");
        //Shape2->Draw("hist same");
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_2P2F_4e.pdf");
        
        
        
        
        
        // SR form 3P1F_4e
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_3P1F_4e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 3P+1F 4e");
        setHistProperties(h1D_m4l_SR_3P1F_4e,lineWidth,1,kRed-7); h1D_m4l_SR_3P1F_4e->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_3P1F_4e, "Estimated from data","L");
        Shape = (TH1D*) h1D_m4l_SR_3P1F_4e->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);leg->AddEntry(Shape, "Smooth","L");
        // Shape->Draw("hist same goff");
        Shape2->Add(Shape2);
        // Shape2->Draw("hist same");//leg->AddEntry(Shape2, "Smooth 2*2P2F","L");
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_3P1F_4e.pdf");
        
        // Estimates from Control Regions
        // SR from 2P2F_4mu
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_2P2F_4mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 2P+2F 4mu");
        setHistProperties(h1D_m4l_SR_2P2F_4mu,lineWidth,1,kRed-7); h1D_m4l_SR_2P2F_4mu->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_2P2F_4mu, "Estimated from data","L");
        Shape2 = (TH1D*) h1D_m4l_SR_2P2F_4mu->Clone("Shape2");
        Shape2->Smooth();
        setHistProperties(Shape2,lineWidth,1,kGreen-7);//leg->AddEntry(Shape2, "Smooth ","L");
        // Shape2->Draw("hist same");
        
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_2P2F_4mu.pdf");
        
        
        
        
        
        
        // SR form 3P1F_4mu
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_3P1F_4mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 3P+1F 4mu");
        setHistProperties(h1D_m4l_SR_3P1F_4mu,lineWidth,1,kRed-7); h1D_m4l_SR_3P1F_4mu->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_3P1F_4mu, "Estimated from data","L");
        
        Shape = (TH1D*) h1D_m4l_SR_3P1F_4mu->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same goff");
        Shape2->Add(Shape2);//leg->AddEntry(Shape2, "Smooth 2*2P2F","L");
        //Shape2->Draw("hist same");
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_3P1F_4mu.pdf");
        
        
        
        // Estimates from Control Regions
        // SR from 2P2F 2e2mu
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_2P2F_2e2mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 2P+2F 2e2mu");
        setHistProperties(h1D_m4l_SR_2P2F_2e2mu,lineWidth,1,kRed-7); h1D_m4l_SR_2P2F_2e2mu->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_2P2F_2e2mu, "Estimated from data","L");
        Shape2 = (TH1D*) h1D_m4l_SR_2P2F_2e2mu->Clone("Shape2");
        Shape2->Smooth();
        setHistProperties(Shape2,lineWidth,1,kGreen-7);//leg->AddEntry(Shape2, "Smooth","L");
        //Shape2->Draw("hist same");
        
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_2P2F_2e2mu.pdf");
        
        
        
        // SR form 3P1F_2e2mu
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_3P1F_2e2mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 3P+1F 2e2mu");
        setHistProperties(h1D_m4l_SR_3P1F_2e2mu,lineWidth,1,kRed-7); h1D_m4l_SR_3P1F_2e2mu->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_3P1F_2e2mu, "Estimated from data","L");
        
        Shape = (TH1D*) h1D_m4l_SR_3P1F_2e2mu->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);//leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same");
        Shape2->Add(Shape2);
        // Shape2->Draw("hist same");leg->AddEntry(Shape2, "Smooth 2*2P2F","L");
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_3P1F_2e2mu.pdf");
        
        
        
        
        
        
        
        // Estimates from Control Regions
        // SR from 2P2F 2mu2e
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_2P2F_2mu2e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 2P+2F 2mu2e");
        setHistProperties(h1D_m4l_SR_2P2F_2mu2e,lineWidth,1,kRed-7); h1D_m4l_SR_2P2F_2mu2e->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_2P2F_2mu2e, "Estimated from data","L");
        Shape2 = (TH1D*) h1D_m4l_SR_2P2F_2mu2e->Clone("Shape2");
        Shape2->Smooth();
        setHistProperties(Shape2,lineWidth,1,kGreen-7);//leg->AddEntry(Shape2, "Smooth","L");
        //Shape2->Draw("hist same");
        
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_2P2F_2mu2e.pdf");
        
        
        
        // SR form 3P1F_2mu2e
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_3P1F_2mu2e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Z+X from 3P+1F 2mu2e");
        setHistProperties(h1D_m4l_SR_3P1F_2mu2e,lineWidth,1,kRed-7); h1D_m4l_SR_3P1F_2mu2e->Draw("same E1 goff"); leg->AddEntry(h1D_m4l_SR_3P1F_2mu2e, "Estimated from data","L");
        
        Shape = (TH1D*) h1D_m4l_SR_3P1F_2mu2e->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);//leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same");
        Shape2->Add(Shape2);
        // Shape2->Draw("hist same");leg->AddEntry(Shape2, "Smooth 2*2P2F","L");
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_3P1F_2mu2e.pdf");
        
        
        
        
        
        // SR total 4e
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_tot_4e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Total Z+X 4e");
        setHistProperties(h1D_m4l_SR_tot_4e,lineWidth,1,kRed-7); h1D_m4l_SR_tot_4e->Draw("same e1 goff"); leg->AddEntry(h1D_m4l_SR_tot_4e, "Estimated from data","L");
        
        TF1 *fsum = new TF1("fsum","landau( 0 )",0,2000);
        fsum->SetParameters(7.07240,132.73,18.0149);
        //double  par4e =0.00426549;
        //fsum->SetParameters(par4e,151.251,36.6082,7.06054,-0.00497032);
        fsum->SetLineColor(kGreen+2);
        fsum->SetLineStyle(1);
        fsum->SetLineWidth(2);
        //leg->AddEntry(fsum, "Fit","L");
        fsum->Draw("same");
        
        
        Shape = (TH1D*) h1D_m4l_SR_tot_4e->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same");
        
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_total_4e.pdf");
        
        // SR total 4mu
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_tot_4mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Total Z+X 4mu");
        setHistProperties(h1D_m4l_SR_tot_4mu,lineWidth,1,kRed-7); h1D_m4l_SR_tot_4mu->Draw("same e1 goff"); leg->AddEntry(h1D_m4l_SR_tot_4mu, "Estimated from data","L");
        
        Shape = (TH1D*) h1D_m4l_SR_tot_4mu->Clone("Shape");
        fsum = new TF1("fsum","landau( 0 )",0,2000);
        
        
        fsum->SetLineColor(kGreen+2);
        fsum->SetLineStyle(1);
        fsum->SetLineWidth(2);
        fsum->SetParameters(7.48631,121.229,18.087);
        
        // double par4mu = 1.84404*h1D_m4l_SR_tot_4mu->Integral(0,h1D_m4l_SR_tot_4mu->GetNbinsX(),"width")/fsum->Integral(50,1100);
        //fsum->SetParameters(par4mu,134.765,22.2293);
        //  cout<<"par4mu = "<<par4mu<<endl;
        // leg->AddEntry(fsum, "Fit","L");
        fsum->Draw("same");
        Shape->Smooth();
        
        
        setHistProperties(Shape,lineWidth,1,kBlue-7);//leg->AddEntry(Shape, "Smooth","L");
        //  Shape->Draw("hist same");
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_total_4mu.pdf");
        
        // SR total 2e2mu
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_tot_2e2mu->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Total Z+X 2e2mu");
        setHistProperties(h1D_m4l_SR_tot_2e2mu,lineWidth,1,kRed-7); h1D_m4l_SR_tot_2e2mu->Draw("same e1 goff"); leg->AddEntry(h1D_m4l_SR_tot_2e2mu, "Estimated from data","L");
        fsum = new TF1("fsum","landau(0)",0,2000);
        fsum->SetParameters(5.23196,128.313,18.5379);
        
        // double par2e2mu = 2.57065*h1D_m4l_SR_tot_2e2mu->Integral(0,h1D_m4l_SR_tot_2e2mu->GetNbinsX(),"width")/fsum->Integral(0,800);
        //fsum->SetParameters(par2e2mu,143.701,24.071);
        // cout<<"par2e2mu = "<<par2e2mu<<endl;
        
        fsum->SetLineColor(kGreen+2);
        fsum->SetLineStyle(1);
        fsum->SetLineWidth(2);
        
        //leg->AddEntry(fsum, "Fit","L");
        fsum->Draw("same");
        
        Shape = (TH1D*) h1D_m4l_SR_tot_2e2mu->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);//leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same");
        
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_total_2e2mu.pdf");
        
        
        
        
        
        
        // SR total 2mu2e
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_tot_2mu2e->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Total Z+X 2mu2e");
        setHistProperties(h1D_m4l_SR_tot_2mu2e,lineWidth,1,kRed-7); h1D_m4l_SR_tot_2mu2e->Draw("same e1 goff"); leg->AddEntry(h1D_m4l_SR_tot_2mu2e, "Estimated from data","L");
        fsum = new TF1("fsum","landau(0)",0,2000);
        fsum->SetParameters(9.73372,123.518,15.5824);
        
        //double par2mu2e = 2.57065*h1D_m4l_SR_tot_2mu2e->Integral(0,h1D_m4l_SR_tot_2mu2e->GetNbinsX(),"width")/fsum->Integral(0,800);
        // fsum->SetParameters(par2mu2e,143.701,24.071);
        //  cout<<"par2mu2e = "<<par2mu2e<<endl;
        
        fsum->SetLineColor(kGreen+2);
        fsum->SetLineStyle(1);
        fsum->SetLineWidth(2);
        
        //leg->AddEntry(fsum, "Fit","L");
        fsum->Draw("same");
        
        Shape = (TH1D*) h1D_m4l_SR_tot_2mu2e->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);//leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same");
        
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_total_2mu2e.pdf");
        
        
        
        
        // SR total
        h1D_dummyCR->GetXaxis()->SetRangeUser(0,800);
        h1D_dummyCR->SetMaximum(2.0*h1D_m4l_SR_tot->GetMaximum());
        h1D_dummyCR->Draw("goff"); cmsPreliminary(c1); leg = new TLegend(leg_xl,leg_yb,leg_xr,leg_yt); setLegendProperties(leg,"Total Z+X");
        setHistProperties(h1D_m4l_SR_tot,lineWidth,1,kRed-7); h1D_m4l_SR_tot->Draw("same e1 goff"); leg->AddEntry(h1D_m4l_SR_tot, "Estimated from data","L");
        
        fsum = new TF1("fsum","landau( 0 )+landau (3)+landau(6)+landau(9) ",50,1100);
        fsum->SetLineColor(kGreen+2);
        fsum->SetLineStyle(1);
        fsum->SetLineWidth(2);
        //fsum->SetParameters(par4e,151.251,36.6082,7.06054,-0.00497032,par4mu,134.765,22.2293,par2e2mu,143.701,24.071);
        //leg->AddEntry(fsum, "Fit","L");
        //fsum->Draw("same");
        
        Shape = (TH1D*) h1D_m4l_SR_tot->Clone("Shape");
        Shape->Smooth();
        setHistProperties(Shape,lineWidth,1,kBlue-7);leg->AddEntry(Shape, "Smooth","L");
        //Shape->Draw("hist same");
        
        leg->Draw("goff"); c1->SaveAs(sPlotsStore+"/SR_total.pdf");
    }
    
    
    /*     cout<<"----MASS4L REFIT---------------------------"<<endl;
     
     cout<<"h1D_m4l_SR_2P2F = "<<SR_mass4Refit[tot_2P2F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_4e = "<<SR_mass4Refit_4e[tot_2P2F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_4mu = "<<SR_mass4Refit_4mu[tot_2P2F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_2P2F_2e2mu = "<<SR_mass4Refit_2e2mu[tot_2P2F]->Integral()<<endl;
     cout<<"--------------------------------"<<endl;
     cout<<"h1D_m4l_SR_3P1F = "<<SR_mass4Refit[tot_3P1F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_4e = "<<SR_mass4Refit_4e[tot_3P1F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_4mu = "<<SR_mass4Refit_4mu[tot_3P1F]->Integral()<<endl;
     cout<<"h1D_m4l_SR_3P1F_2e2mu = "<<SR_mass4Refit_2e2mu[tot_3P1F]->Integral()<<endl;
     cout<<"--------------------------------"<<endl;
     
     cout<<"h1D_m4l_SR_tot = "<<SR_mass4Refit[tot]->Integral()<<endl;
     cout<<"h1D_m4l_SR_tot_4e = "<<SR_mass4Refit_4e[tot]->Integral()<<endl;
     cout<<"h1D_m4l_SR_tot_4mu = "<<SR_mass4Refit_4mu[tot]->Integral()<<endl;
     cout<<"h1D_m4l_SR_tot_2e2mu = "<<SR_mass4Refit_2e2mu[tot]->Integral()<<endl;
     cout<<"--------------------------------"<<endl;*/
    
    for (int i = 1 ; i<= D2plotmass4lDbkg_mass4l_nBin; i++)
        for (int j =1 ; j<=D2plotmass4lDbkg_Dbkg_nBin; j++)
        {
            
            D2plotmass4lDbkg->SetBinContent(i,j,D2plotmass4lDbkg->GetBinContent(i,j)+D2plotmass4lDbkg->Integral()/(1000*D2plotmass4lDbkg_mass4l_nBin*D2plotmass4lDbkg_Dbkg_nBin));
            
        }
    
    
    
    TH2D* D2plotmass4lDbkg_scaled1 =(TH2D *) D2plotmass4lDbkg->Clone("h2D_mass4l_Dbkg_scaledTo1");
    TH2D* D2plotmass4lDbkg_scaledm4lbins =new TH2D("h2D_mass4l_Dbkg_scaledm4lBins","2D_mass4l_Dbkg_scaledm4lBins",D2plotmass4lDbkg_mass4l_nBin, D2plotmass4lDbkg_mass4l_min, D2plotmass4lDbkg_mass4l_max,D2plotmass4lDbkg_Dbkg_nBin, D2plotmass4lDbkg_Dbkg_min, D2plotmass4lDbkg_Dbkg_max);
    
    
    
    
    D2plotmass4lDbkg_scaled1->Scale(1/D2plotmass4lDbkg_scaled1->Integral());
    
    for (int i = 1 ; i<= D2plotmass4lDbkg_mass4l_nBin; i++)
        for (int j =1 ; j<=D2plotmass4lDbkg_Dbkg_nBin; j++)
        {
            //cout<<i<<j<<endl;
            D2plotmass4lDbkg_scaledm4lbins->SetBinContent(i,j,D2plotmass4lDbkg->GetBinContent(i,j)*1/D2plotmass4lDbkg->Integral(i,i,1,D2plotmass4lDbkg_Dbkg_nBin) );
            
        }
    double in = 0;
    for (int i = 1 ; i<= D2plotmass4lDbkg_mass4l_nBin; i++)
        in+=D2plotmass4lDbkg_scaledm4lbins->Integral(i,i,1,D2plotmass4lDbkg_Dbkg_nBin);
    cout<<"Integr 2d  ="<<in<<endl;
    
    cout<<"Integr 2d  ="<<D2plotmass4lDbkg_scaledm4lbins->Integral()<<endl;
    
    D2plotmass4lDbkg_scaled1->Smooth();
    D2plotmass4lDbkg_scaledm4lbins->Smooth();
    
    
    D2plotmass4lDbkg_scaled1->Scale(1/D2plotmass4lDbkg_scaled1->Integral());
    
    
    for (int i = 1 ; i<= D2plotmass4lDbkg_mass4l_nBin; i++)
        for (int j =1 ; j<=D2plotmass4lDbkg_Dbkg_nBin; j++)
        {
            cout<<i<<j<<endl;
            D2plotmass4lDbkg_scaledm4lbins->SetBinContent(i,j,D2plotmass4lDbkg->GetBinContent(i,j)*1/D2plotmass4lDbkg->Integral(i,i,1,D2plotmass4lDbkg_Dbkg_nBin) );
            
        }
    in = 0;
    for (int i = 1 ; i<= D2plotmass4lDbkg_mass4l_nBin; i++)
        in+=D2plotmass4lDbkg_scaledm4lbins->Integral(i,i,1,D2plotmass4lDbkg_Dbkg_nBin);
    cout<<"Integr 2d  ="<<in<<endl;
    
    
    
    
    
    for (int i = mass4lREFIT; i<=mass4lErrREFIT;i++)
    {
        SR_shapes[i]->Smooth();
        SR_shapes[i] ->Scale(1/SR_shapes[i]->Integral());
        
        SR_shapes_4e[i]->Smooth();
        SR_shapes_4e[i]->Scale(1/SR_shapes_4e[i]->Integral());
        
        SR_shapes_4mu[i]->Smooth();
        SR_shapes_4mu[i]->Scale(1/SR_shapes_4mu[i]->Integral());
        
        SR_shapes_2e2mu[i]->Smooth();
        SR_shapes_2e2mu[i]->Scale(1/ SR_shapes_2e2mu[i]->Integral());
        
        
        
        CSR_shapes[i]->Smooth();
        CSR_shapes[i] ->Scale(1/CSR_shapes[i]->Integral());
        
        CSR_shapes_4e[i]->Smooth();
        CSR_shapes_4e[i]->Scale(1/CSR_shapes_4e[i]->Integral());
        
        CSR_shapes_4mu[i]->Smooth();
        CSR_shapes_4mu[i]->Scale(1/CSR_shapes_4mu[i]->Integral());
        
        CSR_shapes_2e2mu[i]->Smooth();
        CSR_shapes_2e2mu[i]->Scale(1/ CSR_shapes_2e2mu[i]->Integral());
        
        
        
    }
    
    
    cout<<"mass4Err_tot_scaled  = " << SR_shapes[mass4lREFIT]->Integral()<<endl;
    
    
    cout<<"SR 2P2F = "<<h1D_m4l_SR_2P2F->Integral();
    cout<<"SR 3P1F = "<<h1D_m4l_SR_3P1F->Integral();
    cout<<"TOtal Z+X = "<<h1D_m4l_SR_tot->Integral();
    
    
    
    
    // store slimmed tree and histograms in .root file
    TString fileNameTag = "estimatesZX"+fileName;
    TString fOption = "RECREATE";
    TFile* fTemplateTree = new TFile(sPlotsStore+"/"+fileNameTag+".root", fOption);
    fTemplateTree->cd();
    
    D2plotmass4lDbkg->Write();
    D2plotmass4lDbkg_scaled1->Write();
    D2plotmass4lDbkg_scaledm4lbins->Write();
    
    for (int i = tot; i<=tot_3P1F ; i++)
    {
        SR_mass4Refit[i] ->Write();
        SR_mass4Refit_4e[i] -> Write();
        SR_mass4Refit_4mu[i] ->Write();
        SR_mass4Refit_2e2mu[i] ->Write();
    }
    
    
    
    
    for (int i = mass4lREFIT; i<=mass4lErrREFIT;i++)
    {
        
        SR_shapes[i] ->Write();
        
        SR_shapes_4e[i]->Write();
        
        SR_shapes_4mu[i]->Write();
        
        SR_shapes_2e2mu[i]->Write();//greska
        
        CSR_shapes[i] ->Write();
        
        CSR_shapes_4e[i]->Write();
        
        CSR_shapes_4mu[i]->Write();
        
        CSR_shapes_2e2mu[i]->Write();//greska
        
        
    }
    
    
    
    
    
    h1D_FRel_EB->SetName("h1D_FRel_EB"); h1D_FRel_EB->Write();
    h1D_FRel_EE->SetName("h1D_FRel_EE"); h1D_FRel_EE->Write();
    h1D_FRmu_EB->SetName("h1D_FRmu_EB"); h1D_FRmu_EB->Write();
    h1D_FRmu_EE->SetName("h1D_FRmu_EE"); h1D_FRmu_EE->Write();
    h1D_m4l_SR_2P2F->SetName("h1D_m4l_SR_2P2F"); h1D_m4l_SR_2P2F->Write();
    h1D_m4l_SR_2P2F_4e->SetName("h1D_m4l_SR_2P2F_4e"); h1D_m4l_SR_2P2F_4e->Write();
    h1D_m4l_SR_2P2F_4mu->SetName("h1D_m4l_SR_2P2F_4mu"); h1D_m4l_SR_2P2F_4mu->Write();
    h1D_m4l_SR_2P2F_2e2mu->SetName("h1D_m4l_SR_2P2F_2e2mu"); h1D_m4l_SR_2P2F_2e2mu->Write();
    h1D_m4l_SR_2P2F_2mu2e->SetName("h1D_m4l_SR_2P2F_2mu2e"); h1D_m4l_SR_2P2F_2mu2e->Write();
    
    h1D_m4l_Add_2P2F->SetName("h1D_m4l_Add_2P2F"); h1D_m4l_Add_2P2F->Write();
    h1D_m4l_Add_2P2F_4e->SetName("h1D_m4l_Add_2P2F_4e"); h1D_m4l_Add_2P2F_4e->Write();
    h1D_m4l_Add_2P2F_4mu->SetName("h1D_m4l_Add_2P2F_4mu"); h1D_m4l_Add_2P2F_4mu->Write();
    h1D_m4l_Add_2P2F_2e2mu->SetName("h1D_m4l_Add_2P2F_2e2mu"); h1D_m4l_Add_2P2F_2e2mu->Write();
    h1D_m4l_Add_2P2F_2mu2e->SetName("h1D_m4l_Add_2P2F_2mu2e"); h1D_m4l_Add_2P2F_2mu2e->Write();
    
    
    h1D_m4l_SR_3P1F->SetName("h1D_m4l_SR_3P1F"); h1D_m4l_SR_3P1F->Write();
    h1D_m4l_SR_3P1F_4e->SetName("h1D_m4l_SR_3P1F_4e"); h1D_m4l_SR_3P1F_4e->Write();
    h1D_m4l_SR_3P1F_4mu->SetName("h1D_m4l_SR_3P1F_4mu"); h1D_m4l_SR_3P1F_4mu->Write();
    h1D_m4l_SR_3P1F_2e2mu->SetName("h1D_m4l_SR_3P1F_2e2mu"); h1D_m4l_SR_3P1F_2e2mu->Write();
    h1D_m4l_SR_3P1F_2mu2e->SetName("h1D_m4l_SR_3P1F_2mu2e"); h1D_m4l_SR_3P1F_2mu2e->Write();
    
    h1D_m4l_SR_tot->SetName("h1D_m4l_SR_tot"); h1D_m4l_SR_tot->Write();
    h1D_m4l_SR_tot_4e->SetName("h1D_m4l_SR_tot_4e"); h1D_m4l_SR_tot_4e->Write();
    h1D_m4l_SR_tot_4mu->SetName("h1D_m4l_SR_tot_4mu"); h1D_m4l_SR_tot_4mu->Write();
    h1D_m4l_SR_tot_2e2mu->SetName("h1D_m4l_SR_tot_2e2mu"); h1D_m4l_SR_tot_2e2mu->Write();
    h1D_m4l_SR_tot_2mu2e->SetName("h1D_m4l_SR_tot_2mu2e"); h1D_m4l_SR_tot_2mu2e->Write();
    
    fTemplateTree->Close();
    
    cout<<"h1D_m4l_SR_2P2F = "<<h1D_m4l_SR_2P2F->Integral()<<endl;
    cout<<"h1D_m4l_SR_2P2F_4e = "<<h1D_m4l_SR_2P2F_4e->Integral()<<endl;
    cout<<"h1D_m4l_SR_2P2F_4mu = "<<h1D_m4l_SR_2P2F_4mu->Integral()<<endl;
    cout<<"h1D_m4l_SR_2P2F_2e2mu = "<<h1D_m4l_SR_2P2F_2e2mu->Integral()<<endl;
    cout<<"--------------------------------"<<endl;
    /* cout<<"h1D_m4l_SR_2P2F_interval = "<<h1D_m4l_SR_2P2F->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_2P2F_interval = "<<h1D_m4l_SR_2P2F_4e->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_2P2F_interval = "<<h1D_m4l_SR_2P2F_4mu->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_2P2F_interval = "<<h1D_m4l_SR_2P2F_2e2mu->Integral(5,5)<<endl;*/
    cout<<"--------------------------------"<<endl;
    cout<<"h1D_m4l_SR_3P1F = "<<h1D_m4l_SR_3P1F->Integral()<<endl;
    cout<<"h1D_m4l_SR_3P1F_4e = "<<h1D_m4l_SR_3P1F_4e->Integral()<<endl;
    cout<<"h1D_m4l_SR_3P1F_4mu = "<<h1D_m4l_SR_3P1F_4mu->Integral()<<endl;
    cout<<"h1D_m4l_SR_3P1F_2e2mu = "<<h1D_m4l_SR_3P1F_2e2mu->Integral()<<endl;
    cout<<"--------------------------------"<<endl;
    /*  cout<<"h1D_m4l_SR_3P1F_interval = "<<h1D_m4l_SR_3P1F->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_3P1F_interval = "<<h1D_m4l_SR_3P1F_4e->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_3P1F_interval = "<<h1D_m4l_SR_3P1F_4mu->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_3P1F_interval = "<<h1D_m4l_SR_3P1F_2e2mu->Integral(5,5)<<endl;*/
    cout<<"--------------------------------"<<endl;
    cout<<"h1D_m4l_SR_tot = "<<h1D_m4l_SR_tot->Integral()<<endl;
    cout<<"h1D_m4l_SR_tot_4e = "<<h1D_m4l_SR_tot_4e->Integral()<<endl;
    cout<<"h1D_m4l_SR_tot_4mu = "<<h1D_m4l_SR_tot_4mu->Integral()<<endl;
    cout<<"h1D_m4l_SR_tot_2e2mu = "<<h1D_m4l_SR_tot_2e2mu->Integral()<<endl;
    cout<<"h1D_m4l_SR_tot_2mu2e = "<<h1D_m4l_SR_tot_2mu2e->Integral()<<endl;
    cout<<"--------------------------------"<<endl;
    /* cout<<"h1D_m4l_SR_tot_interval = "<<h1D_m4l_SR_tot->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_tot_4e_interval = "<<h1D_m4l_SR_tot_4e->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_tot_4mu_interval = "<<h1D_m4l_SR_tot_4mu->Integral(5,5)<<endl;
     cout<<"h1D_m4l_SR_tot_2e2mu_interval = "<<h1D_m4l_SR_tot_2e2mu->Integral(5,5)<<endl;*/
    
}

void PartOrigin(vector<int> *lep_matchedR03_PdgId,vector<int> *lep_matchedR03_MomId,vector<int> *lep_matchedR03_MomMomId,int lep_Hindex[],vector<int> *lep_id,  double FillValue,double FillWeight, TH1D *Hist_prompt, TH1D *Hist_fakes,TH1D *Hist_BDfakes, TH1D *Hist_conv, bool isZX)
{
    //if (lep_matchedR03_PdgId->at(lep_Hindex[2])==lep_id->at(lep_Hindex[2])) cout<<"Mom Id "<<lep_matchedR03_MomId->at(lep_Hindex[2])<<endl;
    if (!isData){
        n_pass++;
        //tau=15
        
        int n=0;
        
        
        
        n+=(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))==15)||(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==15)||(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==15);
        
        
        if (n>=1)
        {
            
            n_nonprompttaus++;
            if (!isZX)
                Hist_fakes->Fill(FillValue,FillWeight);
            
        }
        else
        {
            
            n=0;
            
            n+=lep_matchedR03_PdgId->at(lep_Hindex[2])==22 && (abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==11 || abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==13 || abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==15);
            
            if (n==1) {
                
                n_nonpromptconv++;
                if (!isZX)
                    Hist_conv->Fill(FillValue,FillWeight);
                
            }
            else
            {
                n=0;
                
                n+=lep_matchedR03_PdgId->at(lep_Hindex[2])==22 && !(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==11 || abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==13 || abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==15);
                /*
                 n+=(lep_matchedR03_PdgId->at(lep_Hindex[2])==22)&&((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==11)||(abs(lep_matchedR03_MomId->at(lep_Hindex[i]))==13));
                 n+=(lep_matchedR03_MomId->at(lep_Hindex[2])==22)&&((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))==11)||(abs(lep_matchedR03_PdgId->at(lep_Hindex[i]))==13));
                 n+=(lep_matchedR03_MomMomId->at(lep_Hindex[2])==22)&&((abs(lep_matchedR03_MomId->at(lep_Hindex[2))==11)||(abs(lep_matchedR03_MomId->at(lep_Hindex[i]))==13));
                 */
                if (n>=1)
                {
                    if (!isZX)
                        Hist_fakes->Fill(FillValue,FillWeight);
                    n_nonpromptformphot++;
                }
                else
                {
                    
                    n=0;
                    
                    n+=lep_matchedR03_PdgId->at(lep_Hindex[2])==23;
                    n+=((lep_matchedR03_MomId->at(lep_Hindex[2])==23)||(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))==24)||(lep_matchedR03_MomId->at(lep_Hindex[2])==25))&&(lep_matchedR03_PdgId->at(lep_Hindex[2])==lep_id->at(lep_Hindex[2]));
                    n+=((lep_matchedR03_MomMomId->at(lep_Hindex[2])==23)||(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))==24)||(lep_matchedR03_MomMomId->at(lep_Hindex[2])==25))&&(lep_matchedR03_PdgId->at(lep_Hindex[2])==lep_id->at(lep_Hindex[2]))&&(lep_matchedR03_MomId->at(lep_Hindex[2])==lep_id->at(lep_Hindex[2]));
                    if (n>=1)
                    {
                        if (!isZX)
                            Hist_prompt->Fill(FillValue,FillWeight);
                        n_prompt++;
                        
                    }
                    else
                    {
                        n=0;
                        n+=((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>410)&&(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<436))||((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>10410)&&(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<10433)) || ((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>20412)&&(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<20434));
                        
                        n+=((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>510)&&(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<546))||((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>10510)&&(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<10544)) || ((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>20512)&&(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<20544));
                        
                        
                        
                        n+=((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>410)&&(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))<436))||((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>10410)&&(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))<10433)) || ((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>20412)&&(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))<20434));
                        
                        n+=((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>510)&&(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))<546))||((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>10510)&&(abs(abs(lep_matchedR03_MomId->at(lep_Hindex[2])))<10544)) || ((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>20512)&&(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))<20544));
                        
                        
                        
                        
                        n+=((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>410)&&(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<436))||((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>10410)&&(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<10433)) || ((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>20412)&&(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<20434));
                        
                        n+=((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>510)&&(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<546))||((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>10510)&&(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<10544)) || ((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>20512)&&(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<20544));
                        
                        
                        
                        
                        
                        if (n>=1) {
                            
                            if (!isZX)
                                Hist_BDfakes->Fill(FillValue,FillWeight);
                            
                            n_bdfake++;
                            
                            
                        }
                        else
                        {
                            n=0;
                            
                            n+=((abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))<10)||(abs(lep_matchedR03_PdgId->at(lep_Hindex[2]))>100));
                            n+=((abs(lep_matchedR03_MomId->at(lep_Hindex[2]))<10)||(abs(lep_matchedR03_MomId->at(lep_Hindex[2]))>100));
                            n+=((abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))<10)||(abs(lep_matchedR03_MomMomId->at(lep_Hindex[2]))>100));
                            if (n>=1) {
                                n_fakes++;
                                
                                if (!isZX)
                                    Hist_fakes->Fill(FillValue,FillWeight);
                                
                            }
                            
                            else
                            {
                                n=0;
                                
                                if (!isZX)
                                    Hist_fakes->Fill(FillValue,FillWeight);
                                
                                n_others++;
                            }
                        }
                        
                    }
                }
                
            }
        }
    }
}


//_______________________________________________________________________________________________________________________________________________
void cmsPreliminary(TCanvas* &c){
    c->cd();
    
    TLatex *CMSPrelim = new TLatex();
    CMSPrelim->SetNDC(kTRUE);
    
    CMSPrelim->SetTextSize(0.5*c->GetTopMargin());
    CMSPrelim->SetTextFont(42);
    CMSPrelim->SetTextAlign(31); // align right
    CMSPrelim->DrawLatex(0.93, 0.96,"58.8 fb^{-1} at #sqrt{s} = 13 TeV");
    
    //    CMSPrelim->SetTextSize(0.9*c->GetTopMargin());
    //    CMSPrelim->SetTextFont(62);
    //    CMSPrelim->SetTextAlign(11); // align right
    //    CMSPrelim->DrawLatex(0.27, 0.85, "CMS");
    //
    //    CMSPrelim->SetTextSize(0.7*c->GetTopMargin());
    //    CMSPrelim->SetTextFont(52);
    //    CMSPrelim->SetTextAlign(11);
    //    CMSPrelim->DrawLatex(0.25, 0.8, "Preliminary");
}

//_______________________________________________________________________________________________________________________________________________
void analysisInit() {
    gErrorIgnoreLevel = kWarning;
    gErrorIgnoreLevel = kError;
}

//_______________________________________________________________________________________________________________________________________________
void setCavasAndStyles(TString canvasName, TCanvas* &c, TString stat, double leftMaring, double rightMaring, double bottomMaring, double topMaring){
    // setup environment
    //    gROOT->ProcessLine(".L setTDRStyle.C");
    setTDRStyle();
    // setup canvas
    c = new TCanvas(canvasName,"myPlots",0,0,800,600);
    c->cd(1); c->SetLogy(0);
    gStyle->SetOptStat(stat);
    gStyle->SetPalette(1);
    
    c->GetPad(0)->SetRightMargin(rightMaring);
    c->GetPad(0)->SetLeftMargin(leftMaring);
    c->GetPad(0)->SetTopMargin(topMaring);
    c->GetPad(0)->SetBottomMargin(bottomMaring);
}

//_______________________________________________________________________________________________________________________________________________
int setHistProperties(TH1D* &hist, Width_t lineWidth, Style_t lineStyle, Color_t lineColor, Style_t fillStyle, Color_t fillColor, TString xAxisTitle, TString yAxisTitle){
    if (!hist) return -1;
    // line
    hist->SetLineWidth(lineWidth);
    hist->SetLineStyle(lineStyle);
    hist->SetLineColor(lineColor);
    // fill
    hist->SetFillStyle(fillStyle);
    hist->SetFillColor(fillColor);
    // divisions, offsets, sizes
    hist->GetXaxis()->SetNdivisions(510);
    hist->GetYaxis()->SetNdivisions(510);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    // titles
    if (xAxisTitle!="skip") hist->GetXaxis()->SetTitle(xAxisTitle);
    if (yAxisTitle!="skip") hist->GetYaxis()->SetTitle(yAxisTitle);
    // return
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int setHistProperties2D(TH2D* &hist, TString xAxisTitle, TString yAxisTitle){
    if (!hist) return -1;
    // divisions, offsets, sizes
    hist->GetXaxis()->SetNdivisions(510);
    hist->GetYaxis()->SetNdivisions(510);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    // titles
    if (xAxisTitle!="skip") hist->GetXaxis()->SetTitle(xAxisTitle);
    if (yAxisTitle!="skip") hist->GetYaxis()->SetTitle(yAxisTitle);
    // return
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int setLegendProperties(TLegend* &leg, TString sHeader, Style_t fillStyle, Color_t fillColor){
    // sanity-check
    if (!leg) return -1;
    // titles
    if (sHeader!="skip") leg->SetHeader(sHeader);;
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    // return
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int normaliseHist2D(TH2D* &h2D, double norm){
    if (h2D->Integral()==0) return -1;
    h2D->Scale(norm/h2D->Integral());
    
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int fillEmptyBinsHist2D(TH2D* &h2D, double floor) {
    int nXbins=h2D->GetNbinsX();
    int nYbins=h2D->GetNbinsY();
    for(int i=1; i<=nXbins; i++){
        for(int j=1; j<=nYbins; j++){
            h2D->SetBinContent(i,j,h2D->GetBinContent(i,j)+floor);
        }
    }
    
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int normaliseHist1D(TH1D* &h1D, double norm){
    if (h1D->Integral()==0) return -1;
    h1D->Scale(norm/h1D->Integral());
    
    return 0;
}

//_______________________________________________________________________________________________________________________________________________
int fillEmptyBinsHist1D(TH1D* &h1D, double floor) {
    int nXbins=h1D->GetNbinsX();
    for(int i=1; i<=nXbins; i++){
        h1D->SetBinContent(i,h1D->GetBinContent(i)+floor);
    }
    
    return 0;
}

void Latex(TString capton,TString title [], double value [], double error[], int N)
{
    
    cout<<"\\begin{table}"<<endl;
    cout<<"\\centering"<<endl;
    cout<<"\\caption{"<<capton<<"}"<<endl;
    cout<<"\\label{Labela}"<<endl;
    cout<<"\\begin{tabular}{|c|c|}"<<endl;
    cout<<"\\hline"<<endl;
    cout<<"\\textbf{Type} & \\textbf{Percentage} \\\\ \\hline"<<endl;
    for (int i = 0; i < N; i++)
    {
        cout<<"\\textbf{"<<title[i]<<"} & \\textbf{$"<<value[i]<<"\\pm"<<error[i]<<"$} \\\\ \\hline"<<endl;
    }
    cout<<"\\end{tabular}"<<endl;
    cout<<"\\end{table}"<<endl;
    
    
}

void LatexULerr(TString capton,TString title [], double value [], double Uerror[],double Lerror[], int N)
{
    
    cout<<"\\begin{table}"<<endl;
    cout<<"\\centering"<<endl;
    cout<<"\\caption{"<<capton<<"}"<<endl;
    cout<<"\\label{Labela}"<<endl;
    cout<<"\\begin{tabular}{|c|c|}"<<endl;
    cout<<"\\hline"<<endl;
    cout<<"\\textbf{Type} & \\textbf{Percentage} \\\\ \\hline"<<endl;
    for (int i = 0; i < N; i++)
    {
        cout << fixed;
        cout.precision(5);
        cout<<"\\textbf{"<<title[i]<<"} & \\textbf{$"<<value[i]<<"^{+"<<Uerror[i]<<"}_{-"<<Lerror[i]<<"}$} \\\\ \\hline"<<endl;
    }
    cout<<"\\end{tabular}"<<endl;
    cout<<"\\end{table}"<<endl;
    
    
}

