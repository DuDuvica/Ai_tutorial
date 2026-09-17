#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include <iostream>
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLVUtils.h"
#include <cmath>
#include <vector>
#include <algorithm>
// Draft Z version based on AIW.C (Collins-Soper angles, A0-A7)
// To note that this Macro not only produce Ai form Phowegh or Sherpa MC but also create histogram for Kinematic studies.
// Ai are relevant only in the full-phacespace instead Kinematic plots are interesting also while applying Fiducial cut on the lepton. So a functionality FiducialCut is added to allow comparing basic Kinematic distibution wile applying on not fiducial cut on the leptons.
using namespace std;

bool sherpa = false;
bool test = true; // set to true for quick test with limited events; set to false for full run
bool override = false; // set to true to overwrite existing output file without prompt
bool normXS = true;
bool ifTrueOnly = true;
bool FiducialCut = true; // set to true to apply fiducial cuts at truth level, false to use all events (only relevant if ifTrueOnly=true)
bool FiducialCutEtaonly = false ;
bool FiducialCutCCCF = false ; // set to true to apply fiducial cuts at truth level, false to use all events (only relevant if ifTrueOnly=true)
bool FiducialCutCFonly = true ; // set to true to apply fiducial cuts at truth level, false to use all events (only relevant if ifTrueOnly=true)

// Weighted means with numerator/denominator covariance, including signed MC
// weights. Keep sums so bins with a cancelling denominator are identifiable.
namespace AIZP6 {
void Style(TH1* histogram) {
  histogram->SetStats(false);
  histogram->GetXaxis()->CenterTitle();
  histogram->GetYaxis()->CenterTitle();
  histogram->GetXaxis()->SetTitleSize(0.04);
  histogram->GetYaxis()->SetTitleSize(0.04);
  histogram->GetYaxis()->SetTitleOffset(1.5);
}
struct Moment {
  TH1D *mean, *sumW, *sumWX, *sumW2X, *sumW2X2, *valid;
  Moment(const TString& name, const TString& title, int n, double lo, double hi) {
    mean = new TH1D(name, title, n, lo, hi);
    sumW = new TH1D(name+"_sumW", title, n, lo, hi);
    sumWX = new TH1D(name+"_sumWX", title, n, lo, hi);
    sumW2X = new TH1D(name+"_sumW2X", title, n, lo, hi);
    sumW2X2 = new TH1D(name+"_sumW2X2", title, n, lo, hi);
    valid = new TH1D(name+"_valid", ";bin coordinate;1 = nonzero sum of weights", n, lo, hi);
    mean->Sumw2();
    sumW->Sumw2();
  }
  void Fill(double coordinate, double value, double weight) {
    sumW->Fill(coordinate, weight);
    sumWX->Fill(coordinate, weight*value);
    sumW2X->Fill(coordinate, weight*weight*value);
    sumW2X2->Fill(coordinate, weight*weight*value*value);
  }
  void Write() {
    for (int b=0; b<=mean->GetNbinsX()+1; ++b) {
      const double w = sumW->GetBinContent(b);
      if (w == 0.) continue; // undefined mean: valid mask stays zero
      const double mu = sumWX->GetBinContent(b)/w;
      const double w2 = std::pow(sumW->GetBinError(b), 2);
      const double variance = (sumW2X2->GetBinContent(b)
          - 2.*mu*sumW2X->GetBinContent(b) + mu*mu*w2)/(w*w);
      mean->SetBinContent(b, mu);
      mean->SetBinError(b, std::sqrt(std::max(0., variance)));
      valid->SetBinContent(b, 1.);
    }
    mean->Write(); sumW->Write(); sumWX->Write();
    sumW2X->Write(); sumW2X2->Write(); valid->Write();
  }
};
}

// Analytic basis plot: callable without any ntuple. Load TLVUtils.cxx first.
TCanvas* AIZPlotP6(const char* imageName="AIZ_P6_polynomial.pdf") {
  TH2D *map = new TH2D("P6_polynomial", ";cos#theta_{CS};#phi_{CS};P_{6}",
                       160, -1., 1., 160, 0., 2.*M_PI);
  map->SetDirectory(nullptr);
  std::vector<double> pols;
  for (int ix=1; ix<=map->GetNbinsX(); ++ix)
    for (int iy=1; iy<=map->GetNbinsY(); ++iy) {
      TLVUtils::getAiPolynoms(map->GetXaxis()->GetBinCenter(ix),
                             map->GetYaxis()->GetBinCenter(iy), pols);
      map->SetBinContent(ix, iy, pols[6]);
    }
  map->SetStats(false); map->SetMinimum(-1.); map->SetMaximum(1.);
  TCanvas *canvas = new TCanvas("c_P6_polynomial", "P6 angular basis", 1200, 500);
  canvas->Divide(2,1);
  canvas->cd(1)->SetRightMargin(0.17);
  map->Draw("COLZ");
  canvas->cd(2);
  TLegend *legend = new TLegend(0.57,0.72,0.88,0.88);
  for (int k=0; k<3; ++k) {
    const double phi = (k==0 ? M_PI/2. : (k==1 ? 3.*M_PI/2. : 0.));
    TH1D *slice = new TH1D(Form("P6_slice_%d", k),
        "P_{6} = sin(2#theta) sin#phi;cos#theta_{CS};P_{6}", 200, -1., 1.);
    slice->SetDirectory(nullptr);
    for (int b=1; b<=slice->GetNbinsX(); ++b) {
      TLVUtils::getAiPolynoms(slice->GetBinCenter(b), phi, pols);
      slice->SetBinContent(b, pols[6]);
    }
    slice->SetStats(false); slice->SetMinimum(-1.1); slice->SetMaximum(1.1);
    slice->SetLineColor(k==0 ? kRed+1 : (k==1 ? kBlue+1 : kBlack));
    slice->SetLineWidth(2);
    slice->Draw(k==0 ? "HIST" : "HIST SAME");
    legend->AddEntry(slice, k==0 ? "#phi = #pi/2" : (k==1 ? "#phi = 3#pi/2" : "#phi = 0"), "l");
  }
  legend->Draw();
  if (imageName && imageName[0]) canvas->SaveAs(imageName);
  return canvas;
}

// Macro to plot Ai coefficient from Sherpa and Powheg Z samples
void AIZ(bool isY=false){

  cout << " START AIZ " << endl;
  cout << "with configuration: " << endl;
  cout << " isY: " << isY << endl;
  cout << " ifTrueOnly: " << ifTrueOnly << endl;
  cout << " FiducialCut: " << FiducialCut << endl;
  cout << " FiducialCutEtaonly: " << FiducialCutEtaonly << endl;
  cout << " FiducialCutCCCF: " << FiducialCutCCCF << endl;
  cout << " FiducialCutCFonly: " << FiducialCutCFonly << endl;
  if ( (FiducialCutEtaonly || FiducialCutCCCF || FiducialCutCFonly) && !FiducialCut) {
    cout << " ERROR: FiducialCutEtaonly or FiducialCutCCCF or FiducialCutCFonly cannot be true if FiducialCut is false. Please set FiducialCut to true to apply eta-only fiducial cuts." << endl;
    return;
  }

  double xsecAMI = 0;

  const TString powhegFile = "/data/dust/group/atlas/zai/ntuples/CC_v24A/user.avallier.mc16_13TeV.361106.PP8EG_AZ_Zee.NTP_CC.e3601_s3126_r9364_p4239_vA24_outTree_CC_RecoZAi.root/user.avallier.36910521._000076.outTree_CC_RecoZAi.root";
  const TString sherpaFile = ""; // set to Sherpa Z sample if/when available
  const TString treePath = "HWWTree_ee";
  const TString TruthOnly = "/data/dust/group/atlas/zai/ntuples/CF_v23/unmerged/user.lbayer.mc16_13TeV.361106.PP8EG_AZ_Zee.NTP_CF.e3601_s3126_r10201_p4239_v23_outTree_ZAi.root/user.lbayer.38183008._000199.outTree_ZAi.root"; 

  if (sherpa) {
    cout << "Sherpa sample not define for ZAI tutorial " <<  endl;
    return;
  }

  TString minitree = ifTrueOnly ? TruthOnly : (sherpa ? sherpaFile : powhegFile);

  if (minitree.IsNull() || gSystem->AccessPathName(minitree)) {
    std::cout << " Input file not found for " << (sherpa ? "Sherpa" : "Powheg")
              << " configuration. Please set a valid path in AIZ.C" << std::endl;
    return;
  }

  TFile f1(minitree);
  if (f1.IsZombie()) {
    std::cout << " Unable to open input file " << minitree << std::endl;
    return;
  }

  // MetaData branches are floats in the input ntuples; read them as Float_t
  // to avoid SetBranchAddress type-mismatch errors, then promote to double for
  // calculations that follow.
  Float_t crossSectionPb_f = 0.0f;
  Float_t kFactor_f = 1.0f;
  Float_t filterEff_f = 1.0f;
  double crossSectionPb = 0.0;
  double kFactor = 1.0;
  double filterEff = 1.0;
  bool hasCrossSection = false;

  TTree *meta = static_cast<TTree*>(f1.Get("MetaData"));
  if (!meta) {
    std::cout << " MetaData tree missing, cannot read cross section" << std::endl;
    return;
  }

  if (meta->GetBranch("CrossSection")) {
    meta->SetBranchAddress("CrossSection", &crossSectionPb_f);
    hasCrossSection = true;
  } else if (meta->GetBranch("CrossSection_pb")) {
    meta->SetBranchAddress("CrossSection_pb", &crossSectionPb_f);
    hasCrossSection = true;
  } else if (meta->GetBranch("XSection")) {
    meta->SetBranchAddress("XSection", &crossSectionPb_f);
    hasCrossSection = true;
  }

  if (meta->GetBranch("Kfactor")) {
    meta->SetBranchAddress("Kfactor", &kFactor_f);
  }

  if (meta->GetBranch("FilterEfficiency")) {
    meta->SetBranchAddress("FilterEfficiency", &filterEff_f);
  } else if (meta->GetBranch("Br")) {
    meta->SetBranchAddress("Br", &filterEff_f);
  }

  meta->GetEntry(0);

  // Promote to double after reading to keep downstream code unchanged
  crossSectionPb = crossSectionPb_f;
  kFactor        = kFactor_f;
  filterEff      = filterEff_f;

  if (!hasCrossSection || crossSectionPb <= 0) {
    std::cout << " Missing or invalid cross section in MetaData" << std::endl;
    return;
  }

  if (kFactor <= 0) kFactor = 1.0;
  if (filterEff <= 0) filterEff = 1.0;

  xsecAMI = (crossSectionPb * kFactor * filterEff) / 1000.0; // convert pb to nb
  cout << " Cross section (pb): " << crossSectionPb << endl;
  cout << " K-factor: " << kFactor << endl;
  cout << " Filter efficiency: " << filterEff << endl;
  cout << " Final xsec (nb): " << xsecAMI << endl;
  TString outputName = "Zai_finalbinning";
  if (test) outputName = "test";
  if (sherpa) outputName = outputName+"Sherpa";
  else outputName = outputName+"Powheg";

  TString mode = "_pT";
  if (isY) mode = "_Y";
  if (normXS) mode = mode + "_NormXsec";

  std::cout << "Looking at File " << minitree << std::endl;
  TString prefix = ifTrueOnly ? "Truth_" : "";
  if ( FiducialCut ) prefix = prefix + "Fiducial_";
  if (FiducialCutEtaonly) prefix = prefix + "EtaOnly_";
  if (FiducialCutCCCF) prefix = prefix + "CCCF_";
  if (FiducialCutCFonly) prefix = prefix + "CFonly_";
  TString nameOutput = "AI_Z_"+prefix+outputName+mode+".root";
  TFile* Output = new TFile(nameOutput);
  bool isf = true;

  if (Output->IsZombie()) {
    std::cout << " OUTPUT FILE DO NOT EXIST" << std::endl;
    isf = false;
  }

  if (isf) {
    if (override) {
      Output  = new TFile(nameOutput,"RECREATE");
    } else {
      cout << "--> File " << nameOutput << " exist are you sure you want to override ? if yes put override to true " << endl;
      return;
    }
  } else {
    cout << " Create --> File " << nameOutput  << endl;
    Output  = new TFile(nameOutput,"RECREATE");
  }

  TH1D *hNorm = (TH1D*) f1.Get("CutFlow");
  if (!hNorm or hNorm->IsZombie()) {
    std::cout << " hNorm  DO NOT EXIST" << std::endl;
    return ;
  }

  // double norm=hNorm->GetBinContent(1)/hNorm->GetEntries();
  double denom = hNorm->GetBinContent(41);
  if (denom <= 0) {
    denom = hNorm->GetEntries();
    std::cout << " WARNING: normalization denom <= 0, using GetEntries() = " << denom << std::endl;
  }
  if (denom <= 0) {
    std::cout << " ERROR: normalization denom <= 0" << std::endl;
    return;
  }
  double norm = xsecAMI*1000/denom;

  cout << " Norm is: " << norm << endl;
  cout << " Xsec is:   " <<xsecAMI <<  " TOT ev GEn :  " <<  hNorm->GetBinContent(5) << " SUMM weight tot event   "  << hNorm->GetBinContent(1) << endl;

  // **********************
// Open the MAin TTree and set branch addresses
  // **********************
  cout << " \n Accessing tree " << treePath << " in file " << minitree << endl;
  TTree *tree = (TTree*)f1.Get(treePath);
  if (!tree) {
    std::cout << " Tree " << treePath << " DOES NOT EXIST" << std::endl;
    return;
  }

  const int Nbins = 24; //16; //11 
 //  Double_t bins[Nbins] = { 0., 8., 17., 27., 40., 55., 75., 110., 150., 210., 600.};
//Double_t bins[Nbins] = { 0., 2.5, 5., 8., 12., 15., 18., 25., 30., 40., 55., 75., 110., 150., 210., 600.};
  Double_t bins[]= {0,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0,32.6,36.4,40.4,44.9,50.2,56.4,63.9,73.4,85.4,105.0,132.0,173.0,253.0,600.0};

  if (isY) {
    double width = 0.4;
    for (int k =0 ; k<Nbins ;k++) {
      bins[k] = (width)*k;
      cout << " Y  bins at:  " << k << " bin " << bins[k] << endl;
    }
  }

  TH1D *hctheta = new TH1D("Costheta", "Costheta", 100, -1., 1.0);
  TH1D *hctheta_truth = new TH1D("CosthetaTruth", "CosthetaTruth", 100, -1., 1.0);
  TH1D *hphi = new TH1D("phi", "phi", 100, 0., 2.*M_PI);
  TH1D *hphi_truth = new TH1D("phiTruth", "phiTruth", 100, 0., 2.*M_PI);
  TH1D *Zmass = new TH1D("Zmass", ";m_{ll} [GeV];Events", 120, 60., 120.);
  TH1D *Xs = new TH1D("Xs", "Xs", Nbins-1, bins);
  TH1D *A0 = new TH1D("A0", "A0", Nbins-1, bins);
  TH1D *A1 = new TH1D("A1", "A1", Nbins-1, bins);
  TH1D *A2 = new TH1D("A2", "A2", Nbins-1, bins);
  TH1D *A3 = new TH1D("A3", "A3", Nbins-1, bins);
  TH1D *A4 = new TH1D("A4", "A4", Nbins-1, bins);
  TH1D *A5 = new TH1D("A5", "A5", Nbins-1, bins);
  TH1D *A6 = new TH1D("A6", "A6", Nbins-1, bins);
  TH1D *A7 = new TH1D("A7", "A7", Nbins-1, bins);
  TH1D *Xsw = new TH1D("Xsweighted", "Xsweighted", Nbins-1, bins);
  TH1D *ZMass = new TH1D("ZMass", "ZMass", 100, 50., 150.);

  // Additional truth-level angular distributions vs lepton pT
  const int pt2dBins = 40;
  const double pt2dMax = 200.; // GeV
  const int etaBins = 60;
  const double etaMin = -5.0;
  const double etaMax = 5.0;
  const int phiBins = 64;
  const double phiMin = -M_PI;
  const double phiMax = M_PI;
  const double phiCSMin = 0.;
  const double phiCSMax = 2.*M_PI;
  const int zPtBins = 60;
  const double zPtMax = 300.;
  const int zYBins = 50;
  const double zYMin = -5.0;
  const double zYMax = 5.0;
  TH2D *hCosVspt_el = new TH2D("cosThetaCSTruth_vs_lepPtNeg",
                             ";p_{T}^{truth}(e^{-}) [GeV];cos#theta_{CS}^{truth}",
                             pt2dBins, 0., pt2dMax, 50, -1., 1.);
  TH2D *hCosVspt_pos = new TH2D("cosThetaCSTruth_vs_lepPtPos",
                             ";p_{T}^{truth}(e^{+}) [GeV];cos#theta_{CS}^{truth}",
                             pt2dBins, 0., pt2dMax, 50, -1., 1.);
  TH2D *hPhiVspt_el = new TH2D("phiCSTruth_vs_lepPtNeg",
                             ";p_{T}^{truth}(e^{-}) [GeV];#phi_{CS}^{truth}",
                             pt2dBins, 0., pt2dMax, 50, phiCSMin, phiCSMax);
  TH2D *hPhiVspt_pos = new TH2D("phiCSTruth_vs_lepPtPos",
                             ";p_{T}^{truth}(e^{+}) [GeV];#phi_{CS}^{truth}",
                             pt2dBins, 0., pt2dMax, 50, phiCSMin, phiCSMax);

  // Lepton-lepton eta/phi correlations
  TH2D *hEta_ep_vs_em = new TH2D("eta_ep_vs_em", ";#eta(e^{+});#eta(e^{-})", etaBins, etaMin, etaMax, etaBins, etaMin, etaMax);
  TH2D *hEta_eleading_vs_esubleading = new TH2D("eta_eleading_vs_esubleading", ";|#eta(e_{leading})|;|#eta(e_{subleading})|", etaBins/2, 0, etaMax, etaBins/2, 0, etaMax);
  TH2D *hPhi_ep_vs_em = new TH2D("phi_ep_vs_em", ";#phi(e^{+});#phi(e^{-})", phiBins, phiMin, phiMax, phiBins, phiMin, phiMax);
  TH2D *hPt_ep_vs_em = new TH2D("pt_ep_vs_em", ";p_{T}(e^{+}) [GeV];p_{T}(e^{-}) [GeV]", pt2dBins, 0., pt2dMax, pt2dBins, 0., pt2dMax);
  TH2D *hPt_leading_vs_subleading = new TH2D("pt_leading_vs_subleading", ";p_{T}(leading l) [GeV];p_{T}(subleading l) [GeV]", pt2dBins, 0., pt2dMax, pt2dBins, 0., pt2dMax);

  // Z–lepton correlations
  TH2D *hZptVsEta_m = new TH2D("zPt_vs_eta_m", ";p_{T}(Z) [GeV];#eta(e^{-})", zPtBins, 0., zPtMax, etaBins, etaMin, etaMax);
  TH2D *hZptVsEta_p = new TH2D("zPt_vs_eta_p", ";p_{T}(Z) [GeV];#eta(e^{+})", zPtBins, 0., zPtMax, etaBins, etaMin, etaMax);
  TH2D *hZptVsEta_forward = new TH2D("zPt_vs_eta_forward", ";p_{T}(Z) [GeV];|#eta(forward l)|", zPtBins, 0., zPtMax, etaBins/2, 0, etaMax);
  TH2D *hZptVsEta_central = new TH2D("zPt_vs_eta_central", ";p_{T}(Z) [GeV];|#eta(central l)|", zPtBins, 0., zPtMax, etaBins/2, 0, etaMax);
  TH2D *hZptVsPt_forward = new TH2D("zPt_vs_pt_forward", ";p_{T}(Z) [GeV];p_{T}(forward l) [GeV]", zPtBins, 0., zPtMax, pt2dBins, 0., pt2dMax);
  TH2D *hZptVsPt_central = new TH2D("zPt_vs_pt_central", ";p_{T}(Z) [GeV];p_{T}(central l) [GeV]", zPtBins, 0., zPtMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsPt_m   = new TH2D("zY_vs_pt_m",   ";|y(Z)|;p_{T}(e^{-}) [GeV]", zYBins/2, 0., zYMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsPt_p   = new TH2D("zY_vs_pt_p",   ";|y(Z)|;p_{T}(e^{+}) [GeV]", zYBins/2, 0., zYMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsPt_leading = new TH2D("zY_vs_pt_leading", ";|y(Z)|;p_{T}(leading l) [GeV]", zYBins/2, 0., zYMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsPt_subleading = new TH2D("zY_vs_pt_subleading", ";|y(Z)|;p_{T}(subleading l) [GeV]", zYBins/2, 0., zYMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsPt_forward = new TH2D("zY_vs_pt_forward", ";|y(Z)|;p_{T}(forward l) [GeV]", zYBins/2, 0., zYMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsPt_central = new TH2D("zY_vs_pt_central", ";|y(Z)|;p_{T}(central l) [GeV]", zYBins/2, 0., zYMax, pt2dBins, 0., pt2dMax);
  TH2D *hZYVsEta_leading = new TH2D("zY_vs_eta_leading", ";|y(Z)|;#eta(leading l)", zYBins/2, 0., zYMax, etaBins, etaMin, etaMax);
  TH2D *hZYVsEta_subleading = new TH2D("zY_vs_eta_subleading", ";|y(Z)|;#eta(subleading l)", zYBins/2, 0., zYMax, etaBins, etaMin, etaMax);
  TH2D *hZYVsEta_forward = new TH2D("zY_vs_eta_forward", ";|y(Z)|;|#eta(forward l)|", zYBins/2, 0., zYMax, etaBins/2, 0, etaMax);
  TH2D *hZYVsEta_central = new TH2D("zY_vs_eta_central", ";|y(Z)|;|#eta(central l)|", zYBins/2, 0., zYMax, etaBins/2, 0, etaMax);
  TH2D *hZptVsCostheta = new TH2D("zPt_vs_costheta", ";p_{T}(Z) [GeV];cos#theta_{CS}", zPtBins, 0., zPtMax, 50, -1., 1.);
  TH2D *hZptVsPhi = new TH2D("zPt_vs_phi", ";p_{T}(Z) [GeV];#phi_{CS}", zPtBins, 0., zPtMax, phiBins, phiCSMin, phiCSMax);
  TH2D *hZYVsCostheta = new TH2D("zY_vs_costheta", ";|y(Z)|;cos#theta_{CS}", zYBins/2, 0., zYMax, 50, -1., 1.);
  TH2D *hZYVsPhi = new TH2D("zY_vs_phi", ";|y(Z)|;#phi_{CS}", zYBins/2, 0., zYMax, phiBins, phiCSMin, phiCSMax);

  // Leading and subleading lepton distributions
  TH2D *hEtaVsPt_leading = new TH2D("eta_vs_pt_leading", ";p_{T}(lead) [GeV];#eta(lead)", pt2dBins, 0., pt2dMax, etaBins, etaMin, etaMax);
  TH2D *hEtaVsPt_subleading = new TH2D("eta_vs_pt_subleading", ";p_{T}(sublead) [GeV];#eta(sublead)", pt2dBins, 0., pt2dMax, etaBins, etaMin, etaMax);
  TH2D *hCosthVsPt_leading = new TH2D("costh_vs_pt_leading", ";p_{T}(lead) [GeV];cos#theta_{CS}", pt2dBins, 0., pt2dMax, 50, -1., 1.);
  TH2D *hCosthVsPt_subleading = new TH2D("costh_vs_pt_subleading", ";p_{T}(sublead) [GeV];cos#theta_{CS}", pt2dBins, 0., pt2dMax, 50, -1., 1.);

  // Lepton-lepton angular separation:
  // close-by pairs -> small DeltaR and |DeltaPhi|; back-to-back -> |DeltaPhi| ~ pi and cos(opening) ~ -1.
  // DeltaEta can be large for forward leptons, but should be symmetric around 0 for Z->ll.
  // DeltaR(ll) distribution shape is influenced by lepton pT cuts and Z pT spectrum; low DeltaR region can be depleted by isolation requirements in reconstructed-level analyses.
  // DeltaEta show separation alcong the beam axis and can be sensitive to PDFs and higher-order effects; should be symmetric around 0 for Z->ll.
  // Cosine of opening angle between leptons in lab frame can provide complementary information to DeltaR and DeltaPhi, especially for events where leptons are close in angle but not necessarily back-to-back.
  TH1D *hDeltaR_ll = new TH1D("deltaR_ll", ";#DeltaR(l_{1},l_{2});Events", 80, 0., 8.0);
  TH1D *hDeltaPhi_ll = new TH1D("deltaPhi_ll", ";|#Delta#phi(l_{1},l_{2})|;Events", 64, 0., M_PI);
  TH1D *hDeltaEta_ll = new TH1D("deltaEta_ll", ";|#Delta#eta(l_{1},l_{2})|;Events", 100, 0.0, 10.0);
  TH2D *hDeltaEtaVsDeltaPhi_ll = new TH2D("deltaEta_vs_deltaPhi_ll", ";|#Delta#eta(l_{1},l_{2})|;|#Delta#phi(l_{1},l_{2})|", 100, 0.0, 10.0, 64, 0., M_PI);
  TH2D *hDeltaEtaVsCosth_ll = new TH2D("deltaEta_vs_costh_ll", ";|#Delta#eta(l_{1},l_{2})|;cos#theta_{CS}", 100, 0.0, 10.0, 50, -1., 1.);
  TH2D *hDeltaEtaVsZPt_ll = new TH2D("deltaEta_vs_zPt_ll", ";|#Delta#eta(l_{1},l_{2})|;p_{T}(Z) [GeV]", 100, 0.0, 10.0, zPtBins, 0., zPtMax);
  TH2D *hDeltaEtaVsZY_ll = new TH2D("deltaEta_vs_zY_ll", ";|#Delta#eta(l_{1},l_{2})|;|y(Z)|", 100, 0.0, 10.0, zYBins/2, 0., zYMax);
  TH2D *hDeltaEtaVsLeadPt_ll = new TH2D("deltaEta_vs_leadingPt_ll", ";|#Delta#eta(l_{1},l_{2})|;p_{T}(leading l) [GeV]", 100, 0.0, 10.0, pt2dBins, 0., pt2dMax);
  TH2D *hDeltaEtaVsSubleadPt_ll = new TH2D("deltaEta_vs_subleadingPt_ll", ";|#Delta#eta(l_{1},l_{2})|;p_{T}(subleading l) [GeV]", 100, 0.0, 10.0, pt2dBins, 0., pt2dMax);
 
  // cos opening anle between the two leptons in the lab frame, defined as the cosine of the angle between their three-momenta: cos(opening) = (p1 . p2) / (|p1| |p2|).
  // back-to-back leptons will have cos(opening) close to -1, while collinear leptons will have cos(opening) close to +1. 
  // The distribution of cos(opening) can provide insights into the kinematics of the Z decay and the effects of QCD radiation, especially in boosted regimes where leptons may be close in angle but not necessarily back-to-back.
  // Note: for massless leptons, cos(opening) = tanh(DeltaEta/2)*cos(DeltaPhi) in the limit where the leptons are back-to-back in phi; 
  // however, for boosted Z bosons or when leptons are close in angle, this relationship does not hold and cos(opening) can provide additional insight into the lepton kinematics and event topology. 
  TH1D *hCosOpening_ll = new TH1D("cosOpeningAngle_ll", ";cos(#alpha_{l_{1}l_{2}});Events", 100, -1.0, 1.0);
  TH2D *hCosOpeningVsZPt_ll = new TH2D("cosOpeningAngle_vs_zPt_ll", ";p_{T}(Z) [GeV];cos(#alpha_{l_{1}l_{2}})", zPtBins, 0., zPtMax, 100, -1.0, 1.0);
  TH2D *hCosOpeningVsLeadPt_ll = new TH2D("cosOpeningAngle_vs_leadingPt_ll", ";p_{T}(leading l) [GeV];cos(#alpha_{l_{1}l_{2}})", pt2dBins, 0., pt2dMax, 100, -1.0, 1.0);
  TH2D *hCosOpeningVsSubleadPt_ll = new TH2D("cosOpeningAngle_vs_subleadingPt_ll", ";p_{T}(subleading l) [GeV];cos(#alpha_{l_{1}l_{2}})", pt2dBins, 0., pt2dMax, 100, -1.0, 1.0);
  TH2D *hCosOpeningVsDeltaEta_ll = new TH2D("cosOpeningAngle_vs_deltaEta_ll", ";|#Delta#eta(l_{1},l_{2})|;cos(#alpha_{l_{1}l_{2}})", 100, 0.0, 10.0, 100, -1.0, 1.0);

  // cosTheta slices in pT bins (GeV) for each lepton
  TH1D *hCosLneg0_5   = new TH1D("cosThetaCSTruth_neg_pt_elto5",   "cosThetaCSTruth_neg_pt_elto5", 50, -1., 1.);
  TH1D *hCosLneg5_20  = new TH1D("cosThetaCSTruth_neg_pt5to20",  "cosThetaCSTruth_neg_pt5to20", 50, -1., 1.);
  TH1D *hCosLneg20_40 = new TH1D("cosThetaCSTruth_neg_pt20to40", "cosThetaCSTruth_neg_pt20to40", 50, -1., 1.);
  TH1D *hCosLneg40_80 = new TH1D("cosThetaCSTruth_neg_pt40to80", "cosThetaCSTruth_neg_pt40to80", 50, -1., 1.);
  TH1D *hCosLpos0_5   = new TH1D("cosThetaCSTruth_pos_pt_elto5",   "cosThetaCSTruth_pos_pt_elto5", 50, -1., 1.);
  TH1D *hCosLpos5_20  = new TH1D("cosThetaCSTruth_pos_pt5to20",  "cosThetaCSTruth_pos_pt5to20", 50, -1., 1.);
  TH1D *hCosLpos20_40 = new TH1D("cosThetaCSTruth_pos_pt20to40", "cosThetaCSTruth_pos_pt20to40", 50, -1., 1.);
  TH1D *hCosLpos40_80 = new TH1D("cosThetaCSTruth_pos_pt40to80", "cosThetaCSTruth_pos_pt40to80", 50, -1., 1.);

  // P6-sensitive diagnostics use the same selected events and CS convention.
  // A6_moment = 5<P6> is an Ai estimator ONLY with full angular acceptance.
  TH1D *hP6 = new TH1D("P6", ";P_{6};Weighted events", 100, -1.000001, 1.000001);
  TH1D *hAbsP6 = new TH1D("absP6", ";|P_{6}|;Weighted events", 50, 0., 1.000001);
  TH2D *hP6ForwardCharge = new TH2D("P6_vs_forwardCharge_CF",
      ";Charge of forward lepton (CF only);P_{6}", 2, -2., 2., 100, -1.000001, 1.000001);
  hP6->Sumw2(); hAbsP6->Sumw2(); hP6ForwardCharge->Sumw2();
  const char* p6Axes[] = {"deltaEta", "zPt", "absZY", "acoplanarity"};
  const char* p6Labels[] = {"|#Delta#eta_{ee}|", "p_{T}(Z) [GeV]", "|y(Z)|", "#pi-|#Delta#phi_{ee}|"};
  const int p6NBins[] = {40, 60, 25, 32};
  const double p6Max[] = {10., 300., 5., M_PI};
  std::vector<TH2D*> p6Correlations;
  std::vector<AIZP6::Moment*> p6Moments, p6Signs, p6LeverArms;
  for (int j=0; j<4; ++j) {
    p6Correlations.push_back(new TH2D(Form("P6_vs_%s", p6Axes[j]),
        Form(";%s;P_{6}", p6Labels[j]), p6NBins[j], 0., p6Max[j], 100, -1.000001, 1.000001));
    p6Correlations.back()->Sumw2();
    p6Moments.push_back(new AIZP6::Moment(Form("A6_moment_vs_%s", p6Axes[j]),
        Form("Selected-sample moment;%s;5#LT P_{6}#GT", p6Labels[j]), p6NBins[j], 0., p6Max[j]));
    p6Signs.push_back(new AIZP6::Moment(Form("P6_signAsym_vs_%s", p6Axes[j]),
        Form("Selected-sample sign asymmetry;%s;#LT sign(P_{6})#GT", p6Labels[j]), p6NBins[j], 0., p6Max[j]));
    p6LeverArms.push_back(new AIZP6::Moment(Form("P6_squared_vs_%s", p6Axes[j]),
        Form("Angular lever arm;%s;#LT P_{6}^{2}#GT", p6Labels[j]), p6NBins[j], 0., p6Max[j]));
  }
  Long64_t p6UndefinedPlane = 0;

  Long64_t N = tree->GetEntries();
  if (test) N = std::min<Long64_t>(N, 100000);

  // Tree stores these as Float_t; read with matching types to avoid ROOT warnings
  Float_t lepPtTruth0=-99, lepPhiTruth0 =-99, lepEtaTruth0=-99;
  Float_t lepPtTruth1=-99, lepPhiTruth1 =-99, lepEtaTruth1=-99;
  double_t mcEventWeight=1.0;
  Float_t cosThetaCSTruth=0, phiCSTruth=0 ;
    Float_t         lepID0=1;
   Float_t         lepID1=-1;

  tree->SetBranchAddress("mcEventWeight", &mcEventWeight);
  tree->SetBranchAddress("lepPtTruth0", &lepPtTruth0);
  tree->SetBranchAddress("lepPhiTruth0", &lepPhiTruth0);
  tree->SetBranchAddress("lepEtaTruth0", &lepEtaTruth0);
  tree->SetBranchAddress("lepPtTruth1", &lepPtTruth1);
  tree->SetBranchAddress("lepPhiTruth1", &lepPhiTruth1);
  tree->SetBranchAddress("lepEtaTruth1", &lepEtaTruth1);
  if (!ifTrueOnly) {
      tree->SetBranchAddress("cosThetaCSTruth", &cosThetaCSTruth);
      tree->SetBranchAddress("phiCSTruth", &phiCSTruth);
      tree->SetBranchAddress("lepID0", &lepID0);
      tree->SetBranchAddress("lepID1", &lepID1);
  }
  TLorentzVector em, ep, z;

  // Beam energy per proton at 13 TeV collision energy.
  // This value must be in the same units as the lepton four-vectors (GeV).
  const double ebeamGeV = 6500.0;

  for (Long64_t i = 0; i < N; i++) {

    tree->GetEntry(i);
    // Initialize four-vectors to avoid carrying over values from previous events
    // in case of early continue statements.
      em.SetPtEtaPhiM(0, 0, 0, 0);
      ep.SetPtEtaPhiM(0, 0, 0, 0);
      z.SetPtEtaPhiM(0, 0, 0, 0);
    const int id0 = lround(lepID0);
    const int id1 = lround(lepID1);

    if (id0 * id1 > 0) {
        cout << " WARNING: Event " << i << " has leptons with same charge (lepID0 = " << lepID0 << ", lepID1 = " << lepID1 << "). Skipping event." << endl;
        // Skip events where both leptons have the same charge
        continue;
    }

    
    // We build em = e- and ep = e+ explicitly for each event so that downstream
    // angular definitions are unambiguous. In particular, TLVUtils::getCSFAngles
    // expects the first lepton argument to be associated with charge1.
    //  e−  == 11 e+ == -11 
    if ( id0 == 11 && id1 == -11){ 
         ep.SetPtEtaPhiM(lepPtTruth1/1000.0, lepEtaTruth1, lepPhiTruth1, 0);
         em.SetPtEtaPhiM(lepPtTruth0/1000.0, lepEtaTruth0, lepPhiTruth0, 0);
    }else if (id0 == -11 && id1 == 11){
      ep.SetPtEtaPhiM(lepPtTruth0/1000.0, lepEtaTruth0, lepPhiTruth0, 0);
      em.SetPtEtaPhiM(lepPtTruth1/1000.0, lepEtaTruth1, lepPhiTruth1, 0);
    }else {
        if (ifTrueOnly) {
          // IF truthOnly Ntuple already selected and the first lep is positive
          //   m_miniOutTree->m_lepPtTruth0 = float(pos.Pt());
          //   m_miniOutTree->m_lepPtTruth1 = float(neg.Pt());
          em.SetPtEtaPhiM(lepPtTruth1/1000.0, lepEtaTruth1, lepPhiTruth1, 0);
          ep.SetPtEtaPhiM(lepPtTruth0/1000.0, lepEtaTruth0, lepPhiTruth0, 0);
        } else {
          // Skip events where leptons are not identified as e+e-
          cout << " WARNING: Event " << i << " has unexpected lepton IDs: lepID0 = " << lepID0 << ", lepID1 = " << lepID1 << ". Skipping event." << endl;
          continue; 
        }
    }

    if ( ifTrueOnly && FiducialCut ) {
      if (FiducialCutEtaonly) {
        // Apply only eta cuts for fiducial selection in CC
        if (fabs(em.Eta()) > 2.5 || fabs(ep.Eta()) > 2.5) {
          // Skip events where leptons do not pass eta cuts
          if ( i%10000 == 0 )  cout << " WARNING: Event " << i << " fails eta-only fiducial cuts: eta_el = " << em.Eta()
                << "; eta_pos = " << ep.Eta() << ". Skipping event." << endl;
          continue;
        }
      } else if (FiducialCutCCCF)   {
      // CF fiducial selection:
      // both leptons must satisfy pT > 25 GeV and at least one lepton must be central (|eta| < 2.5).
      // Boundary |eta| = 2.5 is treated as CF (forward).
        if (em.Pt() < 25.0 || ep.Pt() < 25.0 || (fabs(em.Eta()) >= 2.5 && fabs(ep.Eta()) >= 2.5)) {
          if ( i%10000 == 0 )  cout << " WARNING: Event " << i << " fails CF fiducial cuts: pt_el = " << em.Pt()
                << " GeV, eta_el = " << em.Eta() << "; pt_pos = " << ep.Pt()
                << " GeV, eta_pos = " << ep.Eta() << ". Skipping event." << endl;
          continue;
        }
      } else if (FiducialCutCFonly) {
        // CF-only fiducial selection:
        // both leptons must satisfy pT > 25 GeV and the pair must be one central
        // lepton (|eta| < 2.5) plus one forward lepton (|eta| >= 2.5).
        if (em.Pt() < 25.0 || ep.Pt() < 25.0 ||
            !((fabs(em.Eta()) < 2.5 && fabs(ep.Eta()) >= 2.5) ||
              (fabs(em.Eta()) >= 2.5 && fabs(ep.Eta()) < 2.5))) {
          if ( i%10000 == 0 )  cout << " WARNING: Event " << i << " fails CF-only fiducial cuts: pt_el = " << em.Pt()
                << " GeV, eta_el = " << em.Eta() << "; pt_pos = " << ep.Pt()
                << " GeV, eta_pos = " << ep.Eta() << ". Skipping event." << endl;
          continue;
        }
      }
       else {
         // Apply full fiducial cuts (CC) on both pT and eta
        if (em.Pt() < 25.0 || fabs(em.Eta()) > 2.5 || ep.Pt() < 25.0 || fabs(ep.Eta()) > 2.5) {
          // Skip events where leptons do not pass fiducial cuts
         if ( i%10000 == 0 )  cout << " WARNING: Event " << i << " fails fiducial cuts: pt_el = " << em.Pt() << " GeV, eta_el = " << em.Eta()
               << "; pt_pos = " << ep.Pt() << " GeV, eta_pos = " << ep.Eta() << ". Skipping event." << endl;
          continue;
        }
      }
    }

    z = em + ep;

    if (z.M() < 66. || z.M() > 116.) {
       // cout << " WARNING: Event " << i << " has dilepton mass outside Z window: m_ll = " << z.M() << " GeV. Skipping event." << endl;
        // Skip events outside the Z mass window
    continue;
    }
    ZMass->Fill(z.M(), mcEventWeight);

    if (i%10000 == 0) {
      cout << "Processed " << i <<"/" << N << " events" << endl;
      cout << " Z pT is: " << z.Pt() << " Z rapidity is: " << z.Rapidity() << endl;
      cout << " lepton 1 (neg) pT: " << em.Pt() << " eta: " << em.Eta() << " phi: " << em.Phi() << endl;
      cout << " Event weight: " << mcEventWeight << endl;
    }

    // Use the exact Collins-Soper implementation from TLVUtils so AIZ and
    // helper-library conventions remain identical:
    //  - same sign conventions for cos(theta)
    //  - same charge handling (first lepton carries charge1)
    //  - same boosted-frame construction for phi
    // We pass em (negative lepton) as lep1 with charge1 = -1 and ep as lep2.
    double costheta = 0.0;
    double phi = 0.0;
    TLVUtils::getCSFAngles(em, -1, ep, ebeamGeV, costheta, phi); 

    // Guard against non-finite outputs before filling histograms.
    if (!std::isfinite(costheta) || !std::isfinite(phi)) continue;

    double weight = 1.0;
    if (normXS) weight = mcEventWeight; //weight *= norm;
    const double phiCSTruthWrapped = TVector2::Phi_0_2pi(phiCSTruth);

    Zmass->Fill(z.M(), weight);
    hctheta->Fill(costheta, weight);
    hctheta_truth->Fill(cosThetaCSTruth, weight);
    hphi_truth->Fill(phiCSTruthWrapped, weight);

    hphi->Fill(phi, weight);
    hZptVsCostheta->Fill(z.Pt(), costheta, weight);
    hZptVsPhi->Fill(z.Pt(), phi, weight);
    hZYVsCostheta->Fill(fabs(z.Rapidity()), costheta, weight);
    hZYVsPhi->Fill(fabs(z.Rapidity()), phi, weight);

    // Determine leading and subleading leptons
    double pt_leading = (em.Pt() > ep.Pt()) ? em.Pt() : ep.Pt();
    double pt_subleading = (em.Pt() > ep.Pt()) ? ep.Pt() : em.Pt();
    double eta_leading = (em.Pt() > ep.Pt()) ? em.Eta() : ep.Eta();
    double eta_subleading = (em.Pt() > ep.Pt()) ? ep.Eta() : em.Eta();

    hEtaVsPt_leading->Fill(pt_leading, eta_leading, weight);
    hEtaVsPt_subleading->Fill(pt_subleading, eta_subleading, weight);
    hCosthVsPt_leading->Fill(pt_leading, costheta, weight);
    hCosthVsPt_subleading->Fill(pt_subleading, costheta, weight);

    // Pairwise angular observables between the two leptons.
    const double dEta_ll = fabs(em.Eta() - ep.Eta());
    const double dPhi_ll = fabs(atan2(sin(em.Phi() - ep.Phi()), cos(em.Phi() - ep.Phi())));
    const double dR_ll = sqrt(dEta_ll*dEta_ll + dPhi_ll*dPhi_ll);
    const double cosOpening_ll = em.Vect().Unit().Dot(ep.Vect().Unit());

    hDeltaR_ll->Fill(dR_ll, weight);
    hDeltaPhi_ll->Fill(dPhi_ll, weight);
    hDeltaEta_ll->Fill(dEta_ll, weight);
    hDeltaEtaVsDeltaPhi_ll->Fill(dEta_ll, dPhi_ll, weight);
    hDeltaEtaVsCosth_ll->Fill(dEta_ll, costheta, weight);
    hDeltaEtaVsZPt_ll->Fill(dEta_ll, z.Pt(), weight);
    hDeltaEtaVsZY_ll->Fill(dEta_ll, fabs(z.Rapidity()), weight);
    hDeltaEtaVsLeadPt_ll->Fill(dEta_ll, pt_leading, weight);
    hDeltaEtaVsSubleadPt_ll->Fill(dEta_ll, pt_subleading, weight);
    hCosOpening_ll->Fill(cosOpening_ll, weight);
    hCosOpeningVsZPt_ll->Fill(z.Pt(), cosOpening_ll, weight);
    hCosOpeningVsLeadPt_ll->Fill(pt_leading, cosOpening_ll, weight);
    hCosOpeningVsSubleadPt_ll->Fill(pt_subleading, cosOpening_ll, weight);
    hCosOpeningVsDeltaEta_ll->Fill(dEta_ll, cosOpening_ll, weight);

    // Lepton-lepton eta/phi correlations
    hEta_ep_vs_em->Fill(ep.Eta(), em.Eta(), weight);
    hEta_eleading_vs_esubleading->Fill(abs(eta_leading), abs(eta_subleading), weight);
    hPhi_ep_vs_em->Fill(ep.Phi(), em.Phi(), weight);
    hPt_ep_vs_em->Fill(ep.Pt(), em.Pt(), weight);
    hPt_leading_vs_subleading->Fill(pt_leading, pt_subleading, weight);

    // Z–lepton correlations
    hZptVsEta_m->Fill(z.Pt(), em.Eta(), weight);
    hZptVsEta_p->Fill(z.Pt(), ep.Eta(), weight);
    hZYVsPt_m->Fill(fabs(z.Rapidity()), em.Pt(), weight);
    hZYVsPt_p->Fill(fabs(z.Rapidity()), ep.Pt(), weight);
    hZYVsPt_leading->Fill(fabs(z.Rapidity()), pt_leading, weight);
    hZYVsPt_subleading->Fill(fabs(z.Rapidity()), pt_subleading, weight);
    hZYVsEta_leading->Fill(fabs(z.Rapidity()), eta_leading, weight);
    hZYVsEta_subleading->Fill(fabs(z.Rapidity()), eta_subleading, weight);

    double eta_fwd = (fabs(em.Eta()) >= fabs(ep.Eta())) ? fabs(em.Eta()) : fabs(ep.Eta());
    double eta_cen = (fabs(em.Eta()) >= fabs(ep.Eta())) ? fabs(ep.Eta()) : fabs(em.Eta());
    double pt_fwd = (fabs(em.Eta()) >= fabs(ep.Eta())) ? em.Pt() : ep.Pt();
    double pt_cen = (fabs(em.Eta()) >= fabs(ep.Eta())) ? ep.Pt() : em.Pt();
    hZptVsEta_forward->Fill(z.Pt(), eta_fwd, weight);
    hZptVsEta_central->Fill(z.Pt(), eta_cen, weight);
    hZptVsPt_forward->Fill(z.Pt(), pt_fwd, weight);
    hZptVsPt_central->Fill(z.Pt(), pt_cen, weight);
    hZYVsEta_forward->Fill(fabs(z.Rapidity()), eta_fwd, weight);
    hZYVsEta_central->Fill(fabs(z.Rapidity()), eta_cen, weight);
    hZYVsPt_forward->Fill(fabs(z.Rapidity()), pt_fwd, weight);
    hZYVsPt_central->Fill(fabs(z.Rapidity()), pt_cen, weight);

    // Additional truth-level angular correlations
    if (ifTrueOnly) {
      // Truth-only ntuple: cosThetaCSTruth/phiCSTruth not loaded, use the one obtained from the lepton four-vectors instead
      hCosVspt_el->Fill(em.Pt(), costheta, weight);
      hCosVspt_pos->Fill(ep.Pt(), costheta, weight);
      hPhiVspt_el->Fill(em.Pt(), phi, weight);
      hPhiVspt_pos->Fill(ep.Pt(), phi, weight);
    } else {
      hCosVspt_el->Fill(em.Pt(), cosThetaCSTruth, weight);
      hCosVspt_pos->Fill(ep.Pt(), cosThetaCSTruth, weight);
      hPhiVspt_el->Fill(em.Pt(), phiCSTruthWrapped, weight);
      hPhiVspt_pos->Fill(ep.Pt(), phiCSTruthWrapped, weight);
    }

    // cosTheta slices by lepton pT
    if (em.Pt() >= 0. && em.Pt() < 5.)      hCosLneg0_5->Fill(costheta, weight);
    else if (em.Pt() < 20.)             hCosLneg5_20->Fill(costheta, weight);
    else if (em.Pt() < 40.)             hCosLneg20_40->Fill(costheta, weight);
    else if (em.Pt() < 80.)             hCosLneg40_80->Fill(costheta, weight);

    if (ep.Pt() >= 0. && ep.Pt() < 5.)      hCosLpos0_5->Fill(costheta, weight);
    else if (ep.Pt() < 20.)             hCosLpos5_20->Fill(costheta, weight);
    else if (ep.Pt() < 40.)             hCosLpos20_40->Fill(costheta, weight);
    else if (ep.Pt() < 80.)             hCosLpos40_80->Fill(costheta, weight);

    if(isY){
      Xs->Fill(fabs(z.Rapidity()));
      Xsw->Fill(fabs(z.Rapidity()), weight);
    } else {
      Xs->Fill(z.Pt());
      Xsw->Fill(z.Pt(), weight);
    }

    // Build the angular basis with TLVUtils so the polynomial definitions are
    // exactly shared with the common utility implementation.
    std::vector<double> aipols;
    TLVUtils::getAiPolynoms(costheta, phi, aipols);

    // At zero Z transverse momentum the hadron plane (and phi) is undefined.
    // Exclude it from the new diagnostics and count it explicitly.
    if (z.Pt() > 1.e-9 && std::isfinite(weight)) {
      const double p6 = aipols[6];
      const double signP6 = (p6 > 0.) ? 1. : ((p6 < 0.) ? -1. : 0.);
      hP6->Fill(p6, weight);
      hAbsP6->Fill(fabs(p6), weight);
      const double coordinates[] = {dEta_ll, z.Pt(), fabs(z.Rapidity()), M_PI-dPhi_ll};
      for (int j=0; j<4; ++j) {
        p6Correlations[j]->Fill(coordinates[j], p6, weight);
        p6Moments[j]->Fill(coordinates[j], 5.*p6, weight);
        p6Signs[j]->Fill(coordinates[j], signP6, weight);
        p6LeverArms[j]->Fill(coordinates[j], p6*p6, weight);
      }
      const bool negForward = fabs(em.Eta()) >= 2.5;
      const bool posForward = fabs(ep.Eta()) >= 2.5;
      if (negForward != posForward)
        hP6ForwardCharge->Fill(negForward ? -1. : 1., p6, weight);
    } else if (z.Pt() <= 1.e-9) {
      ++p6UndefinedPlane;
    }

    // Keep the historical AIZ histogram convention (scaled basis) for output
    // compatibility. The raw TLVUtils basis terms are transformed as follows:
    //   A0_scaled = 20/3 * p0 + 2/3
    //   A1_scaled = 5    * p1
    //   A2_scaled = 10   * p2
    //   A3_scaled = 4    * p3
    //   A4_scaled = 4    * p4
    //   A5_scaled = 5    * p5
    //   A6_scaled = 4    * p6 (legacy; standard A6 = 5<P6>)
    //   A7_scaled = 4    * p7
    // where p{i} == aipols[i] returned by TLVUtils::getAiPolynoms.
    const double a0Scaled = (20./3.) * aipols[0] + (2./3.);
    const double a1Scaled = 5. * aipols[1];
    const double a2Scaled = 10. * aipols[2];
    const double a3Scaled = 4. * aipols[3];
    const double a4Scaled = 4. * aipols[4];
    const double a5Scaled = 5. * aipols[5];
    const double a6Scaled = 4. * aipols[6];
    const double a7Scaled = 4. * aipols[7];

    if(isY){
      A0->Fill(fabs(z.Rapidity()), a0Scaled * weight);
      A1->Fill(fabs(z.Rapidity()), a1Scaled * weight);
      A2->Fill(fabs(z.Rapidity()), a2Scaled * weight);
      A3->Fill(fabs(z.Rapidity()), a3Scaled * weight);
      A4->Fill(fabs(z.Rapidity()), a4Scaled * weight);
      A5->Fill(fabs(z.Rapidity()), a5Scaled * weight);
      A6->Fill(fabs(z.Rapidity()), a6Scaled * weight);
      A7->Fill(fabs(z.Rapidity()), a7Scaled * weight);
    } else {
      A0->Fill(z.Pt(), a0Scaled * weight);
      A1->Fill(z.Pt(), a1Scaled * weight);
      A2->Fill(z.Pt(), a2Scaled * weight);
      A3->Fill(z.Pt(), a3Scaled * weight);
      A4->Fill(z.Pt(), a4Scaled * weight);
      A5->Fill(z.Pt(), a5Scaled * weight);
      A6->Fill(z.Pt(), a6Scaled * weight);
      A7->Fill(z.Pt(), a7Scaled * weight);
    }
  }

  Output->cd();
  AIZP6::Style(hP6); AIZP6::Style(hAbsP6); AIZP6::Style(hP6ForwardCharge);
  hP6ForwardCharge->GetXaxis()->SetBinLabel(1, "e^{-} forward");
  hP6ForwardCharge->GetXaxis()->SetBinLabel(2, "e^{+} forward");
  hP6->Write(); hAbsP6->Write(); hP6ForwardCharge->Write();
  TH1D p6Skipped("P6_undefinedPlane", ";Reason;Unweighted events", 1, 0., 1.);
  p6Skipped.SetDirectory(nullptr);
  p6Skipped.GetXaxis()->SetBinLabel(1, "pT(Z) <= 1e-9 GeV");
  p6Skipped.SetBinContent(1, p6UndefinedPlane);
  p6Skipped.Write();
  for (int j=0; j<4; ++j) {
    AIZP6::Style(p6Correlations[j]);
    AIZP6::Style(p6Moments[j]->mean);
    AIZP6::Style(p6Signs[j]->mean);
    AIZP6::Style(p6LeverArms[j]->mean);
    p6Correlations[j]->Write();
    p6Moments[j]->Write(); p6Signs[j]->Write(); p6LeverArms[j]->Write();
  }
  TString plotPrefix = nameOutput;
  plotPrefix.ReplaceAll(".root", "");
  TCanvas *cP6Polynomial = AIZPlotP6((plotPrefix+"_P6_polynomial.pdf").Data());
  cP6Polynomial->Write();
  TCanvas *cP6 = new TCanvas("c_P6_observables", "P6-sensitive observables", 1400, 1000);
  cP6->Divide(3,3);
  for (int pad=1; pad<=9; ++pad) {
    cP6->cd(pad)->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15); gPad->SetRightMargin(0.14);
  }
  cP6->cd(1); hP6->Draw("E");
  cP6->cd(2); hAbsP6->Draw("E");
  cP6->cd(3); hP6ForwardCharge->Draw("COLZ");
  for (int j=0; j<3; ++j) {
    cP6->cd(4+j); p6Correlations[j]->Draw("COLZ");
    cP6->cd(7+j); p6Moments[j]->mean->SetStats(false); p6Moments[j]->mean->Draw("E");
  }
  cP6->Write();
  cP6->SaveAs((plotPrefix+"_P6_observables.pdf").Data());
  TCanvas *cP6Sensitivity = new TCanvas("c_P6_sensitivity", "P6 asymmetry and lever arm", 1600, 800);
  cP6Sensitivity->Divide(4,2);
  for (int pad=1; pad<=8; ++pad) {
    cP6Sensitivity->cd(pad)->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.15);
  }
  for (int j=0; j<4; ++j) {
    cP6Sensitivity->cd(1+j); p6Signs[j]->mean->SetStats(false); p6Signs[j]->mean->Draw("E");
    cP6Sensitivity->cd(5+j); p6LeverArms[j]->mean->SetStats(false); p6LeverArms[j]->mean->Draw("E");
  }
  cP6Sensitivity->Write();
  cP6Sensitivity->SaveAs((plotPrefix+"_P6_sensitivity.pdf").Data());

  A0->Divide(Xsw);
  A1->Divide(Xsw);
  A2->Divide(Xsw);
  A3->Divide(Xsw);
  A4->Divide(Xsw);
  A5->Divide(Xsw);
  A6->Divide(Xsw);
  A7->Divide(Xsw);

  hctheta->Write();
  hctheta_truth->Write();
  hphi->Write();
  hphi_truth->Write();
  A0->Write();
  A1->Write();
  A2->Write();
  A3->Write();
  A4->Write();
  A5->Write();
  A6->Write();
  A7->Write();
  Xsw->Write("sigma");
  Xs->Write("sigma_unweighted");
  Zmass->Write();
  
  // Comparison plot: from lepton vs truth-level cosThetaCS
  TH1D *hctheta_norm = (TH1D*)hctheta->Clone("Costheta_norm");
  TH1D *hctheta_truth_norm = (TH1D*)hctheta_truth->Clone("CosthetaTruth_norm");
  if (hctheta_norm->Integral() != 0) hctheta_norm->Scale(1.0 / hctheta_norm->Integral());
  if (hctheta_truth_norm->Integral() != 0) hctheta_truth_norm->Scale(1.0 / hctheta_truth_norm->Integral());

  hctheta_norm->SetLineColor(kRed+1);
  hctheta_norm->SetLineWidth(2);
  hctheta_truth_norm->SetLineColor(kBlue+1);
  hctheta_truth_norm->SetLineStyle(2);
  hctheta_truth_norm->SetLineWidth(2);
  hctheta_norm->SetTitle("cos#theta_{CS}: from lepton vs truth;cos#theta_{CS};Normalized entries");

  TH1D *hctheta_ratio = (TH1D*)hctheta_norm->Clone("Costheta_ratio");
  hctheta_ratio->SetTitle(";cos#theta_{CS};From lepton / Truth");
  hctheta_ratio->Divide(hctheta_truth_norm);
  hctheta_ratio->SetLineColor(kBlack);
  hctheta_ratio->SetLineWidth(2);
  hctheta_ratio->SetMarkerStyle(20);
  hctheta_ratio->SetMarkerSize(0.7);

  TCanvas *c_costheta_comp = new TCanvas("c_costheta_comp", "cosThetaCS comparison", 800, 800);
  c_costheta_comp->Divide(1,2,0,0);

  // Top pad: normalized overlay
  TPad *pad_top = (TPad*)c_costheta_comp->cd(1);
  pad_top->SetPad(0.0, 0.30, 1.0, 1.0);
  pad_top->SetBottomMargin(0.02);
  hctheta_norm->Draw("hist");
  hctheta_truth_norm->Draw("hist same");
  TLegend *leg_costheta = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg_costheta->AddEntry(hctheta_norm, "from lepton costheta", "l");
  leg_costheta->AddEntry(hctheta_truth_norm, "truth costheta", "l");
  leg_costheta->Draw();

  // Bottom pad: ratio
  TPad *pad_bottom = (TPad*)c_costheta_comp->cd(2);
  pad_bottom->SetPad(0.0, 0.0, 1.0, 0.30);
  pad_bottom->SetTopMargin(0.05);
  pad_bottom->SetBottomMargin(0.30);
  hctheta_ratio->SetStats(false);
  hctheta_ratio->GetYaxis()->SetTitleOffset(0.5);
  hctheta_ratio->GetYaxis()->SetTitleSize(0.08);
  hctheta_ratio->GetYaxis()->SetLabelSize(0.08);
  hctheta_ratio->GetXaxis()->SetTitleSize(0.10);
  hctheta_ratio->GetXaxis()->SetLabelSize(0.10);
  hctheta_ratio->GetYaxis()->SetNdivisions(506);
  hctheta_ratio->Draw("ep");

  c_costheta_comp->Write();

  // Comparison plot: reconstructed vs truth phiCS
  TH1D *hphi_norm = (TH1D*)hphi->Clone("phi_norm");
  TH1D *hphi_truth_norm = (TH1D*)hphi_truth->Clone("phiTruth_norm");
  if (hphi_norm->Integral() != 0) hphi_norm->Scale(1.0 / hphi_norm->Integral());
  if (hphi_truth_norm->Integral() != 0) hphi_truth_norm->Scale(1.0 / hphi_truth_norm->Integral());

  hphi_norm->SetLineColor(kRed+1);
  hphi_norm->SetLineWidth(2);
  hphi_truth_norm->SetLineColor(kBlue+1);
  hphi_truth_norm->SetLineStyle(2);
  hphi_truth_norm->SetLineWidth(2);
  hphi_norm->SetTitle("#phi_{CS}: from lepton vs truth;#phi_{CS};Normalized entries");

  TH1D *hphi_ratio = (TH1D*)hphi_norm->Clone("phi_ratio");
  hphi_ratio->SetTitle(";#phi_{CS};From lepton / Truth");
  hphi_ratio->Divide(hphi_truth_norm);
  hphi_ratio->SetLineColor(kBlack);
  hphi_ratio->SetLineWidth(2);
  hphi_ratio->SetMarkerStyle(20);
  hphi_ratio->SetMarkerSize(0.7);

  TCanvas *c_phi_comp = new TCanvas("c_phi_comp", "phiCS comparison", 800, 800);
  c_phi_comp->Divide(1,2,0,0);

  // Top pad: normalized overlay
  TPad *pad_top_phi = (TPad*)c_phi_comp->cd(1);
  pad_top_phi->SetPad(0.0, 0.30, 1.0, 1.0);
  pad_top_phi->SetBottomMargin(0.02);
  hphi_norm->Draw("hist");
  hphi_truth_norm->Draw("hist same");
  TLegend *leg_phi = new TLegend(0.55, 0.72, 0.88, 0.88);
  leg_phi->AddEntry(hphi_norm, "from lepton phi", "l");
  leg_phi->AddEntry(hphi_truth_norm, "truth phi", "l");
  leg_phi->Draw();

  // Bottom pad: ratio
  TPad *pad_bottom_phi = (TPad*)c_phi_comp->cd(2);
  pad_bottom_phi->SetPad(0.0, 0.0, 1.0, 0.30);
  pad_bottom_phi->SetTopMargin(0.05);
  pad_bottom_phi->SetBottomMargin(0.30);
  hphi_ratio->SetStats(false);
  hphi_ratio->GetYaxis()->SetTitleOffset(0.5);
  hphi_ratio->GetYaxis()->SetTitleSize(0.08);
  hphi_ratio->GetYaxis()->SetLabelSize(0.08);
  hphi_ratio->GetXaxis()->SetTitleSize(0.10);
  hphi_ratio->GetXaxis()->SetLabelSize(0.10);
  hphi_ratio->GetYaxis()->SetNdivisions(506);
  hphi_ratio->Draw("ep");

  c_phi_comp->Write();

  // Write additional truth-level histograms
  hCosVspt_el->Write();
  hCosVspt_pos->Write();
  hPhiVspt_el->Write();
  hPhiVspt_pos->Write();
  hEta_ep_vs_em->Write();
  hEta_eleading_vs_esubleading->Write();
  hPhi_ep_vs_em->Write();
  hPt_ep_vs_em->Write();
  hPt_leading_vs_subleading->Write();
  hZptVsEta_m->Write();
  hZptVsEta_p->Write();
  hZptVsEta_forward->Write();
  hZptVsEta_central->Write();
  hZptVsPt_forward->Write();
  hZptVsPt_central->Write();
  hZYVsPt_m->Write();
  hZYVsPt_p->Write();
  hZYVsPt_leading->Write();
  hZYVsPt_subleading->Write();
  hZYVsPt_forward->Write();
  hZYVsPt_central->Write();
  hZYVsEta_leading->Write();
  hZYVsEta_subleading->Write();
  hZYVsEta_forward->Write();
  hZYVsEta_central->Write();
  hCosLneg0_5->Write();
  hCosLneg5_20->Write();
  hCosLneg20_40->Write();
  hCosLneg40_80->Write();
  hCosLpos0_5->Write();
  hCosLpos5_20->Write();
  hCosLpos20_40->Write();
  hCosLpos40_80->Write();
  hZptVsCostheta->Write();
  hZptVsPhi->Write();
  hZYVsCostheta->Write();
  hZYVsPhi->Write();
  hEtaVsPt_leading->Write();
  hEtaVsPt_subleading->Write();
  hCosthVsPt_leading->Write();
  hCosthVsPt_subleading->Write();
  hDeltaR_ll->Write();
  hDeltaPhi_ll->Write();
  hDeltaEta_ll->Write();
  hDeltaEtaVsDeltaPhi_ll->Write();
  hDeltaEtaVsCosth_ll->Write();
  hDeltaEtaVsZPt_ll->Write();
  hDeltaEtaVsZY_ll->Write();
  hDeltaEtaVsLeadPt_ll->Write();
  hDeltaEtaVsSubleadPt_ll->Write();
  hCosOpening_ll->Write();
  hCosOpeningVsZPt_ll->Write();
  hCosOpeningVsLeadPt_ll->Write();
  hCosOpeningVsSubleadPt_ll->Write();
  hCosOpeningVsDeltaEta_ll->Write();
  Output->Close();
}
