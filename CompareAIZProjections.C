#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include <iostream>
#include <algorithm>

// MAcro to post Processing ooutput from AIZ.C or AIW.C 

void CompareAIZProjections(const TString& inFile = "AI_Z_Truth_Zai_finalbinningPowheg_Y_NormXsec.root", bool isY = false) {
  TFile* f = TFile::Open(inFile, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "ERROR: cannot open file " << inFile << std::endl;
    return;
  }
  bool isFiducial = false;
 bool EtaOnly = false;
 bool CFonly = false;
 bool CFonlyEta = false;
 bool CCCF = false;

  if ( inFile.Contains("Fiducial") ) {
    isFiducial = true;
    CFonlyEta = inFile.Contains("CFonlyEta");
    EtaOnly = !CFonlyEta && inFile.Contains("EtaOnly");
    CFonly = !CFonlyEta && inFile.Contains("CFonly");
    CCCF = inFile.Contains("CCCF");
    std::cout << "Comparing projections from " << inFile
              << " with fiducial cuts applied. ETA only is " << EtaOnly
              << " CFonly is: " << CFonly
              << " CFonlyEta is: " << CFonlyEta
              << " CCCF is: " << CCCF << std::endl;
   } else {
    isFiducial = false;
  }
  
    std::cout << "Comparing projections from " << inFile << " with isY = " << isY <<  " Is Fiducial: " << isFiducial << std::endl;

  TH2D* hCosthLead = static_cast<TH2D*>(f->Get("costh_vs_pt_leading"));
  TH2D* hCosthSub  = static_cast<TH2D*>(f->Get("costh_vs_pt_subleading"));
  TH2D* hDEtaLead  = static_cast<TH2D*>(f->Get("deltaEta_vs_leadingPt_ll"));
  TH2D* hDEtaSub   = static_cast<TH2D*>(f->Get("deltaEta_vs_subleadingPt_ll"));
  TH2D* hDEtaCosth = static_cast<TH2D*>(f->Get("deltaEta_vs_costh_ll"));
  TH2D* hPtEpVsEm  = static_cast<TH2D*>(f->Get("pt_ep_vs_em"));
  TH2D* hEtaLeadVsSublead = static_cast<TH2D*>(f->Get("eta_eleading_vs_esubleading"));
  TH2D* hPtLeadVsSublead = static_cast<TH2D*>(f->Get("pt_leading_vs_subleading"));
  TH2D* hBosonVsObsM = static_cast<TH2D*>(f->Get(isY ? "zY_vs_pt_m" : "zPt_vs_eta_m"));
  TH2D* hBosonVsObsP = static_cast<TH2D*>(f->Get(isY ? "zY_vs_pt_p" : "zPt_vs_eta_p"));
  TH2D* hBosonVsCostheta = static_cast<TH2D*>(f->Get(isY ? "zY_vs_costheta" : "zPt_vs_costheta"));
  TH2D* hBosonVsPhi      = static_cast<TH2D*>(f->Get(isY ? "zY_vs_phi" : "zPt_vs_phi"));
  TH2D* hDEtaVsBoson     = static_cast<TH2D*>(f->Get(isY ? "deltaEta_vs_zY_ll" : "deltaEta_vs_zPt_ll"));

  // Optional inputs for additional Y(Z)-sliced leading/subleading overlays.
  TH2D* hZYVsPtLead  = static_cast<TH2D*>(f->Get("zY_vs_pt_leading"));
  TH2D* hZYVsPtSub   = static_cast<TH2D*>(f->Get("zY_vs_pt_subleading"));
  TH2D* hZYVsEtaLead = static_cast<TH2D*>(f->Get("zY_vs_eta_leading"));
  TH2D* hZYVsEtaSub  = static_cast<TH2D*>(f->Get("zY_vs_eta_subleading"));

  if (!hCosthLead || !hCosthSub || !hDEtaLead || !hDEtaSub || !hDEtaCosth || !hPtEpVsEm || !hPtLeadVsSublead || !hBosonVsObsM || !hBosonVsObsP || !hBosonVsCostheta || !hBosonVsPhi || !hDEtaVsBoson) {
    std::cout << "ERROR: one or more required histograms are missing in " << inFile << std::endl;
    std::cout << "Needed: costh_vs_pt_leading, costh_vs_pt_subleading, "
              << "deltaEta_vs_leadingPt_ll, deltaEta_vs_subleadingPt_ll, deltaEta_vs_costh_ll, pt_ep_vs_em, pt_leading_vs_subleading, "
              << (isY ? "zY_vs_pt_m, zY_vs_pt_p, zY_vs_costheta, zY_vs_phi, deltaEta_vs_zY_ll"
                       : "zPt_vs_eta_m, zPt_vs_eta_p, zPt_vs_costheta, zPt_vs_phi, deltaEta_vs_zPt_ll")
              << std::endl;
    f->Close();
    return;
  }

  gStyle->SetOptStat(0);

  TCanvas* cEtaLeadSublead2D = nullptr;
  TCanvas* cPtLeadSublead2D = nullptr;

  if (hEtaLeadVsSublead) {
    cEtaLeadSublead2D = new TCanvas("c_eta_eleading_vs_esubleading_2d", "|#eta(leading)| vs |#eta(subleading)|", 900, 800);
    cEtaLeadSublead2D->cd();
    cEtaLeadSublead2D->SetLogz();
    cEtaLeadSublead2D->SetGridx();
    cEtaLeadSublead2D->SetGridy();
    hEtaLeadVsSublead->SetTitle("|#eta(e_{leading})| vs |#eta(e_{subleading})|;|#eta(e_{leading})|;|#eta(e_{subleading})|");
    hEtaLeadVsSublead->Draw("colz");
    hEtaLeadVsSublead->SetStats(0);

 
  } else {
    std::cout << "WARNING: missing eta_eleading_vs_esubleading. Skipping 2D eta leading/subleading canvas." << std::endl;
  }

  cPtLeadSublead2D = new TCanvas("c_pt_leading_vs_subleading_2d", "p_{T}(leading) vs p_{T}(subleading)", 900, 800);
  cPtLeadSublead2D->cd();
  cPtLeadSublead2D->SetLogz();
  cPtLeadSublead2D->SetGridx();
  cPtLeadSublead2D->SetGridy();
  hPtLeadVsSublead->SetTitle("p_{T}(leading) vs p_{T}(subleading);p_{T}(leading l) [GeV];p_{T}(subleading l) [GeV]");
  hPtLeadVsSublead->Draw("colz");

  // 1) Compare cos(theta_CS) projections from leading/subleading maps.
  TH1D* pCosthLead = hCosthLead->ProjectionY("pCosth_leading");
  TH1D* pCosthSub  = hCosthSub->ProjectionY("pCosth_subleading");
  if (pCosthLead->Integral() > 0) pCosthLead->Scale(1.0 / pCosthLead->Integral());
  if (pCosthSub->Integral() > 0)  pCosthSub->Scale(1.0 / pCosthSub->Integral());

  pCosthLead->SetLineColor(kRed + 1);
  pCosthLead->SetLineWidth(2);
  pCosthSub->SetLineColor(kBlue + 1);
  pCosthSub->SetLineWidth(2);
  pCosthSub->SetLineStyle(2);

  TCanvas* cCosth = new TCanvas("c_compare_costh", "cos#theta_{CS} projections", 900, 800);
  TPad* padCosthTop = new TPad("padCosthTop", "padCosthTop", 0.0, 0.30, 1.0, 1.0);
  TPad* padCosthBot = new TPad("padCosthBot", "padCosthBot", 0.0, 0.00, 1.0, 0.30);
  padCosthTop->SetBottomMargin(0.02);
  padCosthBot->SetTopMargin(0.04);
  padCosthBot->SetBottomMargin(0.30);
  cCosth->cd();
  padCosthTop->Draw();
  padCosthBot->Draw();

  padCosthTop->cd();
  pCosthLead->SetTitle("Projection: cos#theta_{CS} from leading/subleading;p_{T} integrated cos#theta_{CS};Normalized entries");
  pCosthLead->Draw("hist");
  pCosthSub->Draw("hist same");

  TLegend* legCosth = new TLegend(0.55, 0.75, 0.88, 0.88);
  legCosth->AddEntry(pCosthLead, "from costh_vs_pt_leading", "l");
  legCosth->AddEntry(pCosthSub, "from costh_vs_pt_subleading", "l");
  legCosth->Draw();

  TH1D* rCosth = static_cast<TH1D*>(pCosthLead->Clone("ratio_costh_lead_over_sub"));
  rCosth->SetTitle(";cos#theta_{CS};Leading / Subleading");
  rCosth->Divide(pCosthSub);
  rCosth->SetLineColor(kBlack);
  rCosth->SetMarkerColor(kBlack);
  rCosth->SetMarkerStyle(20);
  rCosth->SetMarkerSize(0.7);
  rCosth->GetYaxis()->SetRangeUser(0.5, 1.5);

  padCosthBot->cd();
  rCosth->GetYaxis()->SetTitleSize(0.09);
  rCosth->GetYaxis()->SetLabelSize(0.08);
  rCosth->GetYaxis()->SetTitleOffset(0.55);
  rCosth->GetYaxis()->SetNdivisions(505);
  rCosth->GetXaxis()->SetTitleSize(0.11);
  rCosth->GetXaxis()->SetLabelSize(0.10);
  rCosth->Draw("ep");

  // 2) Compare |DeltaEta| projections from leading/subleading pT maps.
  TH1D* pDEtaLead = hDEtaLead->ProjectionX("pDeltaEta_leadingPt");
  TH1D* pDEtaSub  = hDEtaSub->ProjectionX("pDeltaEta_subleadingPt");
  if (pDEtaLead->Integral() > 0) pDEtaLead->Scale(1.0 / pDEtaLead->Integral());
  if (pDEtaSub->Integral() > 0)  pDEtaSub->Scale(1.0 / pDEtaSub->Integral());

  pDEtaLead->SetLineColor(kRed + 1);
  pDEtaLead->SetLineWidth(2);
  pDEtaSub->SetLineColor(kBlue + 1);
  pDEtaSub->SetLineWidth(2);
  pDEtaSub->SetLineStyle(2);

  TCanvas* cDEtaPt = new TCanvas("c_compare_deltaeta_pt", "|DeltaEta| projections", 900, 800);
  TPad* padDEtaTop = new TPad("padDEtaTop", "padDEtaTop", 0.0, 0.30, 1.0, 1.0);
  TPad* padDEtaBot = new TPad("padDEtaBot", "padDEtaBot", 0.0, 0.00, 1.0, 0.30);
  padDEtaTop->SetBottomMargin(0.02);
  padDEtaBot->SetTopMargin(0.04);
  padDEtaBot->SetBottomMargin(0.30);
  cDEtaPt->cd();
  padDEtaTop->Draw();
  padDEtaBot->Draw();

  padDEtaTop->cd();
  pDEtaLead->SetTitle("Projection: |#Delta#eta| from leading/subleading p_{T} maps;|#Delta#eta(l_{1},l_{2})|;Normalized entries");
  pDEtaLead->Draw("hist");
  pDEtaSub->Draw("hist same");

  TLegend* legDEta = new TLegend(0.55, 0.75, 0.88, 0.88);
  legDEta->AddEntry(pDEtaLead, "from deltaEta_vs_leadingPt_ll", "l");
  legDEta->AddEntry(pDEtaSub, "from deltaEta_vs_subleadingPt_ll", "l");
  legDEta->Draw();

  TH1D* rDEta = static_cast<TH1D*>(pDEtaLead->Clone("ratio_deltaeta_lead_over_sub"));
  rDEta->SetTitle(";|#Delta#eta(l_{1},l_{2})|;Leading / Subleading");
  rDEta->Divide(pDEtaSub);
  rDEta->SetLineColor(kBlack);
  rDEta->SetMarkerColor(kBlack);
  rDEta->SetMarkerStyle(20);
  rDEta->SetMarkerSize(0.7);
  rDEta->GetYaxis()->SetRangeUser(0.5, 1.5);

  padDEtaBot->cd();
  rDEta->GetYaxis()->SetTitleSize(0.09);
  rDEta->GetYaxis()->SetLabelSize(0.08);
  rDEta->GetYaxis()->SetTitleOffset(0.55);
  rDEta->GetYaxis()->SetNdivisions(505);
  rDEta->GetXaxis()->SetTitleSize(0.11);
  rDEta->GetXaxis()->SetLabelSize(0.10);
  rDEta->Draw("ep");

  // 3) Projections of deltaEta_vs_costh.
  TH1D* pDEtaFromCosth = hDEtaCosth->ProjectionX("pDeltaEta_from_deltaEta_vs_costh");
  TH1D* pCosthFromDEta = hDEtaCosth->ProjectionY("pCosth_from_deltaEta_vs_costh");
  if (pDEtaFromCosth->Integral() > 0) pDEtaFromCosth->Scale(1.0 / pDEtaFromCosth->Integral());
  if (pCosthFromDEta->Integral() > 0) pCosthFromDEta->Scale(1.0 / pCosthFromDEta->Integral());

  pDEtaFromCosth->SetLineColor(kGreen + 2);
  pDEtaFromCosth->SetLineWidth(2);
  pCosthFromDEta->SetLineColor(kMagenta + 1);
  pCosthFromDEta->SetLineWidth(2);

  TCanvas* cDEtaCosth = new TCanvas("c_projection_deltaeta_costh", "deltaEta_vs_costh projections", 1100, 500);
  cDEtaCosth->Divide(2, 1);

  cDEtaCosth->cd(1);
  pDEtaFromCosth->SetTitle("ProjectionX of deltaEta_vs_costh_ll;|#Delta#eta(l_{1},l_{2})|;Normalized entries");
  pDEtaFromCosth->Draw("hist");

  cDEtaCosth->cd(2);
  pCosthFromDEta->SetTitle("ProjectionY of deltaEta_vs_costh_ll;cos#theta_{CS};Normalized entries");
  pCosthFromDEta->Draw("hist");

  // 4) pT-sliced projections in absolute bins for leading/subleading.
  // Requested bins: 0-5, 5-10, 10-20, 20-30, 30-40, 50-30000 GeV.
  const int nPtSlices = 6;
  const double ptEdgesLow[nPtSlices]  = {0., 5., 10., 20., 30., 50.};
  const double ptEdgesHigh[nPtSlices] = {5., 10., 20., 30., 40., 30000.};

  TCanvas* cCosthPtSlices = new TCanvas("c_compare_costh_pt_slices", "cos#theta_{CS} projections in p_{T} slices", 1400, 900);
  cCosthPtSlices->Divide(3, 2);

  TCanvas* cDEtaPtSlices = new TCanvas("c_compare_deltaeta_pt_slices", "|#Delta#eta| projections in p_{T} slices", 1400, 900);
  cDEtaPtSlices->Divide(3, 2);

  for (int i = 0; i < nPtSlices; ++i) {
    const double ptLo = ptEdgesLow[i];
    const double ptHi = ptEdgesHigh[i];

    // For costh_vs_pt: pT is X axis, so slice in X and project Y.
    int xBinLoCosth = hCosthLead->GetXaxis()->FindBin(ptLo + 1e-6);
    int xBinHiCosth = (ptHi >= hCosthLead->GetXaxis()->GetXmax())
                        ? hCosthLead->GetXaxis()->GetNbins()
                        : hCosthLead->GetXaxis()->FindBin(ptHi - 1e-6);
    xBinLoCosth = std::max(1, xBinLoCosth);
    xBinHiCosth = std::max(xBinLoCosth, std::min(hCosthLead->GetXaxis()->GetNbins(), xBinHiCosth));

    TH1D* pCosthLeadSlice = hCosthLead->ProjectionY(Form("pCosth_lead_ptbin_%d", i), xBinLoCosth, xBinHiCosth);
    TH1D* pCosthSubSlice  = hCosthSub->ProjectionY(Form("pCosth_sub_ptbin_%d", i), xBinLoCosth, xBinHiCosth);

    if (pCosthLeadSlice->Integral() > 0) pCosthLeadSlice->Scale(1.0 / pCosthLeadSlice->Integral());
    if (pCosthSubSlice->Integral() > 0)  pCosthSubSlice->Scale(1.0 / pCosthSubSlice->Integral());

    pCosthLeadSlice->SetLineColor(kRed + 1);
    pCosthLeadSlice->SetLineWidth(2);
    pCosthSubSlice->SetLineColor(kBlue + 1);
    pCosthSubSlice->SetLineWidth(2);
    pCosthSubSlice->SetLineStyle(2);

    cCosthPtSlices->cd(i + 1);
    pCosthLeadSlice->SetTitle(Form("%.0f < p_{T} < %.0f GeV;cos#theta_{CS};Normalized entries", ptLo, ptHi));
    pCosthLeadSlice->Draw("hist");
    pCosthSubSlice->Draw("hist same");
    TLegend* legCosthSlice = new TLegend(0.48, 0.74, 0.88, 0.88);
    legCosthSlice->AddEntry(pCosthLeadSlice, "leading", "l");
    legCosthSlice->AddEntry(pCosthSubSlice, "subleading", "l");
    legCosthSlice->Draw();

    // For deltaEta_vs_pt: pT is Y axis, so slice in Y and project X.
    int yBinLoDEta = hDEtaLead->GetYaxis()->FindBin(ptLo + 1e-6);
    int yBinHiDEta = (ptHi >= hDEtaLead->GetYaxis()->GetXmax())
                       ? hDEtaLead->GetYaxis()->GetNbins()
                       : hDEtaLead->GetYaxis()->FindBin(ptHi - 1e-6);
    yBinLoDEta = std::max(1, yBinLoDEta);
    yBinHiDEta = std::max(yBinLoDEta, std::min(hDEtaLead->GetYaxis()->GetNbins(), yBinHiDEta));

    TH1D* pDEtaLeadSlice = hDEtaLead->ProjectionX(Form("pDeltaEta_lead_ptbin_%d", i), yBinLoDEta, yBinHiDEta);
    TH1D* pDEtaSubSlice  = hDEtaSub->ProjectionX(Form("pDeltaEta_sub_ptbin_%d", i), yBinLoDEta, yBinHiDEta);

    if (pDEtaLeadSlice->Integral() > 0) pDEtaLeadSlice->Scale(1.0 / pDEtaLeadSlice->Integral());
    if (pDEtaSubSlice->Integral() > 0)  pDEtaSubSlice->Scale(1.0 / pDEtaSubSlice->Integral());

    pDEtaLeadSlice->SetLineColor(kRed + 1);
    pDEtaLeadSlice->SetLineWidth(2);
    pDEtaSubSlice->SetLineColor(kBlue + 1);
    pDEtaSubSlice->SetLineWidth(2);
    pDEtaSubSlice->SetLineStyle(2);

    cDEtaPtSlices->cd(i + 1);
    pDEtaLeadSlice->SetTitle(Form("%.0f < p_{T} < %.0f GeV;|#Delta#eta(l_{1},l_{2})|;Normalized entries", ptLo, ptHi));
    pDEtaLeadSlice->Draw("hist");
    pDEtaSubSlice->Draw("hist same");
    TLegend* legDEtaSlice = new TLegend(0.48, 0.74, 0.88, 0.88);
    legDEtaSlice->AddEntry(pDEtaLeadSlice, "leading", "l");
    legDEtaSlice->AddEntry(pDEtaSubSlice, "subleading", "l");
    legDEtaSlice->Draw();
  }

  // 5) X/Y sliced projections for deltaEta_vs_costh in custom bins.
  // delta_eta bins: [0,2,4,6,8,100]
  // cos(theta) bins: [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]
  const int nDEtaEdges = 8;
  const double dEtaEdges[nDEtaEdges] = {0., 0.5,1 ,1.5, 2., 2.5, 4.,100.};
  const int nCosEdges = 11;
  const double cosEdges[nCosEdges] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.};

  // ProjectionY (cos(theta)) in delta_eta slices (X-binned).
  TCanvas* cCosthInDEtaSlices = new TCanvas("c_projectionY_costh_in_deltaeta_slices", "ProjectionY in delta_eta slices", 1400, 900);
  cCosthInDEtaSlices->Divide(3, 2);
  for (int i = 0; i < nDEtaEdges - 1; ++i) {
    const double xLo = dEtaEdges[i];
    const double xHi = dEtaEdges[i + 1];
    int xBinLo = hDEtaCosth->GetXaxis()->FindBin(xLo + 1e-6);
    int xBinHi = (xHi >= hDEtaCosth->GetXaxis()->GetXmax())
                   ? hDEtaCosth->GetXaxis()->GetNbins()
                   : hDEtaCosth->GetXaxis()->FindBin(xHi - 1e-6);
    xBinLo = std::max(1, xBinLo);
    xBinHi = std::max(xBinLo, std::min(hDEtaCosth->GetXaxis()->GetNbins(), xBinHi));

    TH1D* pCosthSlice = hDEtaCosth->ProjectionY(Form("pCosth_dEtaSlice_%d", i), xBinLo, xBinHi);
    if (pCosthSlice->Integral() > 0) pCosthSlice->Scale(1.0 / pCosthSlice->Integral());
    pCosthSlice->SetLineColor(kAzure + 2);
    pCosthSlice->SetLineWidth(2);

    cCosthInDEtaSlices->cd(i + 1);
    pCosthSlice->SetTitle(Form("%.1f < |#Delta#eta| < %.1f;cos#theta_{CS};Normalized entries", xLo, xHi));
    pCosthSlice->Draw("hist");
  }

  // ProjectionX (delta_eta) in cos(theta) slices (Y-binned).
  TCanvas* cDEtaInCosthSlices = new TCanvas("c_projectionX_deltaeta_in_costh_slices", "ProjectionX in cos(theta) slices", 1600, 900);
  cDEtaInCosthSlices->Divide(4, 3);
  for (int i = 0; i < nCosEdges - 1; ++i) {
    const double yLo = cosEdges[i];
    const double yHi = cosEdges[i + 1];
    int yBinLo = hDEtaCosth->GetYaxis()->FindBin(yLo + 1e-6);
    int yBinHi = (yHi >= hDEtaCosth->GetYaxis()->GetXmax())
                   ? hDEtaCosth->GetYaxis()->GetNbins()
                   : hDEtaCosth->GetYaxis()->FindBin(yHi - 1e-6);
    yBinLo = std::max(1, yBinLo);
    yBinHi = std::max(yBinLo, std::min(hDEtaCosth->GetYaxis()->GetNbins(), yBinHi));

    TH1D* pDEtaSlice = hDEtaCosth->ProjectionX(Form("pDEta_costhSlice_%d", i), yBinLo, yBinHi);
    if (pDEtaSlice->Integral() > 0) pDEtaSlice->Scale(1.0 / pDEtaSlice->Integral());
    pDEtaSlice->SetLineColor(kOrange + 7);
    pDEtaSlice->SetLineWidth(2);

    cDEtaInCosthSlices->cd(i + 1);
    pDEtaSlice->SetTitle(Form("%.1f < cos#theta_{CS} < %.1f;|#Delta#eta(l_{1},l_{2})|;Normalized entries", yLo, yHi));
    pDEtaSlice->Draw("hist");
  }

  // 6) pt_ep_vs_em projections in pT bins.
  // Project pT(e-) for bins of pT(e+) and vice versa.
  TCanvas* cPtEpSlices = new TCanvas("c_projection_pT_em_in_pT_ep_slices", "p_{T}(e^{-}) projections in p_{T}(e^{+}) slices", 1400, 900);
  cPtEpSlices->Divide(3, 2);

  TCanvas* cPtEmSlices = new TCanvas("c_projection_pT_ep_in_pT_em_slices", "p_{T}(e^{+}) projections in p_{T}(e^{-}) slices", 1400, 900);
  cPtEmSlices->Divide(3, 2);

  for (int i = 0; i < nPtSlices; ++i) {
    const double ptLo = ptEdgesLow[i];
    const double ptHi = ptEdgesHigh[i];

    // For pt_ep_vs_em: pT(e+) is X axis, pT(e-) is Y axis
    // ProjectionY (pT of e-) for bins of pT(e+)
    int xBinLoPt = hPtEpVsEm->GetXaxis()->FindBin(ptLo + 1e-6);
    int xBinHiPt = (ptHi >= hPtEpVsEm->GetXaxis()->GetXmax())
                     ? hPtEpVsEm->GetXaxis()->GetNbins()
                     : hPtEpVsEm->GetXaxis()->FindBin(ptHi - 1e-6);
    xBinLoPt = std::max(1, xBinLoPt);
    xBinHiPt = std::max(xBinLoPt, std::min(hPtEpVsEm->GetXaxis()->GetNbins(), xBinHiPt));

    TH1D* pPtEmSlice = hPtEpVsEm->ProjectionY(Form("pPt_em_epbin_%d", i), xBinLoPt, xBinHiPt);
    if (pPtEmSlice->Integral() > 0) pPtEmSlice->Scale(1.0 / pPtEmSlice->Integral());
    pPtEmSlice->SetLineColor(kViolet - 5);
    pPtEmSlice->SetLineWidth(2);

    cPtEpSlices->cd(i + 1);
    pPtEmSlice->SetTitle(Form("%.0f < p_{T}(e^{+}) < %.0f GeV;p_{T}(e^{-}) [GeV];Normalized entries", ptLo, ptHi));
    pPtEmSlice->Draw("hist");

    // ProjectionX (pT of e+) for bins of pT(e-)
    int yBinLoPt = hPtEpVsEm->GetYaxis()->FindBin(ptLo + 1e-6);
    int yBinHiPt = (ptHi >= hPtEpVsEm->GetYaxis()->GetXmax())
                     ? hPtEpVsEm->GetYaxis()->GetNbins()
                     : hPtEpVsEm->GetYaxis()->FindBin(ptHi - 1e-6);
    yBinLoPt = std::max(1, yBinLoPt);
    yBinHiPt = std::max(yBinLoPt, std::min(hPtEpVsEm->GetYaxis()->GetNbins(), yBinHiPt));

    TH1D* pPtEpSlice = hPtEpVsEm->ProjectionX(Form("pPt_ep_embin_%d", i), yBinLoPt, yBinHiPt);
    if (pPtEpSlice->Integral() > 0) pPtEpSlice->Scale(1.0 / pPtEpSlice->Integral());
    pPtEpSlice->SetLineColor(kCyan + 2);
    pPtEpSlice->SetLineWidth(2);

    cPtEmSlices->cd(i + 1);
    pPtEpSlice->SetTitle(Form("%.0f < p_{T}(e^{-}) < %.0f GeV;p_{T}(e^{+}) [GeV];Normalized entries", ptLo, ptHi));
    pPtEpSlice->Draw("hist");
  }

  // 7) pt_leading_vs_subleading projections in pT bins.
  // Project pT(subleading) for bins of pT(leading) and vice versa.
  TCanvas* cPtLeadSlices = new TCanvas("c_projection_pT_sublead_in_pT_lead_slices", "p_{T}(subleading) projections in p_{T}(leading) slices", 1400, 900);
  cPtLeadSlices->Divide(3, 2);

  TCanvas* cPtSubleadSlices = new TCanvas("c_projection_pT_lead_in_pT_sublead_slices", "p_{T}(leading) projections in p_{T}(subleading) slices", 1400, 900);
  cPtSubleadSlices->Divide(3, 2);

  for (int i = 0; i < nPtSlices; ++i) {
    const double ptLo = ptEdgesLow[i];
    const double ptHi = ptEdgesHigh[i];

    // For pt_leading_vs_subleading: pT(leading) is X axis, pT(subleading) is Y axis
    // ProjectionY (pT of subleading) for bins of pT(leading)
    int xBinLoLead = hPtLeadVsSublead->GetXaxis()->FindBin(ptLo + 1e-6);
    int xBinHiLead = (ptHi >= hPtLeadVsSublead->GetXaxis()->GetXmax())
                       ? hPtLeadVsSublead->GetXaxis()->GetNbins()
                       : hPtLeadVsSublead->GetXaxis()->FindBin(ptHi - 1e-6);
    xBinLoLead = std::max(1, xBinLoLead);
    xBinHiLead = std::max(xBinLoLead, std::min(hPtLeadVsSublead->GetXaxis()->GetNbins(), xBinHiLead));

    TH1D* pPtSubleadSlice = hPtLeadVsSublead->ProjectionY(Form("pPt_sublead_leadbin_%d", i), xBinLoLead, xBinHiLead);
    if (pPtSubleadSlice->Integral() > 0) pPtSubleadSlice->Scale(1.0 / pPtSubleadSlice->Integral());
    pPtSubleadSlice->SetLineColor(kSpring + 5);
    pPtSubleadSlice->SetLineWidth(2);

    cPtLeadSlices->cd(i + 1);
    pPtSubleadSlice->SetTitle(Form("%.0f < p_{T}(lead) < %.0f GeV;p_{T}(sublead) [GeV];Normalized entries", ptLo, ptHi));
    pPtSubleadSlice->Draw("hist");

    // ProjectionX (pT of leading) for bins of pT(subleading)
    int yBinLoSublead = hPtLeadVsSublead->GetYaxis()->FindBin(ptLo + 1e-6);
    int yBinHiSublead = (ptHi >= hPtLeadVsSublead->GetYaxis()->GetXmax())
                         ? hPtLeadVsSublead->GetYaxis()->GetNbins()
                         : hPtLeadVsSublead->GetYaxis()->FindBin(ptHi - 1e-6);
    yBinLoSublead = std::max(1, yBinLoSublead);
    yBinHiSublead = std::max(yBinLoSublead, std::min(hPtLeadVsSublead->GetYaxis()->GetNbins(), yBinHiSublead));

    TH1D* pPtLeadSlice = hPtLeadVsSublead->ProjectionX(Form("pPt_lead_subleadbin_%d", i), yBinLoSublead, yBinHiSublead);
    if (pPtLeadSlice->Integral() > 0) pPtLeadSlice->Scale(1.0 / pPtLeadSlice->Integral());
    pPtLeadSlice->SetLineColor(kTeal + 2);
    pPtLeadSlice->SetLineWidth(2);

    cPtSubleadSlices->cd(i + 1);
    pPtLeadSlice->SetTitle(Form("%.0f < p_{T}(sublead) < %.0f GeV;p_{T}(lead) [GeV];Normalized entries", ptLo, ptHi));
    pPtLeadSlice->Draw("hist");
  }
  // 8) pt_ep_vs_em overlay: compare pT(e+) and pT(e-) distributions in pT bins with ratio panels.
  TCanvas* cPtEpEmOverlay = new TCanvas("c_overlay_pT_ep_vs_em_in_pT_slices", "p_{T}(e^{+}) vs p_{T}(e^{-}) overlay in p_{T} slices", 1000, 1500);
  cPtEpEmOverlay->Divide(2, 6);

  for (int i = 0; i < nPtSlices; ++i) {
    const double ptLo = ptEdgesLow[i];
    const double ptHi = ptEdgesHigh[i];

    // For pt_ep_vs_em: pT(e+) is X axis, pT(e-) is Y axis
    // ProjectionY (pT of e-) for bins of pT(e+)
    int xBinLoPt = hPtEpVsEm->GetXaxis()->FindBin(ptLo + 1e-6);
    int xBinHiPt = (ptHi >= hPtEpVsEm->GetXaxis()->GetXmax())
                     ? hPtEpVsEm->GetXaxis()->GetNbins()
                     : hPtEpVsEm->GetXaxis()->FindBin(ptHi - 1e-6);
    xBinLoPt = std::max(1, xBinLoPt);
    xBinHiPt = std::max(xBinLoPt, std::min(hPtEpVsEm->GetXaxis()->GetNbins(), xBinHiPt));

    TH1D* pPtEmSliceOvl = hPtEpVsEm->ProjectionY(Form("pPt_em_epbin_ovl_%d", i), xBinLoPt, xBinHiPt);
    if (pPtEmSliceOvl->Integral() > 0) pPtEmSliceOvl->Scale(1.0 / pPtEmSliceOvl->Integral());
    pPtEmSliceOvl->SetLineColor(kViolet - 5);
    pPtEmSliceOvl->SetLineWidth(2);

    // ProjectionX (pT of e+) for bins of pT(e-)
    int yBinLoPt = hPtEpVsEm->GetYaxis()->FindBin(ptLo + 1e-6);
    int yBinHiPt = (ptHi >= hPtEpVsEm->GetYaxis()->GetXmax())
                     ? hPtEpVsEm->GetYaxis()->GetNbins()
                     : hPtEpVsEm->GetYaxis()->FindBin(ptHi - 1e-6);
    yBinLoPt = std::max(1, yBinLoPt);
    yBinHiPt = std::max(yBinLoPt, std::min(hPtEpVsEm->GetYaxis()->GetNbins(), yBinHiPt));

    TH1D* pPtEpSliceOvl = hPtEpVsEm->ProjectionX(Form("pPt_ep_embin_ovl_%d", i), yBinLoPt, yBinHiPt);
    if (pPtEpSliceOvl->Integral() > 0) pPtEpSliceOvl->Scale(1.0 / pPtEpSliceOvl->Integral());
    pPtEpSliceOvl->SetLineColor(kCyan + 2);
    pPtEpSliceOvl->SetLineWidth(2);
    pPtEpSliceOvl->SetLineStyle(2);

    // Top pad: overlay distributions
    cPtEpEmOverlay->cd(2*i + 1);
    pPtEmSliceOvl->SetTitle(Form("%.0f < p_{T} < %.0f GeV;p_{T} [GeV];Normalized entries", ptLo, ptHi));
    pPtEmSliceOvl->Draw("hist");
    pPtEpSliceOvl->Draw("hist same");
    
    TLegend* legPtOverlay = new TLegend(0.48, 0.74, 0.88, 0.88);
    legPtOverlay->AddEntry(pPtEmSliceOvl, "p_{T}(e^{-}) from e^{+} bin", "l");
    legPtOverlay->AddEntry(pPtEpSliceOvl, "p_{T}(e^{+}) from e^{-} bin", "l");
    legPtOverlay->Draw();
    
    // Bottom pad: ratio
    cPtEpEmOverlay->cd(2*i + 2);
    TH1D* rPtOverlay = static_cast<TH1D*>(pPtEmSliceOvl->Clone(Form("ratio_pt_em_over_ep_%d", i)));
    rPtOverlay->SetTitle(";p_{T} [GeV];p_{T}(e^{-}) / p_{T}(e^{+})");
    rPtOverlay->Divide(pPtEpSliceOvl);
    rPtOverlay->SetLineColor(kBlack);
    rPtOverlay->SetMarkerColor(kBlack);
    rPtOverlay->SetMarkerStyle(20);
    rPtOverlay->SetMarkerSize(0.6);
    rPtOverlay->GetYaxis()->SetRangeUser(0.5, 1.5);
    rPtOverlay->GetYaxis()->SetTitleSize(0.10);
    rPtOverlay->GetYaxis()->SetLabelSize(0.09);
    rPtOverlay->GetYaxis()->SetTitleOffset(0.45);
    rPtOverlay->GetXaxis()->SetTitleSize(0.10);
    rPtOverlay->GetXaxis()->SetLabelSize(0.09);
    rPtOverlay->Draw("ep");
  }

  // 9) pt_leading_vs_subleading overlay: compare pT(leading) and pT(subleading) distributions in pT bins.
  TCanvas* cPtLeadSubleadOverlay = new TCanvas("c_overlay_pT_lead_vs_sublead_in_pT_slices", "p_{T}(leading) vs p_{T}(subleading) overlay in p_{T} slices", 1400, 900);
  cPtLeadSubleadOverlay->Divide(3, 2);

  for (int i = 0; i < nPtSlices; ++i) {
    const double ptLo = ptEdgesLow[i];
    const double ptHi = ptEdgesHigh[i];

    // For pt_leading_vs_subleading: pT(leading) is X axis, pT(subleading) is Y axis
    // ProjectionY (pT of subleading) for bins of pT(leading)
    int xBinLoLead = hPtLeadVsSublead->GetXaxis()->FindBin(ptLo + 1e-6);
    int xBinHiLead = (ptHi >= hPtLeadVsSublead->GetXaxis()->GetXmax())
                       ? hPtLeadVsSublead->GetXaxis()->GetNbins()
                       : hPtLeadVsSublead->GetXaxis()->FindBin(ptHi - 1e-6);
    xBinLoLead = std::max(1, xBinLoLead);
    xBinHiLead = std::max(xBinLoLead, std::min(hPtLeadVsSublead->GetXaxis()->GetNbins(), xBinHiLead));

    TH1D* pPtSubleadOvl = hPtLeadVsSublead->ProjectionY(Form("pPt_sublead_leadbin_ovl_%d", i), xBinLoLead, xBinHiLead);
    if (pPtSubleadOvl->Integral() > 0) pPtSubleadOvl->Scale(1.0 / pPtSubleadOvl->Integral());
    pPtSubleadOvl->SetLineColor(kSpring + 5);
    pPtSubleadOvl->SetLineWidth(2);

    // ProjectionX (pT of leading) for bins of pT(subleading)
    int yBinLoSublead = hPtLeadVsSublead->GetYaxis()->FindBin(ptLo + 1e-6);
    int yBinHiSublead = (ptHi >= hPtLeadVsSublead->GetYaxis()->GetXmax())
                         ? hPtLeadVsSublead->GetYaxis()->GetNbins()
                         : hPtLeadVsSublead->GetYaxis()->FindBin(ptHi - 1e-6);
    yBinLoSublead = std::max(1, yBinLoSublead);
    yBinHiSublead = std::max(yBinLoSublead, std::min(hPtLeadVsSublead->GetYaxis()->GetNbins(), yBinHiSublead));

    TH1D* pPtLeadOvl = hPtLeadVsSublead->ProjectionX(Form("pPt_lead_subleadbin_ovl_%d", i), yBinLoSublead, yBinHiSublead);
    if (pPtLeadOvl->Integral() > 0) pPtLeadOvl->Scale(1.0 / pPtLeadOvl->Integral());
    pPtLeadOvl->SetLineColor(kTeal + 2);
    pPtLeadOvl->SetLineWidth(2);
    //pPtLeadOvl->SetLineStyle(2);

    cPtLeadSubleadOverlay->cd(i + 1);
    pPtSubleadOvl->SetTitle(Form("%.0f < p_{T} < %.0f GeV;p_{T} [GeV];Normalized entries", ptLo, ptHi));
    pPtSubleadOvl->Draw("hist");
    pPtLeadOvl->Draw("hist same");
    
    TLegend* legPtLeadSublead = new TLegend(0.48, 0.74, 0.88, 0.88);
    legPtLeadSublead->AddEntry(pPtSubleadOvl, "p_{T}(sublead) from lead bin", "l");
    legPtLeadSublead->AddEntry(pPtLeadOvl, "p_{T}(lead) from sublead bin", "l");
    legPtLeadSublead->Draw();
  }

  // 10) Boson-variable sliced projections.
  // isY = false: pT(Z) bins [GeV] = [0,5,10,15,35,55,80,2000]
  // isY = true : |y(Z)| bins      = [0,0.4,0.8,1.6,2.4,3.2,4,8]
  const int nBosonSlices = 7;
  const double bosonLow[nBosonSlices]  = {0., 5., 10., 15., 35., 55., 80.};
  const double bosonHigh[nBosonSlices] = {5., 10., 15., 35., 55., 80., 2000.};
  const double yLow[nBosonSlices]      = {0.0, 0.4, 0.8, 1.6, 2.4, 3.2, 4.0};
  const double yHigh[nBosonSlices]     = {0.4, 0.8, 1.6, 2.4, 3.2, 4.0, 8.0};

  const char* bosonLabel = isY ? "|y(Z)|" : "p_{T}(Z)";
  const char* bosonUnit = isY ? "" : " GeV";
  const char* obsMinusLabel = isY ? "p_{T}(e^{-}) [GeV]" : "#eta(e^{-})";
  const char* obsPlusLabel = isY ? "p_{T}(e^{+}) [GeV]" : "#eta(e^{+})";

  TCanvas* cPtLeadSubleadInYSlices = nullptr;
  TCanvas* cEtaLeadSubleadInYSlices = nullptr;

  TCanvas* cBosonCosthSlices = new TCanvas(
    isY ? "c_projection_costheta_in_zY_slices" : "c_projection_costheta_in_zpt_slices",
    isY ? "cos#theta_{CS} in |y(Z)| slices" : "cos#theta_{CS} in p_{T}(Z) slices",
    1600, 900
  );
  cBosonCosthSlices->Divide(4, 2);

  TCanvas* cBosonPhiSlices = new TCanvas(
    isY ? "c_projection_phi_in_zY_slices" : "c_projection_phi_in_zpt_slices",
    isY ? "#phi_{CS} in |y(Z)| slices" : "#phi_{CS} in p_{T}(Z) slices",
    1600, 900
  );
  cBosonPhiSlices->Divide(4, 2);

  TCanvas* cBosonObsMSlices = new TCanvas(
    isY ? "c_projection_pt_m_in_zY_slices" : "c_projection_eta_m_in_zpt_slices",
    isY ? "p_{T}(e^{-}) in |y(Z)| slices" : "#eta(e^{-}) in p_{T}(Z) slices",
    1600, 900
  );
  cBosonObsMSlices->Divide(4, 2);

  TCanvas* cBosonObsPSlices = new TCanvas(
    isY ? "c_projection_pt_p_in_zY_slices" : "c_projection_eta_p_in_zpt_slices",
    isY ? "p_{T}(e^{+}) in |y(Z)| slices" : "#eta(e^{+}) in p_{T}(Z) slices",
    1600, 900
  );
  cBosonObsPSlices->Divide(4, 2);

  TCanvas* cDEtaInBosonSlices = new TCanvas(
    isY ? "c_projection_deltaeta_in_zY_slices" : "c_projection_deltaeta_in_zpt_slices",
    isY ? "|#Delta#eta| in |y(Z)| slices" : "|#Delta#eta| in p_{T}(Z) slices",
    1600, 900
  );
  cDEtaInBosonSlices->Divide(4, 2);

  TH1D* pDEtaBosonSlices[nBosonSlices] = {nullptr};
  const int sliceColors[nBosonSlices] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kAzure + 2, kBlack};

  for (int i = 0; i < nBosonSlices; ++i) {
    const double bosonLo = isY ? yLow[i] : bosonLow[i];
    const double bosonHi = isY ? yHigh[i] : bosonHigh[i];

    int xBinLoBoson = hBosonVsCostheta->GetXaxis()->FindBin(bosonLo + 1e-6);
    int xBinHiBoson = (bosonHi >= hBosonVsCostheta->GetXaxis()->GetXmax())
                        ? hBosonVsCostheta->GetXaxis()->GetNbins()
                        : hBosonVsCostheta->GetXaxis()->FindBin(bosonHi - 1e-6);
    xBinLoBoson = std::max(1, xBinLoBoson);
    xBinHiBoson = std::max(xBinLoBoson, std::min(hBosonVsCostheta->GetXaxis()->GetNbins(), xBinHiBoson));

    TH1D* pCosthInBoson = hBosonVsCostheta->ProjectionY(Form("pCostheta_bosonbin_%d", i), xBinLoBoson, xBinHiBoson);
    if (pCosthInBoson->Integral() > 0) pCosthInBoson->Scale(1.0 / pCosthInBoson->Integral());
    pCosthInBoson->SetLineColor(kBlue + 1);
    pCosthInBoson->SetLineWidth(2);

    cBosonCosthSlices->cd(i + 1);
    pCosthInBoson->SetTitle(Form("%.1f < %s < %.1f%s;cos#theta_{CS};Normalized entries", bosonLo, bosonLabel, bosonHi, bosonUnit));
    pCosthInBoson->Draw("hist");

    TH1D* pPhiInBoson = hBosonVsPhi->ProjectionY(Form("pPhi_bosonbin_%d", i), xBinLoBoson, xBinHiBoson);
    if (pPhiInBoson->Integral() > 0) pPhiInBoson->Scale(1.0 / pPhiInBoson->Integral());
    pPhiInBoson->SetLineColor(kRed + 1);
    pPhiInBoson->SetLineWidth(2);

    cBosonPhiSlices->cd(i + 1);
    pPhiInBoson->SetTitle(Form("%.1f < %s < %.1f%s;#phi_{CS};Normalized entries", bosonLo, bosonLabel, bosonHi, bosonUnit));
    pPhiInBoson->Draw("hist");

    TH1D* pObsMInBoson = hBosonVsObsM->ProjectionY(Form("pObsM_bosonbin_%d", i), xBinLoBoson, xBinHiBoson);
    if (pObsMInBoson->Integral() > 0) pObsMInBoson->Scale(1.0 / pObsMInBoson->Integral());
    pObsMInBoson->SetLineColor(kGreen + 2);
    pObsMInBoson->SetLineWidth(2);

    cBosonObsMSlices->cd(i + 1);
    pObsMInBoson->SetTitle(Form("%.1f < %s < %.1f%s;%s;Normalized entries", bosonLo, bosonLabel, bosonHi, bosonUnit, obsMinusLabel));
    pObsMInBoson->Draw("hist");

    TH1D* pObsPInBoson = hBosonVsObsP->ProjectionY(Form("pObsP_bosonbin_%d", i), xBinLoBoson, xBinHiBoson);
    if (pObsPInBoson->Integral() > 0) pObsPInBoson->Scale(1.0 / pObsPInBoson->Integral());
    pObsPInBoson->SetLineColor(kMagenta + 1);
    pObsPInBoson->SetLineWidth(2);

    cBosonObsPSlices->cd(i + 1);
    pObsPInBoson->SetTitle(Form("%.1f < %s < %.1f%s;%s;Normalized entries", bosonLo, bosonLabel, bosonHi, bosonUnit, obsPlusLabel));
    pObsPInBoson->Draw("hist");

    // In deltaEta_vs_{zPt|zY}, boson variable is on Y axis.
    int yBinLoBoson = hDEtaVsBoson->GetYaxis()->FindBin(bosonLo + 1e-6);
    int yBinHiBoson = (bosonHi >= hDEtaVsBoson->GetYaxis()->GetXmax())
                        ? hDEtaVsBoson->GetYaxis()->GetNbins()
                        : hDEtaVsBoson->GetYaxis()->FindBin(bosonHi - 1e-6);
    yBinLoBoson = std::max(1, yBinLoBoson);
    yBinHiBoson = std::max(yBinLoBoson, std::min(hDEtaVsBoson->GetYaxis()->GetNbins(), yBinHiBoson));

    TH1D* pDEtaInBoson = hDEtaVsBoson->ProjectionX(Form("pDeltaEta_bosonbin_%d", i), yBinLoBoson, yBinHiBoson);
    if (pDEtaInBoson->Integral() > 0) pDEtaInBoson->Scale(1.0 / pDEtaInBoson->Integral());
    pDEtaInBoson->SetLineColor(sliceColors[i]);
    pDEtaInBoson->SetLineWidth(2);
    pDEtaBosonSlices[i] = pDEtaInBoson;

    cDEtaInBosonSlices->cd(i + 1);
    pDEtaInBoson->SetTitle(Form("%.1f < %s < %.1f%s;|#Delta#eta(l_{1},l_{2})|;Normalized entries", bosonLo, bosonLabel, bosonHi, bosonUnit));
    pDEtaInBoson->Draw("hist");
  }

  // Overlay all deltaEta boson slices on one canvas with distinct colors.
  TCanvas* cDEtaBosonOverlay = new TCanvas(
    isY ? "c_overlay_deltaeta_in_zY_slices" : "c_overlay_deltaeta_in_zpt_slices",
    isY ? "Overlay |#Delta#eta| in |y(Z)| slices" : "Overlay |#Delta#eta| in p_{T}(Z) slices",
    1000, 800
  );

  TPad* padDEtaOverlayTop = new TPad("padDEtaOverlayTop", "padDEtaOverlayTop", 0.0, 0.30, 1.0, 1.0);
  TPad* padDEtaOverlayBot = new TPad("padDEtaOverlayBot", "padDEtaOverlayBot", 0.0, 0.00, 1.0, 0.30);
  padDEtaOverlayTop->SetBottomMargin(0.02);
  padDEtaOverlayBot->SetTopMargin(0.05);
  padDEtaOverlayBot->SetBottomMargin(0.32);
  cDEtaBosonOverlay->cd();
  padDEtaOverlayTop->Draw();
  padDEtaOverlayBot->Draw();

  double maxDeltaEtaOverlay = 0.0;
  for (int i = 0; i < nBosonSlices; ++i) {
    if (!pDEtaBosonSlices[i]) continue;
    maxDeltaEtaOverlay = std::max(maxDeltaEtaOverlay, pDEtaBosonSlices[i]->GetMaximum());
  }

  const double legX1 = isY ? 0.14 : 0.56;
  const double legY1 = isY ? 0.60 : 0.58;
  const double legX2 = isY ? 0.46 : 0.88;
  const double legY2 = 0.88;
  TLegend* legDEtaBosonOverlay = new TLegend(legX1, legY1, legX2, legY2);
  legDEtaBosonOverlay->SetBorderSize(0);
  legDEtaBosonOverlay->SetFillStyle(0);

  padDEtaOverlayTop->cd();
  for (int i = 0; i < nBosonSlices; ++i) {
    if (!pDEtaBosonSlices[i]) continue;
    if (i == 0) {
      pDEtaBosonSlices[i]->SetTitle(Form("Overlay: |#Delta#eta| slices in %s;%s;Normalized entries", bosonLabel, "|#Delta#eta(l_{1},l_{2})|"));
      if (maxDeltaEtaOverlay > 0.0) pDEtaBosonSlices[i]->SetMaximum(1.25 * maxDeltaEtaOverlay);
      pDEtaBosonSlices[i]->Draw("hist");
    } else {
      pDEtaBosonSlices[i]->Draw("hist same");
    }
    const double bosonLo = isY ? yLow[i] : bosonLow[i];
    const double bosonHi = isY ? yHigh[i] : bosonHigh[i];
    legDEtaBosonOverlay->AddEntry(pDEtaBosonSlices[i], Form("%.1f < %s < %.1f%s", bosonLo, bosonLabel, bosonHi, bosonUnit), "l");
  }
  legDEtaBosonOverlay->Draw();

  // Ratio panel: each slice divided by the first slice.
  padDEtaOverlayBot->cd();
  TH1D* pDEtaRef = pDEtaBosonSlices[0];
  TH1D* pFirstRatio = nullptr;
  for (int i = 1; i < nBosonSlices; ++i) {
    if (!pDEtaBosonSlices[i] || !pDEtaRef) continue;
    TH1D* pRatio = static_cast<TH1D*>(pDEtaBosonSlices[i]->Clone(Form("ratio_deltaeta_boson_slice_%d", i)));
    pRatio->Divide(pDEtaRef);
    pRatio->SetLineColor(sliceColors[i]);
    pRatio->SetMarkerColor(sliceColors[i]);
    pRatio->SetMarkerStyle(20);
    pRatio->SetMarkerSize(0.55);
    pRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    pRatio->GetYaxis()->SetTitle("Slice / Ref");
    pRatio->GetYaxis()->SetTitleSize(0.10);
    pRatio->GetYaxis()->SetLabelSize(0.09);
    pRatio->GetYaxis()->SetTitleOffset(0.45);
    pRatio->GetYaxis()->SetNdivisions(505);
    pRatio->GetXaxis()->SetTitle("|#Delta#eta(l_{1},l_{2})|");
    pRatio->GetXaxis()->SetTitleSize(0.12);
    pRatio->GetXaxis()->SetLabelSize(0.10);
    if (!pFirstRatio) {
      pFirstRatio = pRatio;
      pFirstRatio->Draw("ep");
    } else {
      pRatio->Draw("ep same");
    }
  }

  if (pFirstRatio && pDEtaRef) {
    const double xMin = pDEtaRef->GetXaxis()->GetXmin();
    const double xMax = pDEtaRef->GetXaxis()->GetXmax();
    TLine* unityLine = new TLine(xMin, 1.0, xMax, 1.0);
    unityLine->SetLineColor(kGray + 2);
    unityLine->SetLineStyle(2);
    unityLine->SetLineWidth(2);
    unityLine->Draw("same");

    const double refLo = isY ? yLow[0] : bosonLow[0];
    const double refHi = isY ? yHigh[0] : bosonHigh[0];
    TLegend* legRatioInfo = new TLegend(0.54, 0.72, 0.88, 0.90);
    legRatioInfo->SetBorderSize(0);
    legRatioInfo->SetFillStyle(0);
    legRatioInfo->SetTextSize(0.08);
    legRatioInfo->AddEntry((TObject*)0, Form("Ref: %.1f < %s < %.1f%s", refLo, bosonLabel, refHi, bosonUnit), "");
    legRatioInfo->AddEntry((TObject*)0, "All curves: slice / ref", "");
    legRatioInfo->Draw();
  }

  // 11) Requested overlays in Y(Z) slices: pT(leading/subleading) and eta(leading/subleading).
  if (isY) {
    if (hZYVsPtLead && hZYVsPtSub) {
      cPtLeadSubleadInYSlices = new TCanvas(
        "c_overlay_pT_lead_vs_sublead_in_yZ_slices",
        "p_{T}(leading) vs p_{T}(subleading) in |y(Z)| slices",
        1600, 900
      );
      cPtLeadSubleadInYSlices->Divide(4, 2);

      for (int i = 0; i < nBosonSlices; ++i) {
        const double yLo = yLow[i];
        const double yHi = yHigh[i];

        int xBinLo = hZYVsPtLead->GetXaxis()->FindBin(yLo + 1e-6);
        int xBinHi = (yHi >= hZYVsPtLead->GetXaxis()->GetXmax())
                       ? hZYVsPtLead->GetXaxis()->GetNbins()
                       : hZYVsPtLead->GetXaxis()->FindBin(yHi - 1e-6);
        xBinLo = std::max(1, xBinLo);
        xBinHi = std::max(xBinLo, std::min(hZYVsPtLead->GetXaxis()->GetNbins(), xBinHi));

        TH1D* pPtLeadInY = hZYVsPtLead->ProjectionY(Form("pPt_lead_yZbin_%d", i), xBinLo, xBinHi);
        TH1D* pPtSubInY  = hZYVsPtSub->ProjectionY(Form("pPt_sublead_yZbin_%d", i), xBinLo, xBinHi);

        if (pPtLeadInY->Integral() > 0) pPtLeadInY->Scale(1.0 / pPtLeadInY->Integral());
        if (pPtSubInY->Integral() > 0)  pPtSubInY->Scale(1.0 / pPtSubInY->Integral());

        pPtLeadInY->SetLineColor(kRed + 1);
        pPtLeadInY->SetLineWidth(2);
        pPtSubInY->SetLineColor(kBlue + 1);
        pPtSubInY->SetLineWidth(2);
        pPtSubInY->SetLineStyle(2);

        cPtLeadSubleadInYSlices->cd(i + 1);
        pPtLeadInY->SetTitle(Form("%.1f < |y(Z)| < %.1f;p_{T} [GeV];Normalized entries", yLo, yHi));
        pPtLeadInY->Draw("hist");
        pPtSubInY->Draw("hist same");

        TLegend* legPtInY = new TLegend(0.48, 0.74, 0.88, 0.88);
        legPtInY->AddEntry(pPtLeadInY, "p_{T}(leading)", "l");
        legPtInY->AddEntry(pPtSubInY, "p_{T}(subleading)", "l");
        legPtInY->Draw();

        // Annotate mean pT values in each |y(Z)| slice panel.
        TLatex* meanText = new TLatex();
        meanText->SetNDC();
        meanText->SetTextSize(0.035);
        meanText->SetTextFont(42);
        const TString leadMeanLabel = (pPtLeadInY->GetEntries() > 0)
                                      ? Form("<#it{p}_{T}> lead = %.1f GeV", pPtLeadInY->GetMean())
                                      : "<#it{p}_{T}> lead = n/a";
        const TString subMeanLabel = (pPtSubInY->GetEntries() > 0)
                                     ? Form("<#it{p}_{T}> sublead = %.1f GeV", pPtSubInY->GetMean())
                                     : "<#it{p}_{T}> sublead = n/a";
        meanText->DrawLatex(0.48, 0.66, leadMeanLabel);
        meanText->DrawLatex(0.48, 0.60, subMeanLabel);
      }
    } else {
      std::cout << "WARNING: missing zY_vs_pt_leading and/or zY_vs_pt_subleading. Skipping pT(leading/subleading) in Y(Z) slices." << std::endl;
    }

    if (hZYVsEtaLead && hZYVsEtaSub) {
      cEtaLeadSubleadInYSlices = new TCanvas(
        "c_overlay_eta_lead_vs_sublead_in_yZ_slices",
        "#eta(leading) vs #eta(subleading) in |y(Z)| slices",
        1600, 900
      );
      cEtaLeadSubleadInYSlices->Divide(4, 2);

      for (int i = 0; i < nBosonSlices; ++i) {
        const double yLo = yLow[i];
        const double yHi = yHigh[i];

        int xBinLo = hZYVsEtaLead->GetXaxis()->FindBin(yLo + 1e-6);
        int xBinHi = (yHi >= hZYVsEtaLead->GetXaxis()->GetXmax())
                       ? hZYVsEtaLead->GetXaxis()->GetNbins()
                       : hZYVsEtaLead->GetXaxis()->FindBin(yHi - 1e-6);
        xBinLo = std::max(1, xBinLo);
        xBinHi = std::max(xBinLo, std::min(hZYVsEtaLead->GetXaxis()->GetNbins(), xBinHi));

        TH1D* pEtaLeadInY = hZYVsEtaLead->ProjectionY(Form("pEta_lead_yZbin_%d", i), xBinLo, xBinHi);
        TH1D* pEtaSubInY  = hZYVsEtaSub->ProjectionY(Form("pEta_sublead_yZbin_%d", i), xBinLo, xBinHi);

        if (pEtaLeadInY->Integral() > 0) pEtaLeadInY->Scale(1.0 / pEtaLeadInY->Integral());
        if (pEtaSubInY->Integral() > 0)  pEtaSubInY->Scale(1.0 / pEtaSubInY->Integral());

        pEtaLeadInY->SetLineColor(kRed + 1);
        pEtaLeadInY->SetLineWidth(2);
        pEtaSubInY->SetLineColor(kBlue + 1);
        pEtaSubInY->SetLineWidth(2);
        pEtaSubInY->SetLineStyle(2);

        cEtaLeadSubleadInYSlices->cd(i + 1);
        pEtaLeadInY->SetTitle(Form("%.1f < |y(Z)| < %.1f;#eta;Normalized entries", yLo, yHi));
        pEtaLeadInY->Draw("hist");
        pEtaSubInY->Draw("hist same");

        TLegend* legEtaInY = new TLegend(0.48, 0.74, 0.88, 0.88);
        legEtaInY->AddEntry(pEtaLeadInY, "#eta(leading)", "l");
        legEtaInY->AddEntry(pEtaSubInY, "#eta(subleading)", "l");
        legEtaInY->Draw();
      }
    } else {
      std::cout << "WARNING: missing zY_vs_eta_leading and/or zY_vs_eta_subleading. Skipping eta(leading/subleading) in Y(Z) slices." << std::endl;
    }
  }

  // Save comparison plots next to the input file for quick checks.
  TString outPrefix;
  if (isFiducial) outPrefix += "fiducial_";
  if (CFonlyEta) outPrefix += "cfonlyeta_";
  if (EtaOnly) outPrefix += "etaonly_";
  if (CFonly) outPrefix += "cfonly_";
  if (CCCF) outPrefix += "cccf_";
  auto prefixed = [&](const char* name) { return outPrefix + name; };

  if (cEtaLeadSublead2D) cEtaLeadSublead2D->SaveAs(prefixed("compare_2D_eta_eleading_vs_esubleading.pdf").Data());
  cPtLeadSublead2D->SaveAs(prefixed("compare_2D_pt_leading_vs_subleading.pdf").Data());
  cCosth->SaveAs(prefixed("compare_projection_costh_lead_vs_sublead.pdf").Data());
  cDEtaPt->SaveAs(prefixed("compare_projection_deltaeta_lead_vs_sublead.pdf").Data());
  cDEtaCosth->SaveAs(prefixed("compare_projection_deltaeta_vs_costh.pdf").Data());
  cCosthPtSlices->SaveAs(prefixed("compare_projection_costh_pt_slices_lead_vs_sublead.pdf").Data());
  cDEtaPtSlices->SaveAs(prefixed("compare_projection_deltaeta_pt_slices_lead_vs_sublead.pdf").Data());
  cCosthInDEtaSlices->SaveAs(prefixed("compare_projectionY_costh_in_deltaeta_slices.pdf").Data());
  cDEtaInCosthSlices->SaveAs(prefixed("compare_projectionX_deltaeta_in_costh_slices.pdf").Data());
  cPtEpSlices->SaveAs(prefixed("compare_projection_pT_em_in_pT_ep_slices.pdf").Data());
  cPtEmSlices->SaveAs(prefixed("compare_projection_pT_ep_in_pT_em_slices.pdf").Data());
  cPtLeadSlices->SaveAs(prefixed("compare_projection_pT_sublead_in_pT_lead_slices.pdf").Data());
  cPtSubleadSlices->SaveAs(prefixed("compare_projection_pT_lead_in_pT_sublead_slices.pdf").Data());
  cPtEpEmOverlay->SaveAs(prefixed("compare_projection_overlay_pT_ep_vs_em_in_pT_slices.pdf").Data());
  cPtLeadSubleadOverlay->SaveAs(prefixed("compare_projection_overlay_pT_lead_vs_sublead_in_pT_slices.pdf").Data());
  cBosonCosthSlices->SaveAs(prefixed(isY ? "compare_projection_costheta_in_yZ_slices.pdf" : "compare_projection_costheta_in_pTZ_slices.pdf").Data());
  cBosonPhiSlices->SaveAs(prefixed(isY ? "compare_projection_phi_in_yZ_slices.pdf" : "compare_projection_phi_in_pTZ_slices.pdf").Data());
  cBosonObsMSlices->SaveAs(prefixed(isY ? "compare_projection_pT_m_in_yZ_slices.pdf" : "compare_projection_eta_m_in_pTZ_slices.pdf").Data());
  cBosonObsPSlices->SaveAs(prefixed(isY ? "compare_projection_pT_p_in_yZ_slices.pdf" : "compare_projection_eta_p_in_pTZ_slices.pdf").Data());
  cDEtaInBosonSlices->SaveAs(prefixed(isY ? "compare_projection_deltaeta_in_yZ_slices.pdf" : "compare_projection_deltaeta_in_pTZ_slices.pdf").Data());
  cDEtaBosonOverlay->SaveAs(prefixed(isY ? "compare_projection_overlay_deltaeta_in_yZ_slices.pdf" : "compare_projection_overlay_deltaeta_in_pTZ_slices.pdf").Data());
  if (cPtLeadSubleadInYSlices) cPtLeadSubleadInYSlices->SaveAs(prefixed("compare_projection_overlay_pT_lead_vs_sublead_in_yZ_slices.pdf").Data());
  if (cEtaLeadSubleadInYSlices) cEtaLeadSubleadInYSlices->SaveAs(prefixed("compare_projection_overlay_eta_lead_vs_sublead_in_yZ_slices.pdf").Data());

  std::cout << "Saved plots:" << std::endl;
  if (cEtaLeadSublead2D) std::cout << "  " << prefixed("compare_2D_eta_eleading_vs_esubleading.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_2D_pt_leading_vs_subleading.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_costh_lead_vs_sublead.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_deltaeta_lead_vs_sublead.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_deltaeta_vs_costh.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_costh_pt_slices_lead_vs_sublead.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_deltaeta_pt_slices_lead_vs_sublead.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projectionY_costh_in_deltaeta_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projectionX_deltaeta_in_costh_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_pT_em_in_pT_ep_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_pT_ep_in_pT_em_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_pT_sublead_in_pT_lead_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_pT_lead_in_pT_sublead_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_overlay_pT_ep_vs_em_in_pT_slices.pdf") << std::endl;
  std::cout << "  " << prefixed("compare_projection_overlay_pT_lead_vs_sublead_in_pT_slices.pdf") << std::endl;
  std::cout << "  " << prefixed(isY ? "compare_projection_costheta_in_yZ_slices.pdf" : "compare_projection_costheta_in_pTZ_slices.pdf") << std::endl;
  std::cout << "  " << prefixed(isY ? "compare_projection_phi_in_yZ_slices.pdf" : "compare_projection_phi_in_pTZ_slices.pdf") << std::endl;
  std::cout << "  " << prefixed(isY ? "compare_projection_pT_m_in_yZ_slices.pdf" : "compare_projection_eta_m_in_pTZ_slices.pdf") << std::endl;
  std::cout << "  " << prefixed(isY ? "compare_projection_pT_p_in_yZ_slices.pdf" : "compare_projection_eta_p_in_pTZ_slices.pdf") << std::endl;
  std::cout << "  " << prefixed(isY ? "compare_projection_deltaeta_in_yZ_slices.pdf" : "compare_projection_deltaeta_in_pTZ_slices.pdf") << std::endl;
  std::cout << "  " << prefixed(isY ? "compare_projection_overlay_deltaeta_in_yZ_slices.pdf" : "compare_projection_overlay_deltaeta_in_pTZ_slices.pdf") << std::endl;
  if (cPtLeadSubleadInYSlices) std::cout << "  " << prefixed("compare_projection_overlay_pT_lead_vs_sublead_in_yZ_slices.pdf") << std::endl;
  if (cEtaLeadSubleadInYSlices) std::cout << "  " << prefixed("compare_projection_overlay_eta_lead_vs_sublead_in_yZ_slices.pdf") << std::endl;

  //f->Close();
}
