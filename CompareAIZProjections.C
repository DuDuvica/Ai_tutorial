#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include <iostream>
#include <algorithm>

void CompareAIZProjections(const TString& inFile = "AI_Z_Truth_Zai_finalbinningPowheg_pT_NormXsec.root") {
  TFile* f = TFile::Open(inFile, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "ERROR: cannot open file " << inFile << std::endl;
    return;
  }

  TH2D* hCosthLead = static_cast<TH2D*>(f->Get("costh_vs_pt_leading"));
  TH2D* hCosthSub  = static_cast<TH2D*>(f->Get("costh_vs_pt_subleading"));
  TH2D* hDEtaLead  = static_cast<TH2D*>(f->Get("deltaEta_vs_leadingPt_ll"));
  TH2D* hDEtaSub   = static_cast<TH2D*>(f->Get("deltaEta_vs_subleadingPt_ll"));
  TH2D* hDEtaCosth = static_cast<TH2D*>(f->Get("deltaEta_vs_costh_ll"));

  if (!hCosthLead || !hCosthSub || !hDEtaLead || !hDEtaSub || !hDEtaCosth) {
    std::cout << "ERROR: one or more required histograms are missing in " << inFile << std::endl;
    std::cout << "Needed: costh_vs_pt_leading, costh_vs_pt_subleading, "
              << "deltaEta_vs_leadingPt_ll, deltaEta_vs_subleadingPt_ll, deltaEta_vs_costh_ll"
              << std::endl;
    f->Close();
    return;
  }

  gStyle->SetOptStat(0);

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

  // Save comparison plots next to the input file for quick checks.
  cCosth->SaveAs("compare_projection_costh_lead_vs_sublead.pdf");
  cDEtaPt->SaveAs("compare_projection_deltaeta_lead_vs_sublead.pdf");
  cDEtaCosth->SaveAs("compare_projection_deltaeta_vs_costh.pdf");
  cCosthPtSlices->SaveAs("compare_projection_costh_pt_slices_lead_vs_sublead.pdf");
  cDEtaPtSlices->SaveAs("compare_projection_deltaeta_pt_slices_lead_vs_sublead.pdf");
  cCosthInDEtaSlices->SaveAs("compare_projectionY_costh_in_deltaeta_slices.pdf");
  cDEtaInCosthSlices->SaveAs("compare_projectionX_deltaeta_in_costh_slices.pdf");

  std::cout << "Saved plots:" << std::endl;
  std::cout << "  compare_projection_costh_lead_vs_sublead.pdf" << std::endl;
  std::cout << "  compare_projection_deltaeta_lead_vs_sublead.pdf" << std::endl;
  std::cout << "  compare_projection_deltaeta_vs_costh.pdf" << std::endl;
  std::cout << "  compare_projection_costh_pt_slices_lead_vs_sublead.pdf" << std::endl;
  std::cout << "  compare_projection_deltaeta_pt_slices_lead_vs_sublead.pdf" << std::endl;
  std::cout << "  compare_projectionY_costh_in_deltaeta_slices.pdf" << std::endl;
  std::cout << "  compare_projectionX_deltaeta_in_costh_slices.pdf" << std::endl;

  //f->Close();
}
