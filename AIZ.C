#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <math.h>
// Draft Z version based on AIW.C (Collins-Soper angles, A0-A7)
using namespace std;

bool sherpa = true;
bool test = false;
bool overight = false;
bool normXS = true;

// Macro to plot Ai coefficient from Sherpa and Powheg Z samples
void AIZ(bool isY=false){

  double xsecAMI = 0;

  // TODO: update to your Z samples (Sherpa vs Powheg) and tree path
  TString minitree = "mc16_13TeV.Zee.Sherpa_or_Powheg.root";

  if (sherpa) {
    xsecAMI = 1.0; // TODO: set Z cross-section (nb)
  } else {
    xsecAMI = 1.0; // TODO: set Z cross-section (nb)
  }

  TString outputName = "Zai_finalbinning";
  if (test) outputName = "test";
  if (sherpa) outputName = outputName+"Sherpa";
  else outputName = outputName+"Powheg";

  TString mode = "_pT";
  if (isY) mode = "_Y";
  if (normXS) mode = mode + "_NormXsec";

  std::cout << "Looking at File " << minitree << std::endl;

  TFile f1(minitree);
  TString nameOutput = "AI_Z_"+outputName+mode+".root";
  TFile* Output = new TFile(nameOutput);
  bool isf = true;

  if (Output->IsZombie()) {
    std::cout << " OUTPUT FILE DO NOT EXIST" << std::endl;
    isf = false;
  }

  if (isf) {
    if (overight)  {
      Output  = new TFile(nameOutput,"RECREATE");
    } else {
      cout << "--> File " << nameOutput << " exist are you sure you want to overight ? if yes put overight to true " << endl;
      return;
    }
  } else {
    cout << " Create --> File " << nameOutput  << endl;
    Output  = new TFile(nameOutput,"RECREATE");
  }

  TH1D *hNorm = (TH1D*) f1.Get("MicroTree/Ai_Weights");
  if (!hNorm or hNorm->IsZombie()) {
    std::cout << " hNorm  DO NOT EXIST" << std::endl;
    return ;
  }

  // double norm=hNorm->GetBinContent(1)/hNorm->GetEntries();
  double denom = hNorm->GetBinContent(5);
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

  // TODO: update tree path for Z microtree if different
  TTree *tree = (TTree*)f1.Get("MicroTree/HWWTree_ee");
  if (!tree) {
    std::cout << " HWWTree_ee DOES NOT EXIST" << std::endl;
    return;
  }

  const int Nbins = 11.;
  Double_t bins[Nbins] = { 0., 8., 17., 27., 40., 55., 75., 110., 150., 210., 600.};
  if (isY) {
    double width = 0.4;
    for (int k =0 ; k<Nbins ;k++) {
      bins[k] = (width)*k;
      cout << " Y  bins at:  " << k << " bin " << bins[k] << endl;
    }
  }

  TH1D *hctheta = new TH1D("Costheta", "", 100, -1., 1.0);
  TH1D *hphi = new TH1D("phi", "", 100, -M_PI, M_PI);
  TH1D *Xs = new TH1D("Xs", "Xs", Nbins-1, bins);
  TH1D *A0 = new TH1D("A0", "", Nbins-1, bins);
  TH1D *A1 = new TH1D("A1", "", Nbins-1, bins);
  TH1D *A2 = new TH1D("A2", "", Nbins-1, bins);
  TH1D *A3 = new TH1D("A3", "", Nbins-1, bins);
  TH1D *A4 = new TH1D("A4", "", Nbins-1, bins);
  TH1D *A5 = new TH1D("A5", "", Nbins-1, bins);
  TH1D *A6 = new TH1D("A6", "", Nbins-1, bins);
  TH1D *A7 = new TH1D("A7", "", Nbins-1, bins);
  TH1D *Xsw = new TH1D("Xsweighted", "Xsweighted", Nbins-1, bins);

  Long64_t N = tree->GetEntries();
  if (test)  N = 100000;

  Double_t lepPtTruth0, lepPhiTruth0, lepEtaTruth0;
  Double_t lepPtTruth1, lepPhiTruth1, lepEtaTruth1;
  Double_t genWeight;

  tree->SetBranchAddress("genWeight", &genWeight);
  tree->SetBranchAddress("lepPtTruth0", &lepPtTruth0);
  tree->SetBranchAddress("lepPhiTruth0", &lepPhiTruth0);
  tree->SetBranchAddress("lepEtaTruth0", &lepEtaTruth0);
  tree->SetBranchAddress("lepPtTruth1", &lepPtTruth1);
  tree->SetBranchAddress("lepPhiTruth1", &lepPhiTruth1);
  tree->SetBranchAddress("lepEtaTruth1", &lepEtaTruth1);

  TLorentzVector em, ep, z, delta;
  TVector3 pem, pep, pz, rt, pa(0,0,1), deltat, qttrans;

  for (Long64_t i = 0; i < N; i++) {

    tree->GetEntry(i);
    // Convention: angles are defined w.r.t. the negative lepton.
    // z-axis: direction of positive longitudinal Z momentum in lab (sign of Z rapidity).
    // y-axis: normal to the plane spanned by incoming proton momenta (beam axis) and Z momentum.
    // x-axis: completes a right-handed system.
    // Note: when pT(Z) == 0, phi is undefined; we keep events (phi will be arbitrary).
    // Convention: lepTruth1 is negative lepton, lepTruth0 is positive lepton <-- check this in your samples!
    em.SetPtEtaPhiM(lepPtTruth1/1000.0, lepEtaTruth1, lepPhiTruth1, 0);
    ep.SetPtEtaPhiM(lepPtTruth0/1000.0, lepEtaTruth0, lepPhiTruth0, 0);
    z = em + ep;

    delta = em - ep;
    deltat = delta.Vect();
    deltat.SetZ(0);

    if (i%1000000 == 0) cout << "Processed " << i <<"/" << N << " events" << endl;

    double lpe  = em.E()+em.Pz();
    double lme  = em.E()-em.Pz();
    double lpeplus = ep.E()+ep.Pz();
    double lmeplus = ep.E()-ep.Pz();

    double costheta = (1./z.M())*pow((z.M2()+pow(z.Pt(),2)),-0.5)*(lpeplus*lme - lmeplus*lpe);
    costheta *= z.Rapidity()>=0. ? 1:-1;
    if (costheta!=costheta) continue;

    double weight = genWeight;
    if (normXS) weight *= norm;

    hctheta->Fill(costheta, weight);

    pem.SetPtEtaPhi(lepPtTruth1/1000.0, lepEtaTruth1, lepPhiTruth1);
    pep.SetPtEtaPhi(lepPtTruth0/1000.0, lepEtaTruth0, lepPhiTruth0);
    pz = pem + pep;

    rt = pa.Cross(pz);
    if (rt.Mag() != 0) rt = rt*(1./rt.Mag());
    rt.SetZ(0);

    qttrans = z.Vect();
    qttrans.SetZ(0);
    if (qttrans.Mag() != 0) qttrans = qttrans*(1./qttrans.Mag());

    double x = pow((pow(z.M2()+pow(z.Pt(),2),0.5)/z.M()),-1)*(deltat*qttrans);
    double y = deltat*rt*(z.Rapidity()>=0. ? 1:-1);
    double phi = atan2(y,x);

    hphi->Fill(phi, weight);

    if(isY){
      Xs->Fill(fabs(z.Rapidity()));
      Xsw->Fill(fabs(z.Rapidity()), weight);
    } else {
      Xs->Fill(z.Pt());
      Xsw->Fill(z.Pt(), weight);
    }

    double sinth  = -sqrt(max(0.,1.-costheta*costheta));
    double sin2th = 2*costheta*sinth;
    double cosph  = cos(phi);
    double cos2ph = 2*cos(phi)*cos(phi)-1.;
    double sinph  = sqrt(max(0.,1.-cos(phi)*cos(phi)))*(phi>0?1:-1);
    double sin2ph = 2*cos(phi)*sin(phi);

    if(isY){
      A0->Fill(fabs(z.Rapidity()), (20./3.*(0.5-1.5*costheta*costheta)+2./3.) * weight);
      A1->Fill(fabs(z.Rapidity()), (5.*(2.*costheta*sinth*cosph)) * weight);
      A2->Fill(fabs(z.Rapidity()), (10.*(sinth*sinth*cos2ph)) * weight);
      A3->Fill(fabs(z.Rapidity()), (4.*(sinth*cosph)) * weight);
      A4->Fill(fabs(z.Rapidity()), (4.*costheta) * weight);
      A5->Fill(fabs(z.Rapidity()), (5.*(sinth*sinth*sin2ph)) * weight);
      A6->Fill(fabs(z.Rapidity()), (4.*(sin2th*sinph)) * weight);
      A7->Fill(fabs(z.Rapidity()), (4.*(sinth*sinph)) * weight);
    } else {
      A0->Fill(z.Pt(), (20./3.*(0.5-1.5*costheta*costheta)+2./3.) * weight);
      A1->Fill(z.Pt(), (5.*(2.*costheta*sinth*cosph)) * weight);
      A2->Fill(z.Pt(), (10.*(sinth*sinth*cos2ph)) * weight);
      A3->Fill(z.Pt(), (4.*(sinth*cosph)) * weight);
      A4->Fill(z.Pt(), (4.*costheta) * weight);
      A5->Fill(z.Pt(), (5.*(sinth*sinth*sin2ph)) * weight);
      A6->Fill(z.Pt(), (4.*(sin2th*sinph)) * weight);
      A7->Fill(z.Pt(), (4.*(sinth*sinph)) * weight);
    }
  }

  A0->Divide(Xsw);
  A1->Divide(Xsw);
  A2->Divide(Xsw);
  A3->Divide(Xsw);
  A4->Divide(Xsw);
  A5->Divide(Xsw);
  A6->Divide(Xsw);
  A7->Divide(Xsw);

  hctheta->Write();
  hphi->Write();
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
}
