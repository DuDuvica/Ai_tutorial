#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"
#include <math.h>  
// Version 2025 for WAI paper latest plots 
using namespace std;

bool sherpa = true;
bool test =false;
bool testSherpaNEW = true; // Sherpa221 NEW
bool overight = false;
bool normXS = true;
bool newPwg = false;
bool signSwop = false;

// Macro to plot Ai coefficient from Sherpa and Phoweg samples
void AIW(bool isPlus=false, bool isY=false){
  
  double xsecAMI = 0;
  // TString minitreewm("/eos/atlas/atlascerngroupdisk/phys-sm/LowMu2017WZ/histograms/v20200515_laperiob_prod_MicroTree/MicroTree_wminusenuanalysis_MC_13TeV_ai/GenSherpat/mc16_13TeV.564170.Sherpa_221_NNPDF30NNLO_Wminenu.e5340_s3126_r10244_r10210_p3665.root");
  // TString minitreewp("/eos/atlas/atlascerngroupdisk/phys-sm/LowMu2017WZ/histograms/v20200515_laperiob_prod_MicroTree/MicroTree_wplusenuanalysis_MC_13TeV_ai/GenSherpat/mc16_13TeV.464170.Sherpa_221_NNPDF30NNLO_Wplusenu.e5340_s3126_r10244_r10210_p3665.root");
  
 TString minitreewm("mc16_13TeV.564170.Sherpa_221_NNPDF30NNLO_Wminenu.e5340_s3126_r10244_r10210_p3665.root");
 TString minitreewp("mc16_13TeV.464170.Sherpa_221_NNPDF30NNLO_Wplusenu.e5340_s3126_r10244_r10210_p3665.root");

  if (testSherpaNEW) {
    minitreewm = "/data/dust/group/atlas/wai/MicroTrees/histograms/v20250617_dponomar_prod_MicroTree_sherpa_noweight/MicroTree_ai_current_wminusmunuanalysis_MC_13TeV_ai/GenSyst/mc16_13TeV.564156.Sherpa_221_NNPDF30NNLO_Wminmunu.e5340_s3126_r10244_r10210_p3654.root";
    minitreewp = "/data/dust/group/atlas/wai/MicroTrees/histograms/v20250617_dponomar_prod_MicroTree_sherpa_noweight/MicroTree_ai_current_wplusmunuanalysis_MC_13TeV_ai/GenSyst/mc16_13TeV.464156.Sherpa_221_NNPDF30NNLO_Wplusmunu.e5340_s3126_r10244_r10210_p3654.root"; 
  }

 if (sherpa) {
   newPwg = false;
   if (!isPlus) xsecAMI = 8.662; // nb Sherpa221 W- 
   else xsecAMI = 11.259; // nb Sherpa221 W+
   if (testSherpaNEW) signSwop = false; // Sherpa221 NEW is not sign swop
  }
 if (!sherpa) {
   testSherpaNEW = false;
   signSwop = false;
   minitreewm = "/data/dust/user/ruth/af-atlas.merged/Ruth/Wai/HistMaker_Output/histograms/v20200611_ruth_prod_MicroTree/MicroTree_wminusenuanalysis_MC_13TeV_ai/Nominal/mc16_13TeV.361103.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusenu.e3601_s3126_r10244+r11165_p3665.root";
   // "/nfs/dust/atlas/user/ruth/Ruth/Wai/HistMaker_Output/histograms/v20200611_ruth_prod_MicroTree/MicroTree_wminusenuanalysis_MC_13TeV_ai/Nominal/mc16_13TeV.361103.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusenu.e3601_s3126_r10244+r11165_p3665.root";
   if (newPwg) minitreewm = "/data/dust/group/atlas/wai/MicroTrees/histograms/v20240815_dponomar_prod_MicroTree_nom/MicroTree_ai_current_wminusenuanalysis_MC_13TeV_ai/Nominal/mc16_13TeV.361103.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wminusenu.e3601_s3126_r10244+r11165_p3665.root";
   
   minitreewp = "/data/dust/user/ruth/af-atlas.merged/Ruth/Wai/HistMaker_Output/histograms/v20200611_ruth_prod_MicroTree/MicroTree_wplusenuanalysis_MC_13TeV_ai/Nominal/mc16_13TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e3601_s3126_r10244+r11165_p3665.root";
   //"/nfs/dust/atlas/user/ruth/Ruth/Wai/HistMaker_Output/histograms/v20200611_ruth_prod_MicroTree/MicroTree_wplusenuanalysis_MC_13TeV_ai/Nominal/mc16_13TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e3601_s3126_r10244+r11165_p3665.root";
   if (newPwg) minitreewp = "/data/dust/group/atlas/wai/MicroTrees/histograms/v20240815_dponomar_prod_MicroTree_nom/MicroTree_ai_current_wplusenuanalysis_MC_13TeV_ai/Nominal/mc16_13TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e3601_s3126_r10244+r11165_p3665.root";
   
   if (!isPlus) xsecAMI = 8.28; // nb Phoweg W-
   else xsecAMI = 11.3; // nb Phoweg W+
 }
 
 TString outputName = "Wai_finalbinning_newSign";        
 if (test) outputName = "test";                                                                                                                             
 if (!sherpa) {
   outputName = outputName+"Phoweg";  
   if (newPwg) outputName = outputName+"NEW";
 } else { 
   outputName = outputName+"Sherpa221"; // Sherpa221 is the default for Sherpa
   if (testSherpaNEW) outputName = outputName+"NEW";
   if (isPlus && signSwop) outputName = outputName+"swapSign";
 }
 
 TString minitree = minitreewm;
 TString charge = "Wm";
 TString mode ="_pT";
 
 if (isY) mode = "_Y";
 if (normXS) mode = mode + "_NormXsec";  
 
 if (isPlus) {
   minitree = minitreewp;
   charge = "Wp";
 }
 std::cout << "Looking at File " << minitree << std::endl; 
 
 
 TFile f1(minitree);
 //  TFile f2("AI_W"+charge+"_"+outputName+mode+".root","RECREATE");
 // saving the output
 TString nameOutput = "AI_W"+charge+"_"+outputName+mode+".root" ;
 TFile* Output = new TFile(nameOutput);
 bool isf = true ;
 
 if (Output->IsZombie()) {
   std::cout << " OUTPUT FILE DO NOT EXIST" << std::endl;
   isf = false;
 }
 
 if (isf) {
   if (overight)  {
     Output  = new TFile(nameOutput,"RECREATE");
   } else {
     cout << "--> File " << nameOutput << " exist are you sure you want to overight ? if yes put overight to true " << endl;  
     return ;
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
 if (sherpa) denom = hNorm->GetBinContent(1); // Sherpa221 sum of weights is in bin 1
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
 
 
 TTree *tree = (TTree*)f1.Get("MicroTree/HWWTree_ee"); 
 if (!tree) {
   std::cout << " HWWTree_ee DOES NOT EXIST" << std::endl;
   return;
 }
 const int Nbins = 11.;//12; //26;
 //  Double_t bins[Nbins]= {0,2.5,5.0,8.0,11.4,14.9,18.5,22.0,25.5,29.0,32.6,36.4,40.4,44.9,50.2,56.4,63.9,73.4,85.4,105.0,132.0,173.0,253.0,600.0,900.0,1000.0};
 // Double_t bins[Nbins]= { 0. , 8. , 17., 27., 40., 55., 75., 95., 120., 160., 220., 600};
 Double_t bins[Nbins] = { 0., 8., 17., 27., 40., 55., 75., 110., 150., 210., 600.}; // final WAi analyisi
 if (isY) {
   double width = 0.4;
   double high =5.0;
   for (int k =0 ; k<Nbins ;k++) {
     //      bins[k] = (high/(Nbins-1))*k;
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
 
 
 Double_t pt[100], eta[100], phi[100];
 Double_t  costheta;
 Long64_t N = tree->GetEntries();
 if (test)  N = 100000;
 int nchild;
 Double_t lpe, lme, lpnu, lmnu, tanphi;
 
 Double_t        lepPtTruth0;
 Double_t        lepPhiTruth0;
 Double_t        lepEtaTruth0;
 Double_t        lepPtTruth1;
 Double_t        lepPhiTruth1;
 Double_t        lepEtaTruth1;
 Double_t        genWeight;
 tree->SetBranchAddress("genWeight", &genWeight);
 
 tree->SetBranchAddress("lepPtTruth0", &lepPtTruth0);
 tree->SetBranchAddress("lepPhiTruth0", &lepPhiTruth0);
 tree->SetBranchAddress("lepEtaTruth0", &lepEtaTruth0);
 tree->SetBranchAddress("lepPtTruth1", &lepPtTruth1);
 tree->SetBranchAddress("lepPhiTruth1", &lepPhiTruth1);
 tree->SetBranchAddress("lepEtaTruth1", &lepEtaTruth1);
 
 // already selected :
 //Particle 1
 //Antiparticle 0
 
 TLorentzVector e, nu, w, delta;
 TVector3 pe, pnu, pw, rt, pa(0,0,1),deltat, qttrans;
 
 for (Long64_t i = 0; i < N; i++) {
   
   tree->GetEntry(i);
   e.SetPtEtaPhiM(lepPtTruth1/1000.0,lepEtaTruth1,lepPhiTruth1,0);
   nu.SetPtEtaPhiM(lepPtTruth0/1000.0,lepEtaTruth0,lepPhiTruth0,0);
   w = e + nu;
   delta = nu - e;
   deltat = delta.Vect();
   deltat.SetZ(0);
   if (i%1000000 == 0) cout << "Processed " << i <<"/" << N << " events" << endl;
   if (i%1000000 == 0) cout << " GenWeight:   " << genWeight << "  NORM:  " << norm << endl;
   lpe = e.E()+e.Pz();
   lme = e.E()-e.Pz();
   lpnu =  nu.E()+nu.Pz();
   lmnu =  nu.E()-nu.Pz();
   costheta = (1./w.M())*pow((w.M2()+pow(w.Pt(),2)),-0.5)*(lpe*lmnu-lme*lpnu); 
   costheta*=w.Rapidity()>=0. ? 1:-1;
   if (costheta!=costheta) continue;
   //cout<< "costheta "<<costheta<<endl;
   if (isPlus && signSwop) costheta*=-1;
   double weight = genWeight;
   if (normXS) weight *= norm;
   hctheta->Fill(costheta, weight);
   
   pe.SetPtEtaPhi(lepPtTruth1/1000.0,lepEtaTruth1,lepPhiTruth1);
   pnu.SetPtEtaPhi(lepPtTruth0/1000.0,lepEtaTruth0,lepPhiTruth0);
   pw = pe + pnu;
   rt =pa.Cross(pw);
   rt = rt*(1./rt.Mag());
   rt.SetZ(0);
   qttrans = w.Vect();
   qttrans.SetZ(0);
   qttrans = qttrans*(1./qttrans.Mag());
   tanphi = (pow(w.M2()+pow(w.Pt(),2),0.5)/w.M())*(deltat*rt/(deltat*qttrans));
   Double_t x = pow((pow(w.M2()+pow(w.Pt(),2),0.5)/w.M()),-1)*(deltat*qttrans);
   Double_t y = deltat*rt*(w.Rapidity()>=0. ? 1:-1);
   if (isPlus && signSwop) y*=-1;
   Double_t phi = atan2(y,x);
   hphi->Fill(phi, weight);
   
   if(isY){
     if (normXS) {
       Xs->Fill(fabs(w.Rapidity()));
       Xsw->Fill(fabs(w.Rapidity()),genWeight*norm);
     } else {
       Xs->Fill(fabs(w.Rapidity()));
       Xsw->Fill(fabs(w.Rapidity()),genWeight);
     }
   } else {
     if (normXS) {
       Xs->Fill(w.Pt());
       Xsw->Fill(w.Pt(),genWeight*norm);
     } else {
       Xs->Fill(w.Pt());
       Xsw->Fill(w.Pt(),genWeight);
     }
   }
   Double_t sinth     = -sqrt(max(0.,1.-costheta*costheta));
   Double_t sin2th     = 2*costheta*sinth;
   Double_t cosph   = cos(phi);
   Double_t cos2ph  = 2*cos(phi)*cos(phi)-1.;
   Double_t sinph   = sqrt(max(0.,1.-cos(phi)*cos(phi)))*(phi>0?1:-1);
   Double_t sin2ph  = 2*cos(phi)*sin(phi);
   
   if(isY){
     A0->Fill(fabs(w.Rapidity()), (20./3.*(0.5-1.5*costheta*costheta)+2./3.) * weight);
     A1->Fill(fabs(w.Rapidity()), (5.*(2.*costheta*sinth*cos(phi))) * weight);
     A2->Fill(fabs(w.Rapidity()), (10.*(sinth*sinth*cos2ph)) * weight);
     A3->Fill(fabs(w.Rapidity()), (4.*(sinth*cos(phi))) * weight);
     A4->Fill(fabs(w.Rapidity()), (4.*costheta) * weight);
     A5->Fill(fabs(w.Rapidity()), (5. * (sinth*sinth*sin2ph)) * weight);
     A6->Fill(fabs(w.Rapidity()), (4. * (sin2th*sinph)) * weight);
     A7->Fill(fabs(w.Rapidity()), (4. * (sinth*sinph)) * weight);
   }else {
     A0->Fill(w.Pt(), (20./3.*(0.5-1.5*costheta*costheta)+2./3.) * weight);
     A1->Fill(w.Pt(), (5.*(2.*costheta*sinth*cos(phi))) * weight);
     A2->Fill(w.Pt(), (10.*(sinth*sinth*cos2ph)) * weight);
     A3->Fill(w.Pt(), (4.*(sinth*cos(phi))) * weight);
     A4->Fill(w.Pt(), (4.*costheta) * weight);
     A5->Fill(w.Pt(), (5. * (sinth*sinth*sin2ph)) * weight);
     A6->Fill(w.Pt(), (4. * (sin2th*sinph)) * weight);
     A7->Fill(w.Pt(), (4. * (sinth*sinph)) * weight);
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
