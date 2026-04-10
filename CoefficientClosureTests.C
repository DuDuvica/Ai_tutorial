#include "TLVUtils.h"

#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TString.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

namespace {

enum class PhiRelation {
    ExactMatch,
    NegPhi,
    PhiPlusPi,
    NegPhiPlusPi,
    Other
};

constexpr int kNAi = 8;
const std::array<const char*, kNAi> kAiNames = {
    "A0", "A1", "A2", "A3", "A4", "A5", "A6", "A7"};
const std::array<double, kNAi> kAiMin = {
    -6.0, -5.0, -10.0, -4.0, -4.0, -5.0, -4.0, -4.0};
const std::array<double, kNAi> kAiMax = {
    4.0, 5.0, 10.0, 4.0, 4.0, 5.0, 4.0, 4.0};

const char* inputConventionLabel(int inputConvention) {
    switch (inputConvention) {
        case 0:
            return "lep=l, antilep=lbar for both W- and W+";
        case 1:
            return "W-: (lep, antilep), W+: swapped (antilep, lep)";
        case 2:
            return "swapped (antilep, lep) for both W- and W+";
        default:
            return "invalid";
    }
}

TString inputConventionSuffix(int inputConvention) {
    return Form("conv%d", inputConvention);
}

struct CSAngles {
    double costh = 0.0;
    double phi = 0.0;
};

struct EventSummary {
    TLorentzVector lepton;
    TLorentzVector antilepton;
    int wCharge = 0;
    CSAngles tlv;
    CSAngles kd;
    double absCosthDiff = 0.0;
    double absPhiDiff = 0.0;
    bool costhSignFlip = false;
    bool phiSignFlip = false;
    PhiRelation phiRelation = PhiRelation::Other;
};

std::array<double, kNAi> computeScaledAi(double costh, double phi) {
    const double sinth = std::sqrt(std::max(0.0, 1.0 - costh * costh));
    const double sin2th = 2.0 * costh * sinth;
    const double cosph = std::cos(phi);
    const double cos2ph = std::cos(2.0 * phi);
    const double sinph = std::sin(phi);
    const double sin2ph = std::sin(2.0 * phi);

    std::array<double, kNAi> a{};
    a[0] = 20. / 3. * (0.5 - 1.5 * costh * costh) + 2. / 3.;
    a[1] = 5. * (sin2th * cosph);
    a[2] = 10. * (sinth * sinth * cos2ph);
    a[3] = 4. * (sinth * cosph);
    a[4] = 4. * costh;
    a[5] = 5. * (sinth * sinth * sin2ph);
    a[6] = 4. * (sin2th * sinph);
    a[7] = 4. * (sinth * sinph);
    return a;
}

struct CategoryStats {
    int nEvents = 0;
    int failCosth = 0;
    int failPhi = 0;
    int signFlipCosth = 0;
    int signFlipPhi = 0;
    int nExactMatch = 0;
    int nNegPhi = 0;
    int nPhiPlusPi = 0;
    int nNegPhiPlusPi = 0;
    int nOther = 0;
    double maxAbsCosthDiff = 0.0;
    double maxAbsPhiDiff = 0.0;
    EventSummary worstCosth;
    EventSummary worstPhi;
    std::array<EventSummary, 5> sampleOtherEvents{};
    int nStoredOtherEvents = 0;
};

// Keep the exact scalar convention from P4Helpers::deltaPhi:
// deltaPhi(phiA, phiB) = wrapped(phiB - phiA) in [-pi, pi[.
double p4DeltaPhi(double phiA, double phiB) {
    return -std::remainder(-phiA + phiB, 2.0 * M_PI);
}

// Wrap a single angle phi to [-pi, pi[.
// Must use p4DeltaPhi(phi, 0) = -remainder(-phi, 2pi), NOT p4DeltaPhi(0, phi).
// p4DeltaPhi(0, phi) = -remainder(phi, 2pi) negates angles already in (-pi,pi),
// so for phi=0.7 it returns -0.7 and for phi=3.84 it returns +2.44 — both wrong.
double wrapPhi(double phi) { return p4DeltaPhi(phi, 0.0); }

double phiDiff(double lhs, double rhs) { return p4DeltaPhi(rhs, lhs); }

double phiClassTol(double tol) { return std::max(tol, 1e-6); }

// Pi-shift checks are more sensitive to branch-cut numerics around +/-pi.
double phiPiShiftTol(double tol) { return std::max(tol, 1e-4); }

bool isSignFlip(double lhs, double rhs, double tol) {
    return std::fabs(std::fabs(lhs) - std::fabs(rhs)) <= 1e-5;
}

bool phiMatches(double lhs, double rhs, double tol) {
    return std::fabs(phiDiff(lhs, rhs)) <= phiClassTol(tol);
}

bool phiMatchesPiShift(double lhs, double rhs, double tol) {
    const double d = std::fabs(phiDiff(lhs, rhs));
    return std::fabs(d - TMath::Pi()) <= phiPiShiftTol(tol);
}

PhiRelation classifyPhiRelation(double phiTLV, double phiKD, double tol) {
    if (phiMatches(phiTLV, phiKD, tol))
        return PhiRelation::ExactMatch;
    if (phiMatches(phiTLV, -phiKD, tol))
        return PhiRelation::NegPhi;
    if (phiMatchesPiShift(phiTLV, phiKD, tol))
        return PhiRelation::PhiPlusPi;
    if (phiMatches(phiTLV, -phiKD + TMath::Pi(), tol))
        return PhiRelation::NegPhiPlusPi;
    return PhiRelation::Other;
}

const char* phiRelationName(PhiRelation relation) {
    switch (relation) {
        case PhiRelation::ExactMatch:
            return "exact match";
        case PhiRelation::NegPhi:
            return "phi -> -phi";
        case PhiRelation::PhiPlusPi:
            return "phi -> phi + pi";
        case PhiRelation::NegPhiPlusPi:
            return "phi -> -phi + pi";
        case PhiRelation::Other:
            return "other";
    }
    return "other";
}

const char* flippedCoefficients(PhiRelation relation) {
    switch (relation) {
        case PhiRelation::ExactMatch:
            return "none";
        case PhiRelation::NegPhi:
            return "A5, A6, A7";
        case PhiRelation::PhiPlusPi:
            return "A1, A3, A6, A7";
        case PhiRelation::NegPhiPlusPi:
            return "A1, A3, A5";
        case PhiRelation::Other:
            return "event-dependent / not one of the 4 templates";
    }
    return "event-dependent / not one of the 4 templates";
}

void printPhiClassificationSelfCheck(double tol) {
    std::cout << "\n=== Phi classification self-check (P4Helpers convention) ===\n";

    const double phiKD = 0.7;
    struct PhiCase {
        const char* expectedLabel;
        double phiTLV;
        double phiKD;
    };

    const std::array<PhiCase, 4> tests = {{
        {"exact match", phiKD, phiKD},
        {"phi -> -phi", -phiKD, phiKD},
        {"phi -> phi + pi", wrapPhi(phiKD + TMath::Pi()), phiKD},
        {"phi -> -phi + pi", wrapPhi(-phiKD + TMath::Pi()), phiKD},
    }};

    std::cout << std::left << std::setw(18) << "expected"
              << std::setw(18) << "classified"
              << std::setw(16) << "phiTLV"
              << std::setw(16) << "phiKD"
              << std::setw(16) << "wrap(tlv-kd)"
              << std::setw(16) << "wrap(tlv+kd)"
              << "coef sign flips" << "\n";

    for (const auto& t : tests) {
        const PhiRelation rel = classifyPhiRelation(t.phiTLV, t.phiKD, tol);
        std::cout << std::left << std::setw(18) << t.expectedLabel
                  << std::setw(18) << phiRelationName(rel)
                  << std::setw(16) << std::scientific << t.phiTLV
                  << std::setw(16) << t.phiKD
                  << std::setw(16) << phiDiff(t.phiTLV, t.phiKD)
                  << std::setw(16) << wrapPhi(t.phiTLV + t.phiKD)
                  << flippedCoefficients(rel)
                  << std::defaultfloat << "\n";
    }
}

void updateCategoryStats(CategoryStats& stats, const EventSummary& event,
                         double tol) {
    ++stats.nEvents;

    if (event.absCosthDiff > stats.maxAbsCosthDiff) {
        stats.maxAbsCosthDiff = event.absCosthDiff;
        stats.worstCosth = event;
    }
    if (event.absPhiDiff > stats.maxAbsPhiDiff) {
        stats.maxAbsPhiDiff = event.absPhiDiff;
        stats.worstPhi = event;
    }

    if (event.absCosthDiff > tol && !event.costhSignFlip)
        ++stats.failCosth;
    if (event.absPhiDiff > tol && !event.phiSignFlip)
        ++stats.failPhi;
    if (event.costhSignFlip)
        ++stats.signFlipCosth;
    if (event.phiSignFlip)
        ++stats.signFlipPhi;

    switch (event.phiRelation) {
        case PhiRelation::ExactMatch:
            ++stats.nExactMatch;
            break;
        case PhiRelation::NegPhi:
            ++stats.nNegPhi;
            break;
        case PhiRelation::PhiPlusPi:
            ++stats.nPhiPlusPi;
            break;
        case PhiRelation::NegPhiPlusPi:
            ++stats.nNegPhiPlusPi;
            break;
        case PhiRelation::Other:
            ++stats.nOther;
            if (stats.nStoredOtherEvents <
                static_cast<int>(stats.sampleOtherEvents.size())) {
                stats.sampleOtherEvents[stats.nStoredOtherEvents] = event;
                ++stats.nStoredOtherEvents;
            }
            break;
    }
}

TLorentzVector makeMasslessLepton(TRandom3& rng) {
    const double pt = rng.Uniform(5.0, 120.0);
    const double eta = rng.Uniform(-5, 5);
    const double phi = rng.Uniform(-TMath::Pi(), TMath::Pi());

    TLorentzVector lep;
    lep.SetPtEtaPhiM(pt, eta, phi, 0.0);
    return lep;
}

bool isWellDefinedBoson(const TLorentzVector& boson) {
    const double m2 = boson.Mag2();
    const double pt2 = boson.Pt() * boson.Pt();
    const double mass = boson.Mag();
    const double rapidity = abs(boson.Rapidity());
    return std::fabs(boson.Pz()) > 1e-12 && m2 > 1e-12 && (m2 + pt2) > 1e-12 &&
           mass >= 79.0 && mass <= 81.0 &&
           rapidity >= 0.0 && rapidity <= 5.0;
}

CSAngles kinematicDefinitionsAngles(const TLorentzVector& lepton,
                                    const TLorentzVector& antilepton) {
    const TLorentzVector boson = lepton + antilepton;

    const double m2 = boson.Mag2();
    const double pt = boson.Pt();
    const double pt2 = pt * pt;

    const double vplusLm = lepton.E() + lepton.Pz();
    const double vminusLm = lepton.E() - lepton.Pz();
    const double vplusLp = antilepton.E() + antilepton.Pz();
    const double vminusLp = antilepton.E() - antilepton.Pz();

    CSAngles out;
    out.costh = (vplusLm * vminusLp - vplusLp * vminusLm) /
                std::sqrt(m2 * (m2 + pt2));
    out.costh *= boson.Pz() < 0.0 ? -1.0 : 1.0;
    out.costh = std::max(-1.0, std::min(1.0, out.costh));

    double plxCS = 0.0;
    double plyCS = 0.0;
    if (pt < 1e-12) {
        plxCS = lepton.Px();
        plyCS = lepton.Py();
    } else {
        const double m = boson.Mag();
        const double deltaX = lepton.Px() - antilepton.Px();
        const double deltaY = lepton.Py() - antilepton.Py();
        const double qhtX = boson.Px() / pt;
        const double qhtY = boson.Py() / pt;
        const double rhtX = -boson.Py() / pt;
        const double rhtY = boson.Px() / pt;

        plxCS = 0.5 * m / std::sqrt(m2 + pt2) * (deltaX * qhtX + deltaY * qhtY);
        plyCS = 0.5 * (deltaX * rhtX + deltaY * rhtY);
    }

    const double sign = boson.Pz() > 0.0 ? 1.0 : -1.0;
    out.phi = std::atan2(sign * plyCS, plxCS);
    return out;
}

void printEventSummary(const EventSummary& event, int iev) {
    std::cout << "\n[event " << iev << "]\n";
    std::cout << "  category: W" << (event.wCharge > 0 ? "+" : "-") << "\n";
    std::cout << std::setprecision(12)
              << "  lepton: pt=" << event.lepton.Pt()
              << ", eta=" << event.lepton.Eta()
              << ", phi=" << event.lepton.Phi()
              << ", px=" << event.lepton.Px()
              << ", py=" << event.lepton.Py()
              << ", pz=" << event.lepton.Pz()
              << ", E=" << event.lepton.E() << "\n";
    std::cout << "  antilepton: pt=" << event.antilepton.Pt()
              << ", eta=" << event.antilepton.Eta()
              << ", phi=" << event.antilepton.Phi()
              << ", px=" << event.antilepton.Px()
              << ", py=" << event.antilepton.Py()
              << ", pz=" << event.antilepton.Pz()
              << ", E=" << event.antilepton.E() << "\n";
    std::cout << "  TLVUtils:              costh=" << event.tlv.costh
              << ", phi=" << event.tlv.phi << "\n";
    std::cout << "  KinematicDefinitions:  costh=" << event.kd.costh
              << ", phi=" << event.kd.phi << "\n";
    std::cout << "  phi classification:    " << phiRelationName(event.phiRelation)
              << " | coefficient sign flips: "
              << flippedCoefficients(event.phiRelation) << "\n";
    std::cout << "  |delta costh|=" << event.absCosthDiff
              << ", |delta phi|=" << event.absPhiDiff << "\n";
}

void printCategorySummary(const char* name, const CategoryStats& stats) {
    std::cout << "\n[" << name << "]\n";
    std::cout << "events: " << stats.nEvents << "\n";
    std::cout << std::left << std::setw(14) << "Observable"
              << std::setw(18) << "max|diff|"
              << std::setw(12) << "fails"
              << std::setw(14) << "sign flips" << "\n";
    std::cout << std::left << std::setw(14) << "cos(theta)"
              << std::setw(18) << std::scientific << stats.maxAbsCosthDiff
              << std::setw(12) << stats.failCosth
              << std::setw(14) << stats.signFlipCosth << std::defaultfloat << "\n";
    std::cout << std::left << std::setw(14) << "phi"
              << std::setw(18) << std::scientific << stats.maxAbsPhiDiff
              << std::setw(12) << stats.failPhi
              << std::setw(14) << stats.signFlipPhi << std::defaultfloat << "\n";

    std::cout << "phi convention classification:\n";
    std::cout << std::left << std::setw(18) << "relation"
              << std::setw(12) << "events"
              << "coefficients flipping sign" << "\n";
    std::cout << std::left << std::setw(18) << "exact match"
              << std::setw(12) << stats.nExactMatch
              << flippedCoefficients(PhiRelation::ExactMatch) << "\n";
    std::cout << std::left << std::setw(18) << "phi -> -phi"
              << std::setw(12) << stats.nNegPhi
              << flippedCoefficients(PhiRelation::NegPhi) << "\n";
    std::cout << std::left << std::setw(18) << "phi -> phi + pi"
              << std::setw(12) << stats.nPhiPlusPi
              << flippedCoefficients(PhiRelation::PhiPlusPi) << "\n";
    std::cout << std::left << std::setw(18) << "phi -> -phi + pi"
              << std::setw(12) << stats.nNegPhiPlusPi
              << flippedCoefficients(PhiRelation::NegPhiPlusPi) << "\n";
    std::cout << std::left << std::setw(18) << "other"
              << std::setw(12) << stats.nOther
              << flippedCoefficients(PhiRelation::Other) << "\n";

    if (stats.nStoredOtherEvents > 0) {
        std::cout << "sample 'other' events (phi only):\n";
        for (int i = 0; i < stats.nStoredOtherEvents; ++i) {
            const EventSummary& event = stats.sampleOtherEvents[i];
            std::cout << "  [" << i << "] tlv.phi=" << event.tlv.phi
                      << ", kd.phi=" << event.kd.phi
                      << ", wrap(tlv-kd)="
                      << phiDiff(event.tlv.phi, event.kd.phi)
                      << ", wrap(tlv+kd)="
                      << wrapPhi(event.tlv.phi + event.kd.phi)
                      << ", ||wrap(tlv-kd)|-pi|="
                      << std::fabs(std::fabs(phiDiff(event.tlv.phi,
                                                   event.kd.phi)) -
                                   TMath::Pi())
                      << "\n";
        }
    }
}

EventSummary buildEventSummary(const TLorentzVector& lepton,
                               const TLorentzVector& antilepton,
                               int charge1, int inputConvention = 0,
                               double ebeam = 6500., double tol = 1e-5) {
    EventSummary event;
    event.lepton = lepton;
    event.antilepton = antilepton;
    event.wCharge = charge1;

   if (inputConvention == 0) {
        if (charge1 < 0) {
            TLVUtils::getCSFAngles(lepton, charge1, antilepton, ebeam,
                                   event.tlv.costh, event.tlv.phi);
        } else {
            TLVUtils::getCSFAngles(antilepton, charge1, lepton, ebeam,
                                   event.tlv.costh, event.tlv.phi);
        }
    } else if (inputConvention == 1) {
        TLVUtils::getCSFAngles(lepton, charge1, antilepton, ebeam,
                               event.tlv.costh, event.tlv.phi);
    } else if (inputConvention == 2) {
        TLVUtils::getCSFAngles(antilepton, charge1, lepton, ebeam,
                               event.tlv.costh, event.tlv.phi);
    } else {
        std::cout << "Invalid input convention: " << inputConvention
                  << ". Falling back to convention 0.\n";
        if (charge1 < 0) {
            TLVUtils::getCSFAngles(lepton, charge1, antilepton, ebeam,
                                   event.tlv.costh, event.tlv.phi);
        } else {
            TLVUtils::getCSFAngles(antilepton, charge1, lepton, ebeam,
                                   event.tlv.costh, event.tlv.phi);
        }
    }

    event.kd = kinematicDefinitionsAngles(lepton, antilepton);
    event.absCosthDiff = std::fabs(event.tlv.costh - event.kd.costh);
    event.absPhiDiff = std::fabs(phiDiff(event.tlv.phi, event.kd.phi));
    event.costhSignFlip = isSignFlip(event.tlv.costh, event.kd.costh, tol);
    event.phiRelation = classifyPhiRelation(event.tlv.phi, event.kd.phi, tol);
    event.phiSignFlip = event.phiRelation != PhiRelation::ExactMatch &&
                        event.phiRelation != PhiRelation::Other;
    return event;
}

}  // namespace

void RunCoefficientClosureTests(int nEvents = 50000, double tol = 1e-10,
                                int seed = 12345, bool verbose = false,
                                int nVerboseEvents = 5,
                                double ebeam = 6500.,
                                int inputConvention = 0) {
    TRandom3 rngN(seed);
    TRandom3 rngP(seed+33);

    if (nVerboseEvents < 0)
        nVerboseEvents = 0;

    printPhiClassificationSelfCheck(tol);

    int nAccepted = 0;
    int nSkipped = 0;
    CategoryStats statsWMinus;
    CategoryStats statsWPlus;

    TH1D hCosthDiff("hCosthDiff", ";costh_{TLV} - costh_{KD};Events", 200, -2., 2.);
    TH1D hPhiDiff("hPhiDiff", ";wrap(#phi_{TLV} - #phi_{KD});Events", 200, -TMath::Pi(), TMath::Pi());
    TH1D hPhiExact("hPhiExact", ";wrap(#phi_{TLV} - #phi_{KD});Events", 200, -TMath::Pi(), TMath::Pi());
    TH1D hPhiNeg("hPhiNeg", ";wrap(#phi_{TLV} + #phi_{KD});Events", 200, -TMath::Pi(), TMath::Pi());
    TH1D hPhiPlusPi("hPhiPlusPi", ";wrap(#phi_{TLV} - #phi_{KD} - #pi);Events", 200, -TMath::Pi(), TMath::Pi());
    TH1D hPhiNegPlusPi("hPhiNegPlusPi", ";wrap(#phi_{TLV} + #phi_{KD} - #pi);Events", 200, -TMath::Pi(), TMath::Pi());
    TH1D hPhiClass("hPhiClass", ";Classification;Events", 5, 0.5, 5.5);
    TH1D hBosonPt("hBosonPt", ";Boson p_{T} [GeV];Accepted events", 120, 0., 120.);
    TH1D hBosonY("hBosonY", ";Boson Y;Accepted events", 100, -5., 5.);
    TH2D hPhiTLVvsKD("hPhiTLVvsKD", ";#phi_{KD};#phi_{TLV}", 120, -TMath::Pi(), TMath::Pi(),
                     120, -TMath::Pi(), TMath::Pi());

    std::array<TH1D, kNAi> hAiTLV_WMinus;
    std::array<TH1D, kNAi> hAiKD_WMinus;
    std::array<TH1D, kNAi> hAiTLV_WPlus;
    std::array<TH1D, kNAi> hAiKD_WPlus;

    for (int i = 0; i < kNAi; ++i) {
        const TString nameTLVm = Form("hAiTLV_WMinus_%s", kAiNames[i]);
        const TString nameKDm = Form("hAiKD_WMinus_%s", kAiNames[i]);
        const TString nameTLVp = Form("hAiTLV_WPlus_%s", kAiNames[i]);
        const TString nameKDp = Form("hAiKD_WPlus_%s", kAiNames[i]);
        const TString title = Form(";%s;Events", kAiNames[i]);

        hAiTLV_WMinus[i] = TH1D(nameTLVm, title, 120, kAiMin[i], kAiMax[i]);
        hAiKD_WMinus[i] = TH1D(nameKDm, title, 120, kAiMin[i], kAiMax[i]);
        hAiTLV_WPlus[i] = TH1D(nameTLVp, title, 120, kAiMin[i], kAiMax[i]);
        hAiKD_WPlus[i] = TH1D(nameKDp, title, 120, kAiMin[i], kAiMax[i]);

        hAiTLV_WMinus[i].SetLineColor(kBlue + 1);
        hAiTLV_WMinus[i].SetLineWidth(2);
        hAiKD_WMinus[i].SetLineColor(kRed + 1);
        hAiKD_WMinus[i].SetLineWidth(2);

        hAiTLV_WPlus[i].SetLineColor(kBlue + 1);
        hAiTLV_WPlus[i].SetLineWidth(2);
        hAiKD_WPlus[i].SetLineColor(kRed + 1);
        hAiKD_WPlus[i].SetLineWidth(2);

        hAiTLV_WMinus[i].SetStats(0);
        hAiKD_WMinus[i].SetStats(0);
        hAiTLV_WPlus[i].SetStats(0);
        hAiKD_WPlus[i].SetStats(0);
    }

    hPhiClass.GetXaxis()->SetBinLabel(1, "exact");
    hPhiClass.GetXaxis()->SetBinLabel(2, "-phi");
    hPhiClass.GetXaxis()->SetBinLabel(3, "phi+pi");
    hPhiClass.GetXaxis()->SetBinLabel(4, "-phi+pi");
    hPhiClass.GetXaxis()->SetBinLabel(5, "other");

    for (int iev = 0; iev < nEvents; ++iev) {
        const TLorentzVector lep = makeMasslessLepton(rngN);
        const TLorentzVector antilep = makeMasslessLepton(rngP);
        const TLorentzVector boson = lep + antilep;

        if (!isWellDefinedBoson(boson)) {
            ++nSkipped;
            continue;
        }

        ++nAccepted;
        hBosonPt.Fill(boson.Pt());
        hBosonY.Fill(boson.Rapidity());

       
        EventSummary eventWMinus;
        EventSummary eventWPlus;
       
          // W- category: lepton is negatively charged, W+ category: anti-lepton is positively charged
            eventWMinus = buildEventSummary(lep, antilep, -1, inputConvention, ebeam, tol);
            eventWPlus = buildEventSummary(lep, antilep, +1, inputConvention, ebeam, tol);
        if (iev%1000 == 0)
            std::cout << "Processing event " << iev << " / " << nEvents
                      << "   inputConvention " << inputConvention
                      << ", accepted: " << nAccepted
                      << ", skipped: " << nSkipped << ")\n";

        const EventSummary events[2] = {eventWMinus, eventWPlus};
        CategoryStats* stats[2] = {&statsWMinus, &statsWPlus};

        for (int icase = 0; icase < 2; ++icase) {
            const EventSummary& event = events[icase];
            updateCategoryStats(*stats[icase], event, tol);

            const auto aiTLV = computeScaledAi(event.tlv.costh, event.tlv.phi);
            const auto aiKD = computeScaledAi(event.kd.costh, event.kd.phi);
            for (int i = 0; i < kNAi; ++i) {
                if (icase == 0) {
                    hAiTLV_WMinus[i].Fill(aiTLV[i]);
                    hAiKD_WMinus[i].Fill(aiKD[i]);
                } else {
                    hAiTLV_WPlus[i].Fill(aiTLV[i]);
                    hAiKD_WPlus[i].Fill(aiKD[i]);
                }
            }

            const double rawCosthDiff = event.tlv.costh - event.kd.costh;
            const double rawPhiDiff = phiDiff(event.tlv.phi, event.kd.phi);
            const double rawPhiNeg = wrapPhi(event.tlv.phi + event.kd.phi);
            const double rawPhiPlusPi =
                wrapPhi(event.tlv.phi - event.kd.phi - TMath::Pi());
            const double rawPhiNegPlusPi =
                wrapPhi(event.tlv.phi + event.kd.phi - TMath::Pi());

            hCosthDiff.Fill(rawCosthDiff);
            hPhiDiff.Fill(rawPhiDiff);
            hPhiExact.Fill(rawPhiDiff);
            hPhiNeg.Fill(rawPhiNeg);
            hPhiPlusPi.Fill(rawPhiPlusPi);
            hPhiNegPlusPi.Fill(rawPhiNegPlusPi);
            hPhiTLVvsKD.Fill(event.kd.phi, event.tlv.phi);

            switch (event.phiRelation) {
                case PhiRelation::ExactMatch:
                    hPhiClass.Fill(1.0);
                    break;
                case PhiRelation::NegPhi:
                    hPhiClass.Fill(2.0);
                    break;
                case PhiRelation::PhiPlusPi:
                    hPhiClass.Fill(3.0);
                    break;
                case PhiRelation::NegPhiPlusPi:
                    hPhiClass.Fill(4.0);
                    break;
                case PhiRelation::Other:
                    hPhiClass.Fill(5.0);
                    break;
            }
        }

        if (verbose && nAccepted <= nVerboseEvents) {
            std::cout << "\n[W- category]";
            printEventSummary(eventWMinus, iev);
            std::cout << "\n[W+ category]";
            printEventSummary(eventWPlus, iev);
        }
    }

    std::cout << "\n=== Event-by-event Collins-Soper angle comparison ===\n";
    std::cout << "Compared implementations:\n";
    if (inputConvention == 0) {
        std::cout << "  W-: TLVUtils::getCSFAngles(lepton, -1, antilepton, ...), KD(lepton, antilepton)\n";
        std::cout << "  W+: TLVUtils::getCSFAngles(antilepton, +1, lepton, ...), KD(lepton, antilepton)\n";
    } else if (inputConvention == 1) {
        std::cout << "  W-: TLVUtils::getCSFAngles(lepton, -1, antilepton, ...), KD(lepton, antilepton)\n";
        std::cout << "  W+: TLVUtils::getCSFAngles(lepton, +1, antilepton, ...), KD(lepton, antilepton)\n";
    } else if (inputConvention == 2) {
        std::cout << "  W-: TLVUtils::getCSFAngles(antilepton, -1, lepton, ...), KD(lepton, antilepton)\n";
        std::cout << "  W+: TLVUtils::getCSFAngles(antilepton, +1, lepton, ...), KD(lepton, antilepton)\n";
    } else {
        std::cout << "  invalid input convention; TLVUtils falls back to convention 0 internally\n";
        std::cout << "  W-: TLVUtils::getCSFAngles(lepton, -1, antilepton, ...), KD(lepton, antilepton)\n";
        std::cout << "  W+: TLVUtils::getCSFAngles(lepton, +1, antilepton, ...), KD(lepton, antilepton)\n";
    }
    const TString convSuffix = inputConventionSuffix(inputConvention);

    std::cout << "Input convention: "
              << inputConventionLabel(inputConvention)
              << " (" << convSuffix << ")\n\n";
    std::cout << "Generated events: " << nEvents
              << ", accepted: " << nAccepted
              << ", skipped: " << nSkipped
              << ", tolerance: " << std::scientific << tol << std::defaultfloat
              << ", ebeam: " << ebeam
              << ", seed: " << seed << "\n\n";

    std::cout << "\nSign-flip criterion:\n";
    std::cout << "  cos(theta): ||costh_TLV| - |costh_KD|| <= 1e-5\n";
    std::cout << "  phi: classified as phi -> -phi, phi -> phi + pi,"
              << " or phi -> -phi + pi\n";

    printCategorySummary("W- category", statsWMinus);
    printCategorySummary("W+ category", statsWPlus);

    if (statsWMinus.nEvents > 0) {
        std::cout << "\nWorst W- cos(theta) event:";
        printEventSummary(statsWMinus.worstCosth, -1);
        std::cout << "\nWorst W- phi event:";
        printEventSummary(statsWMinus.worstPhi, -1);
    }
    if (statsWPlus.nEvents > 0) {
        std::cout << "\nWorst W+ cos(theta) event:";
        printEventSummary(statsWPlus.worstCosth, -1);
        std::cout << "\nWorst W+ phi event:";
        printEventSummary(statsWPlus.worstPhi, -1);
    }

    gROOT->SetBatch(kTRUE);

    TCanvas cSummary("cSummary", "CS angle validation summary", 1400, 900);
    cSummary.Divide(2, 3);
    cSummary.cd(1);
    hPhiClass.SetStats(0);
    hPhiClass.Draw("hist");
    cSummary.cd(2);
    hPhiTLVvsKD.SetStats(0);
    hPhiTLVvsKD.Draw("colz");
    cSummary.cd(3);
    hCosthDiff.Draw("hist");
    cSummary.cd(4);
    hPhiDiff.Draw("hist");
    cSummary.cd(5);
    hBosonPt.Draw("hist");
    cSummary.cd(6);
    hBosonY.Draw("hist");
    cSummary.SaveAs(Form("CoefficientClosureTests_validation_summary_%s.pdf",
                         convSuffix.Data()));

    TCanvas cHyp("cHyp", "Phi convention hypotheses", 1400, 900);
    cHyp.Divide(2, 2);
    cHyp.cd(1);
    hPhiExact.Draw("hist");
    cHyp.cd(2);
    hPhiNeg.Draw("hist");
    cHyp.cd(3);
    hPhiPlusPi.Draw("hist");
    cHyp.cd(4);
    hPhiNegPlusPi.Draw("hist");
    cHyp.SaveAs(Form("CoefficientClosureTests_phi_hypotheses_%s.pdf",
                     convSuffix.Data()));

    TCanvas cAiWMinus("cAiWMinus", "A0-A7 overlays (W-)", 1600, 900);
    cAiWMinus.Divide(4, 2);
    for (int i = 0; i < kNAi; ++i) {
        cAiWMinus.cd(i + 1);
        const double maxY = std::max(hAiTLV_WMinus[i].GetMaximum(),
                                     hAiKD_WMinus[i].GetMaximum());
        hAiTLV_WMinus[i].SetMaximum(maxY * 1.2);
        hAiTLV_WMinus[i].Draw("hist");
        hAiKD_WMinus[i].Draw("hist same");
        if (i == 0) {
            TLegend* leg = new TLegend(0.55, 0.72, 0.89, 0.89);
            leg->SetBorderSize(0);
            leg->AddEntry(&hAiTLV_WMinus[i], "TLVUtils angles", "l");
            leg->AddEntry(&hAiKD_WMinus[i], "KinematicDefinitions angles", "l");
            leg->Draw();
        }
    }
    cAiWMinus.SaveAs(Form("CoefficientClosureTests_Ai_overlay_WMinus_%s.pdf",
                          convSuffix.Data()));

    TCanvas cAiWPlus("cAiWPlus", "A0-A7 overlays (W+)", 1600, 900);
    cAiWPlus.Divide(4, 2);
    for (int i = 0; i < kNAi; ++i) {
        cAiWPlus.cd(i + 1);
        const double maxY = std::max(hAiTLV_WPlus[i].GetMaximum(),
                                     hAiKD_WPlus[i].GetMaximum());
        hAiTLV_WPlus[i].SetMaximum(maxY * 1.2);
        hAiTLV_WPlus[i].Draw("hist");
        hAiKD_WPlus[i].Draw("hist same");
        if (i == 0) {
            TLegend* leg = new TLegend(0.55, 0.72, 0.89, 0.89);
            leg->SetBorderSize(0);
            leg->AddEntry(&hAiTLV_WPlus[i], "TLVUtils angles", "l");
            leg->AddEntry(&hAiKD_WPlus[i], "KinematicDefinitions angles", "l");
            leg->Draw();
        }
    }
    cAiWPlus.SaveAs(Form("CoefficientClosureTests_Ai_overlay_WPlus_%s.pdf",
                         convSuffix.Data()));

    std::cout << "\nValidation plots written to:\n";
    std::cout << "  CoefficientClosureTests_validation_summary_" << convSuffix
              << ".pdf\n";
    std::cout << "  CoefficientClosureTests_phi_hypotheses_" << convSuffix
              << ".pdf\n";
    std::cout << "  CoefficientClosureTests_Ai_overlay_WMinus_" << convSuffix
              << ".pdf\n";
    std::cout << "  CoefficientClosureTests_Ai_overlay_WPlus_" << convSuffix
              << ".pdf\n";

    std::cout << "\nDone.\n";
}

void CoefficientClosureTests(const char* call = "") {
    if (call && std::strlen(call) > 0) {
        gROOT->ProcessLine(call);
        return;
    }
    RunCoefficientClosureTests();
}
