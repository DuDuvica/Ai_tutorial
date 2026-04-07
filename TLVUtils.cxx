#ifndef TLVUtils_cxx
#define TLVUtils_cxx

#include "TLVUtils.h"

#include <iostream>

namespace TLVUtils {

// masses in MeV
constexpr double ELE_MASS = 0.51099895069;
constexpr double MU_MASS = 105.6583755;

double deltaPhimin(const TLorentzVector &tlv,
                   const vector<TLorentzVector> &vtlv, int &index) {

    double dphimin = 666;
    for (unsigned int i = 0; i < vtlv.size(); i++) {
        double dphi = tlv.DeltaPhi(vtlv[i]);
        if (dphi < dphimin) {
            dphimin = dphi;
            index = i;
        }
    }
    return dphimin;
}

double deltaRmin(const TLorentzVector &tlv, const vector<TLorentzVector> &vtlv,
                 int &index, const float min) {

    double drmin = 666;
    for (unsigned int i = 0; i < vtlv.size(); i++) {
        double dr = tlv.DeltaR(vtlv[i]);
        if (dr < drmin && dr > min) {
            drmin = dr;
            index = i;
        }
    }
    return drmin;
}

double deltaRmin(const vector<TLorentzVector> &vtlv, int &index1, int &index2) {

    double drmin = 666;
    for (unsigned int i = 0; i < vtlv.size(); i++) {
        for (unsigned int j = i + 1; j < vtlv.size(); j++) {
            double dr = vtlv[i].DeltaR(vtlv[j]);
            if (dr < drmin) {
                drmin = dr;
                index1 = i;
                index2 = j;
            }
        }
    }
    return drmin;
}

void remove_if_not_in(vector<int> &values, const vector<int> &indexes) {

    vector<int>::iterator idx = values.begin();
    while (idx != values.end()) {
        vector<int>::const_iterator pos =
            find(indexes.begin(), indexes.end(), *idx);
        if (pos == indexes.end())  // index not in indexes => has to be removed
            idx = values.erase(
                idx);  // then idx points to the element after the erased one
        else
            ++idx;
    }
}

// Generic functions to sort vector<TLorentzVector>

vector<int> rank(const vector<float> &vecfloat, bool decreasing) {

    vector<pair<int, float> > vecpair;
    for (unsigned int i = 0; i < vecfloat.size(); i++)
        vecpair.push_back(pair<int, float>(i, vecfloat[i]));

    if (decreasing)
        sort(vecpair.begin(), vecpair.end(), larger());
    else
        sort(vecpair.begin(), vecpair.end(), smaller());

    vector<int> vecint;
    for (unsigned int i = 0; i < vecfloat.size(); i++)
        vecint.push_back(vecpair[i].first);

    return vecint;
}

vector<int> rankPt(const vector<TLorentzVector> &vectlv, bool decreasing) {

    vector<float> vecfloat;
    for (unsigned int i = 0; i < vectlv.size(); i++) {
        vecfloat.push_back(vectlv[i].Pt());
    }

    return rank(vecfloat, decreasing);
}

vector<int> rankEta(const vector<TLorentzVector> &vectlv, bool decreasing) {

    vector<float> vecfloat;
    for (unsigned int i = 0; i < vectlv.size(); i++)
        vecfloat.push_back(vectlv[i].Eta());

    return rank(vecfloat, decreasing);
}

vector<int> rankMass(const vector<TLorentzVector> &vectlv, bool decreasing) {

    vector<float> vecfloat;
    for (unsigned int i = 0; i < vectlv.size(); i++)
        vecfloat.push_back(vectlv[i].M());

    return rank(vecfloat, decreasing);
}

void getBDFAngles(const TLorentzVector &lep1, const int &charge1,
                  const TLorentzVector &lep2, double &costh, double &phi) {

    const TLorentzVector dilep = lep1 + lep2;

    double cosThetap = (charge1 > 0)
                           ? TMath::TanH((lep2.Eta() - lep1.Eta()) / 2.)
                           : TMath::TanH((lep1.Eta() - lep2.Eta()) / 2.);
    costh = (dilep.Rapidity() < 0) ? -cosThetap : cosThetap;

    double sintheta = TMath::Sqrt(1. - costh * costh);
    double phiacop =
        TMath::Pi() - TVector2::Phi_mpi_pi(lep1.Phi() - lep2.Phi());
    phi = TMath::Tan(phiacop / 2) * sintheta;
}

void getVHFAngles(const TLorentzVector &lep1, const int &charge1,
                  const TLorentzVector &lep2, double &costh, double &phi) {

    TLorentzVector boson = lep1 + lep2;
    TVector3 boostV = boson.BoostVector();
    TLorentzVector lep1_boosted = (charge1 > 0) ? lep2 : lep1;
    lep1_boosted.Boost(-boostV);
    phi = TVector2::Phi_mpi_pi(lep1_boosted.Phi() - boson.Phi());
    double theta = lep1_boosted.Angle(boson.Vect());

    if (charge1 == 0)
        std::cout << "[93mCharge1 est nulle " << charge1 << " AH AH AH AH[0m"
                  << std::endl;

    costh = TMath::Cos(theta);
}

void getCSFAngles(const TLorentzVector &lep1, const int &charge1,
                  const TLorentzVector &lep2, double ebeam, double &costh,
                  double &phi) {

    TLorentzVector boson = lep1 + lep2;
    double Lplus = (charge1 < 0) ? lep1.E() + lep1.Pz() : lep2.E() + lep2.Pz();
    double Lminus = (charge1 < 0) ? lep1.E() - lep1.Pz() : lep2.E() - lep2.Pz();
    double Pplus = (charge1 < 0) ? lep2.E() + lep2.Pz() : lep1.E() + lep1.Pz();
    double Pminus = (charge1 < 0) ? lep2.E() - lep2.Pz() : lep1.E() - lep1.Pz();

    costh = (Lplus * Pminus - Lminus * Pplus);
    costh *= TMath::Abs(boson.Pz());
    costh /= (boson.Mag() * boson.Pz());
    costh /= TMath::Sqrt(boson.Mag2() + boson.Pt() * boson.Pt());

    if (costh > 1)
        costh = 1.;
    if (costh < -1)
        costh = -1.;

    TVector3 boostV = -boson.BoostVector();
    TLorentzVector lep1_boosted = (charge1 < 0) ? lep1 : lep2;
    lep1_boosted.Boost(boostV);

    TVector3 CSAxis, xAxis, yAxis;
    TLorentzVector p1, p2;
    double sign = +1.;
    if (boson.Z() < 0)
        sign = -1.;
    p1.SetXYZM(0., 0., sign * ebeam, 0.938);   // quark (?)
    p2.SetXYZM(0., 0., -sign * ebeam, 0.938);  // antiquark (?)

    p1.Boost(boostV);
    p2.Boost(boostV);
    CSAxis = (p1.Vect().Unit() - p2.Vect().Unit()).Unit();
    yAxis = (p1.Vect().Unit()).Cross(p2.Vect().Unit());
    yAxis = yAxis.Unit();
    xAxis = yAxis.Cross(CSAxis);
    xAxis = xAxis.Unit();

    phi =
        TMath::ATan2(lep1_boosted.Vect() * yAxis, lep1_boosted.Vect() * xAxis);

    //    std::cout << TMath::ACos(CSAxis) << std::endl;
}

void getProjections(const TLorentzVector &proj, const TLorentzVector &axis,
                    double &par, double &perp) {

    par = (proj.Px() * axis.Px() + proj.Py() * axis.Py()) / axis.Pt();
    perp = (proj.Py() * axis.Px() - proj.Px() * axis.Py()) / axis.Pt();
}

void getAiPolynoms(double costh, double phi, std::vector<double> &aipols) {

    aipols.resize(8);

    double theta = TMath::ACos(costh);
    double sintheta = TMath::Sin(theta);
    double sin2theta = TMath::Sin(2.0 * theta);

    double cosphi = TMath::Cos(phi);
    double cos2phi = TMath::Cos(2.0 * phi);
    double sinphi = TMath::Sin(phi);
    double sin2phi = TMath::Sin(2.0 * phi);

    aipols[0] = 0.5 - 1.5 * costh * costh;
    aipols[1] = sin2theta * cosphi;
    aipols[2] = sintheta * sintheta * cos2phi;
    aipols[3] = sintheta * cosphi;
    aipols[4] = costh;
    aipols[5] = sintheta * sintheta * sin2phi;
    aipols[6] = sin2theta * sinphi;
    aipols[7] = sintheta * sinphi;
}

double CosThetaBeamSystem(const TLorentzVector &lep1,
                          const TLorentzVector &lep2) {
    double costheta = TMath::TanH(0.5 * (lep1.Eta() - lep2.Eta()));
    if ((lep1.Pz() + lep2.Pz()) < 0)
        costheta *= -1.;
    return costheta;
}

double CosThetaRestSystem(const TLorentzVector &lep1,
                          const TLorentzVector &lep2) {
    TLorentzVector boson = lep1 + lep2;
    TVector3 boostV = boson.BoostVector();
    TLorentzVector lep1_boosted = lep1;
    lep1_boosted.Boost(-boostV);
    return (TMath::Cos(lep1_boosted.Angle(boson.Vect())));
}

double PhiStar(const TLorentzVector &lep1, const TLorentzVector &lep2) {
    double costheta = CosThetaBeamSystem(lep1, lep2);
    double sintheta = TMath::Sqrt(1. - costheta * costheta);
    double phiacop =
        TMath::Pi() - TVector2::Phi_mpi_pi(lep1.Phi() - lep2.Phi());
    return (TMath::Tan(phiacop / 2) * sintheta);
}

double AiPolSum(const std::vector<double> &ai,
                const std::vector<double> &aimom) {
    static const double norm = 3. / (16. * TMath::Pi());
    double xsec = 1 + aimom[4] * aimom[4];
    for (int imom = 0; imom < 8; imom++)
        xsec += ((imom == 2) ? 0.5 : 1) * ai[imom] * aimom[imom];
    return xsec * norm;
}

float deltaR2(const float &eta1, const float &phi1, const float &eta2,
              const float &phi2) {
    float deta = eta1 - eta2;
    float dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return deta * deta + dphi * dphi;
}

std::optional<TLorentzVector> getMatchedTruthTLV(
    const std::vector<float> *tlep_pt, const std::vector<float> *tlep_eta,
    const std::vector<float> *tlep_phi, const std::vector<int> *tlep_pdgId,
    const float &rec_eta, const float &rec_phi, const int &pdgId, float maxDR) {
    // check that all vectors are of the same size
    if (!tlep_pt || !tlep_eta || !tlep_phi || !tlep_pdgId ||
        tlep_pt->size() != tlep_eta->size() ||
        tlep_pt->size() != tlep_phi->size() ||
        tlep_pt->size() != tlep_pdgId->size())
        return {};

    float minDR = 999.;
    int i_minDR = -1;
    for (unsigned int i = 0; i < tlep_pt->size(); i++) {
        if (tlep_pdgId->at(i) != pdgId)
            continue;
        float dr2 = deltaR2(rec_eta, rec_phi, tlep_eta->at(i), tlep_phi->at(i));
        if (dr2 < minDR) {
            i_minDR = i;
            minDR = dr2;
        }
    }

    minDR = sqrt(minDR);
    if (i_minDR >= 0 && minDR < maxDR) {
        TLorentzVector tlv_truth;
        int i = i_minDR;
        float tlep_m = (fabs(tlep_pdgId->at(i)) == 11) ? ELE_MASS : MU_MASS;
        tlv_truth.SetPtEtaPhiM(tlep_pt->at(i), tlep_eta->at(i), tlep_phi->at(i),
                               tlep_m);
        return tlv_truth;
    } else {
        return {};
    }
}

int getMatchedTruthPDGID(const std::vector<float> *tlep_pt,
                         const std::vector<float> *tlep_eta,
                         const std::vector<float> *tlep_phi,
                         const std::vector<int> *tlep_pdgId,
                         const float &rec_eta, const float &rec_phi,
                         float maxDR) {
    // check that all vectors are of the same size
    if (!tlep_pt || !tlep_eta || !tlep_phi || !tlep_pdgId ||
        tlep_pt->size() != tlep_eta->size() ||
        tlep_pt->size() != tlep_phi->size())  //||
        return {};

    float minDR = 999.;
    int i_minDR = -1;
    for (unsigned int i = 0; i < tlep_pt->size(); i++) {
        float dr2 = deltaR2(rec_eta, rec_phi, tlep_eta->at(i), tlep_phi->at(i));
        if (dr2 < minDR) {
            i_minDR = i;
            minDR = dr2;
        }
    }

    minDR = sqrt(minDR);
    if (i_minDR >= 0 && minDR < maxDR) {
        TLorentzVector tlv_truth;
        int i = i_minDR;
        int pdgID = tlep_pdgId->at(i);
        return pdgID;
    } else {
        return {};
    }
}

}  // namespace TLVUtils

#endif
