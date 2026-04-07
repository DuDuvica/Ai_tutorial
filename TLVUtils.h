#ifndef TLVUtils_h
#define TLVUtils_h

#include <TLorentzVector.h>

#include <optional>

using namespace std;

namespace TLVUtils {

// Functions to compute deltaRs

// min phi  distance between a tlv and a collection of tlvs
// returns the min distance
// index will point to the position of the tlv that minimizes the phi distance
double deltaPhimin(const TLorentzVector& tlv,
                   const vector<TLorentzVector>& vtlv, int& index);

// min distance between a tlv and a collection of tlvs
// returns the min distance
// index will point to the position of the tlv that minimizes the distance
double deltaRmin(const TLorentzVector& tlv, const vector<TLorentzVector>& vtlv,
                 int& index, const float min = -1);

// min distance among a collection of tlvs
// returns the min distance
// index1 and index2 will point to the positions of the tlvs that minimizes the
// distance NB: index2>index1
double deltaRmin(const vector<TLorentzVector>& vtlv, int& index1, int& index2);

// Return decay angles in the Collins-Soper Frame
// Require tlv of both decay-product particles
// Return cos(theta_{CS}) and phi_{CS} (passed by reference in argument)
// Return Ai moment computation is pointer to vector given (passed by reference
// by pointer)
void getCSFAngles(const TLorentzVector& lep1, const int& charge1,
                  const TLorentzVector& lep2, double ebeam, double& costh,
                  double& phi);

// Return decay angles in the Vector boson Helicity Frame
// Require tlv of both decay-product particles
// Return cos(theta_{VH}) and phi_{VH} (passed by reference in argument)
// Return Ai moment computation is pointer to vector given (passed by reference
// by pointer)
void getVHFAngles(const TLorentzVector& lep1, const int& charge1,
                  const TLorentzVector& lep2, double& costh, double& phi);

// calculate angular polynomials
void getAiPolynoms(double costh, double phi, std::vector<double>& aipols);

// Return decay angles in the Beam Direction Frame
// Require tlv of both decay-product particles
// Return cos(theta_{BD}) and phi_{BD} (passed by reference in argument)
void getBDFAngles(const TLorentzVector& tlv1, const int& charge1,
                  const TLorentzVector& tlv2, double& costh, double& phi);

// projections
void getProjections(const TLorentzVector& proj, const TLorentzVector& axis,
                    double& par, double& perp);

// polarisation helper functions
double CosThetaBeamSystem(const TLorentzVector& lep1,
                          const TLorentzVector& lep2);
double CosThetaRestSystem(const TLorentzVector& lep1,
                          const TLorentzVector& lep2);
double PhiStar(const TLorentzVector& lep1, const TLorentzVector& lep2);
double AiPolSum(const std::vector<double>& ai,
                const std::vector<double>& aimom);

// remove from values all elements that are not present in indexes.
void remove_if_not_in(vector<int>& values, const vector<int>& indexes);

// Generic functions to sort vector<TLorentzVector>
vector<int> rank(const vector<float>& vecfloat, bool decreasing = true);
vector<int> rankPt(const vector<TLorentzVector>& vectlv,
                   bool decreasing = true);
vector<int> rankEta(const vector<TLorentzVector>& vectlv,
                    bool decreasing = true);
vector<int> rankMass(const vector<TLorentzVector>& vectlv,
                     bool decreasing = true);

struct larger {
    bool operator()(const pair<int, float>& p1, const pair<int, float>& p2) {
        return p1.second > p2.second;
    }
};

struct smaller {
    bool operator()(const pair<int, float>& p1, const pair<int, float>& p2) {
        return p2.second > p1.second;
    }
};

float deltaR2(const float& eta1, const float& phi1, const float& eta2,
              const float& phi2);

// a function for reco-truth matching
std::optional<TLorentzVector> getMatchedTruthTLV(
    const std::vector<float>* tlep_pt, const std::vector<float>* tlep_eta,
    const std::vector<float>* tlep_phi, const std::vector<int>* tlep_pdgId,
    const float& rec_eta, const float& rec_phi, const int& pdgId,
    float maxDR = 0.2);

// a function for reco-truth matching for electron charge flips
int getMatchedTruthPDGID(const std::vector<float>* tlep_pt,
                         const std::vector<float>* tlep_eta,
                         const std::vector<float>* tlep_phi,
                         const std::vector<int>* tlep_pdgId,
                         const float& rec_eta, const float& rec_phi,
                         float maxDR = 0.2);

}  // namespace TLVUtils

#endif
