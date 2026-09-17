# Explanation of `AIZ.C`

`AIZ.C` is a ROOT macro for a truth-level Z->ee angular-coefficient analysis. It reads a Z electron ntuple, rebuilds the electron and Z four-vectors, computes Collins-Soper angles, derives the angular coefficients `A0` to `A7`, and writes the results plus many diagnostic histograms to a ROOT output file.

## Main Purpose

The macro measures the angular coefficients `Ai` for Z boson production and decay:

```text
Z -> e+ e-
```

The coefficients are evaluated either as a function of:

- Z transverse momentum, `pT(Z)`, with `AIZ(false)`, which is the default.
- Absolute Z rapidity, `|y(Z)|`, with `AIZ(true)`.

## Default Configuration

At the top of the macro, the main switches are:

```cpp
bool sherpa = false;
bool test = true;
bool override = false;
bool normXS = true;
bool ifTrueOnly = true;
bool FiducialCut = true;
bool FiducialCutEtaonly = false;
bool FiducialCutCCCF = false;
bool FiducialCutCFonly = true;
int polynomialIndex = 6;
bool appendPolynomialOutputs = true;
```

With these defaults:

- The macro uses the truth-only Powheg Z->ee ntuple.
- Sherpa is not active because no Sherpa Z sample path is defined.
- It applies the CF-only fiducial selection because `FiducialCutCFonly = true`:
   both leptons must have `pT > 25 GeV`, with one central lepton (`|eta| < 2.5`)
   and one forward lepton (`|eta| >= 2.5`).
- It processes at most 100,000 events because `test = true`. Set `test = false`
   for full statistics.
- It refuses to overwrite an existing output file unless `override = true`.

## Input and Output

The input tree is:

```cpp
HWWTree_ee
```

The macro also reads:

- `MetaData`, to get cross section, k-factor, and filter efficiency.
- `CutFlow`, to get normalization information.

The output file name is built from the configuration. With the current defaults and
`AIZ(false)` in test mode, it is:

```text
AI_Z_Truth_Fiducial_CFonly_testPowheg_pT_NormXsec.root
```

Set `test = false` to replace `testPowheg` with `Zai_finalbinningPowheg`.
For `AIZ(true)`, the output name uses `_Y_` instead of `_pT_`.

## Event Loop

For each event, the macro:

1. Reads truth lepton kinematics:
   - `lepPtTruth0`, `lepEtaTruth0`, `lepPhiTruth0`
   - `lepPtTruth1`, `lepEtaTruth1`, `lepPhiTruth1`
   - `mcEventWeight`

2. Builds two `TLorentzVector` objects:
   - `em` for the negative electron, `e-`
   - `ep` for the positive electron, `e+`

3. Builds the Z boson:

   ```cpp
   z = em + ep;
   ```

4. Applies the Z mass window:

   ```text
   66 GeV < m(e+e-) < 116 GeV
   ```

5. Computes the Collins-Soper angles using `TLVUtils`:

   ```cpp
   TLVUtils::getCSFAngles(em, -1, ep, ebeamGeV, costheta, phi);
   ```

   Here:

   - `costheta` is `cos(theta_CS)`.
   - `phi` is `phi_CS`, normalized to the ATLAS range $[0,2\pi)$. A signed value in $[-\pi,\pi)$ is equivalent after adding $2\pi$ when negative. This does not change the periodic angular basis, but it does set the histogram and plot axes. Laboratory lepton `phi` values remain in ROOT's signed convention.
   - `ebeamGeV = 6500.0`, corresponding to 13 TeV pp collisions.

## Angular Coefficients

After computing `cos(theta_CS)` and `phi_CS`, the macro calls:

```cpp
TLVUtils::getAiPolynoms(costheta, phi, aipols);
```

This returns the angular polynomial basis terms. The macro then converts them into the historical `AIZ.C` coefficient convention:

```cpp
A0 = 20/3 * p0 + 2/3
A1 = 5    * p1
A2 = 10   * p2
A3 = 4    * p3
A4 = 4    * p4
A5 = 5    * p5
A6 = 4    * p6
A7 = 4    * p7
```

Each coefficient is filled into a histogram binned in either `pT(Z)` or `|y(Z)|`.

At the end of the event loop, each coefficient histogram is divided by the weighted event-count histogram:

```cpp
A0->Divide(Xsw);
A1->Divide(Xsw);
...
A7->Divide(Xsw);
```

So each bin contains the weighted average value of the corresponding angular coefficient.

## Configurable Polynomial Diagnostics

The diagnostic plots can use any angular basis polynomial returned by
`TLVUtils::getAiPolynoms`, indexed from 0 through 7:

```cpp
int polynomialIndex = 6;
```

`polynomialIndex = 6` preserves the historical default. Before calling `AIZ`, select
another polynomial, for example:

```cpp
polynomialIndex = 2;
AIZ(true);             // use |y(Z)| binning
```

Alternatively, configure it directly in the `AIZ` call:

```cpp
AIZ(true, 2);          // use P2 and |y(Z)| binning
```

The second argument is optional; existing calls such as `AIZ(true)` continue to use
the global `polynomialIndex` value.

The selected index controls all three related diagnostic PDFs:

- `*_P2_polynomial.pdf`: analytic `P2` basis, selected-event angular map, and slices.
- `*_P2_observables.pdf`: selected-event `P2` distributions and correlations.
- `*_P2_sensitivity.pdf`: selected-event sign asymmetry and squared-polynomial lever arms.

Replace `P2` with the selected index. The same index is also used in the corresponding
ROOT histogram names, such as `P2`, `P2_vs_deltaEta`, and `P2_moment_vs_deltaEta`.
Values outside 0 through 7 are rejected. These configurable diagnostics are separate
from the legacy `A0` through `A7` coefficient histograms, which are always produced
using the historical scaling convention described above.

When `appendPolynomialOutputs = true`, existing output ROOT files are opened in
`UPDATE` mode. This allows repeated runs with different `polynomialIndex` values to
keep all `P0` through `P7` diagnostic objects in one file. Repeating an index updates
that index's keys with `TObject::kOverwrite`, avoiding duplicate cycles for the
polynomial-specific objects. Set it to `false` when a complete replacement of the
output file is desired.

## Histograms Written

The macro writes the main coefficient histograms:

- `A0`, `A1`, `A2`, `A3`
- `A4`, `A5`, `A6`, `A7`
- `sigma`, the weighted event yield histogram
- `sigma_unweighted`, the unweighted event yield histogram

It also writes many diagnostic histograms, including:

- `Costheta`
- `phi`
- `Zmass`
- Z `pT` versus Collins-Soper angles
- Z rapidity versus Collins-Soper angles
- lepton `eta`, `phi`, and `pT` correlations
- leading and subleading lepton distributions
- `DeltaR`, `DeltaPhi`, and `DeltaEta` between the leptons
- opening-angle distributions between the two leptons
- cos(theta) slices in lepton `pT` bins

For the selected `polynomialIndex`, it additionally writes the polynomial diagnostic
histograms and canvases used by the three configurable PDF outputs.

The macro also writes two comparison canvases:

- `c_costheta_comp`
- `c_phi_comp`

These compare angles computed from the lepton four-vectors against truth-angle branches, when those branches are available.

## Important Caveat

With the current default:

```cpp
ifTrueOnly = true;
```

the macro does not read the branches:

```cpp
cosThetaCSTruth
phiCSTruth
```

Those variables remain initialized to zero, but the macro still fills:

```cpp
CosthetaTruth
phiTruth
```

Therefore, in truth-only mode, the comparison canvases `c_costheta_comp` and `c_phi_comp` are probably not meaningful. The angular coefficients themselves are still based on the Collins-Soper angles recomputed from the lepton four-vectors, so the `A0` to `A7` histograms use the intended `TLVUtils` calculation.

## How to Run

Inside ROOT:

```cpp
.L TLVUtils.cxx
.L AIZ.C
AIZ(false);   // coefficients vs pT(Z)
AIZ(true);    // coefficients vs |y(Z)|
```

For a quick test:

```cpp
.L TLVUtils.cxx
.L AIZ.C
test = true;
polynomialIndex = 2;  // optional: use P2 for the diagnostic plots
AIZ(false);
```

If the output file already exists:

```cpp
override = true;
```

must be set before running.

## Short Summary

`AIZ.C` produces a ROOT file containing Z->ee Collins-Soper angular coefficients `A0` to `A7`, binned by either Z `pT` or `|y(Z)|`. It uses `TLVUtils` for the Collins-Soper angle and angular-polynomial definitions, applies a Z mass window and optional fiducial lepton cuts, and saves both the coefficient histograms and a broad set of kinematic control plots.
