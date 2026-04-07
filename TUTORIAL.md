# AIW.C Tutorial (W) + How To Adapt For Z

This tutorial explains how to run `AIW.C`, what the main switches do, what
outputs it produces, and what to change to make a Z-boson version.

## What The Macro Does

`AIW.C` reads a MicroTree ROOT file, computes Collins–Soper angles
(`theta`, `phi`) for W->lnu events, and fills the angular coefficients
`A0–A7` as a function of W `pT` (or rapidity if `isY=true`). It also
produces weighted and unweighted cross-section histograms.

Important: `A0–A7` are `TH1D` objects (1D histograms), not `TH2`.
So each fill must use `Fill(x, weight)`, where the histogram weight is
the product of the angular basis term and the event weight.

Outputs written to a ROOT file:
- `Costheta`, `phi`
- `A0` ... `A7` (normalized by `sigma`)
- `sigma` (weighted cross-section histogram)
- `sigma_unweighted`

## Quick Run

From a ROOT prompt:

```bash
root -l
root [0] .L AIW.C+
root [1] AIW(false, false)  // W- vs pT (default)
root [2] AIW(true,  false)  // W+ vs pT
root [3] AIW(false, true)   // W- vs |Y|
```

Or non-interactive:

```bash
root -l -b -q 'AIW.C+(false,false)'
```

## Key Configuration Flags (Top Of File)

Edit these globals in `AIW.C` before running:
- `sherpa`: use Sherpa samples if `true`, Powheg otherwise.
- `testSherpaNEW`: switches Sherpa input files to the newer 2025 paths.
- `newPwg`: toggles newer Powheg paths.
- `normXS`: if `true`, normalize by AMI cross-sections.
- `signSwop`: optional sign flip for W+ (used for older Sherpa sign issues).
- `overight`: set `true` to overwrite an existing output file.

The input file is selected by:
- `sherpa` vs `!sherpa`
- `isPlus` argument (W+ or W-)
- `testSherpaNEW` / `newPwg`

The output file name is built from:
- charge: `Wm` or `Wp`
- `outputName` string
- mode: `_pT` or `_Y` and `_NormXsec` suffix when enabled

## Common Edits

### 1) Use different binning
Edit the `bins` array near:

```cpp
const int Nbins = 11;
Double_t bins[Nbins] = {0., 8., 17., 27., 40., 55., 75., 110., 150., 210., 600.};
```

### 2) Use rapidity binning
Call `AIW(isPlus, true)` or set `isY=true`.
The macro then builds uniform `|Y|` bins of width 0.4.

### 3) Change normalization
- `xsecAMI` is set per generator and charge.
- `normXS=false` to use raw generator weights only.

## How To Adapt For Z Bosons

Create a copy of the macro (recommended: `AIZ.C`) and make these changes:

### 1) Input file and tree name
Replace W files with Z files and the tree path:

```cpp
TTree *tree = (TTree*)f1.Get("MicroTree/HWWTree_ee");
```

If your Z MicroTree uses a different tree path, update accordingly.

### 2) Use two charged leptons
Replace neutrino usage with the second lepton:

```cpp
TLorentzVector em, ep, z, delta;

em.SetPtEtaPhiM(lepPtTruth1/1000.0, lepEtaTruth1, lepPhiTruth1, 0);
// use the other lepton as e+
ep.SetPtEtaPhiM(lepPtTruth0/1000.0, lepEtaTruth0, lepPhiTruth0, 0);

z = em + ep;

delta = em - ep; // negative minus positive lepton
```

### 3) Update the Collins–Soper `cos(theta)` numerator
For Z:

```cpp
lpe  = em.E() + em.Pz();
lme  = em.E() - em.Pz();
lpeplus = ep.E() + ep.Pz();
lmeplus = ep.E() - ep.Pz();

costheta = (1./z.M()) * pow((z.M2()+pow(z.Pt(),2)),-0.5)
         * (lpeplus*lme - lmeplus*lpe);
```

### 4) Replace W kinematics with Z kinematics
Use `z` in place of `w` in:
- `costheta` normalization
- `rapidity` sign
- `pT` or `|Y|` fills
- `qttrans` computation

### 5) Cross-section normalization
Update `xsecAMI` for Z samples. If not needed, set `normXS=false`.

### 6) Output naming
Change the output prefix to something like:

```cpp
TString nameOutput = "AI_Z_" + outputName + mode + ".root";
```

## Expected Outputs

Example output files for W+:
- `AI_WWp_Wai_finalbinning_newSignSherpa221NEW_pT_NormXsec.root`

Inside the file:
- `A0` ... `A7` (normalized)
- `sigma` and `sigma_unweighted`
- `Costheta`, `phi`

## Notes

- Explicit bug fix (TH1D fill convention): `A0–A7` fills were updated from
  `Fill(x, angular_term, weight)` to `Fill(x, angular_term * weight)` in
  both `AIW.C` and `AIZ.C`.
- This is required because `A0–A7` are `TH1D` histograms and accept
  `Fill(x, weight)` (2 arguments), not 3-argument `TH2`-style fills.
- The macro now weights A0–A7 with the same event weight used for `sigma`.
- `signSwop` is only applied for W+ and should be used only if you know
  the generator sign convention requires it.

