# Ai_tutorial

Quick plotting utilities for `Ai`.

## Quick Start (3 commands)

```bash
conda activate MAC && cd /data/dust/user/ludovica/CraigMacro/Tutorial_AI/Ai_tutorial
root -l -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'override=true;' -e 'test=true;' -e 'ifTrueOnly=true;' -e 'FiducialCut=false;' -e 'FiducialCutEtaonly=false;' -e 'normXS=true;' -e 'AIZ(true);'
root -l -q 'CompareAIZProjections.C("AI_Z_Truth_testPowheg_Y_NormXsec.root", true)'
```

The polynomial diagnostics default to `polynomialIndex = 6`. Set it to any value from
0 to 7 either before calling `AIZ`, or pass it directly as the second argument:

```bash
root -l -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'override=true;' -e 'test=false;' -e 'AIZ(true, 4);'
```

This writes `P2_polynomial.pdf`, `P2_observables.pdf`, and `P2_sensitivity.pdf`
using the selected polynomial. The existing `A0` through `A7` coefficient histograms
are still produced independently.

By default, `appendPolynomialOutputs = true`, so running `AIZ` repeatedly with
different indices keeps all `P0` through `P7` diagnostic histograms and canvases in
the same ROOT file. Repeating an index updates that index's keys instead of creating
duplicate cycles. Set `appendPolynomialOutputs = false` to restore replace-file behavior.

## Files

- `AIW.C`: macro used to produce the reference `Ai` for the MC.
- `AIZ.C`: related plotting macro.
- `TUTORIAL.md`: tutorial notes for running the workflow.

## Reference Input

The reference `Ai` for the MC (`AI_WWm_Wai_finalbinning_newSignPhowegNEW_pT_NormXsec.root`) is produced on NAF in:

`/nfs/dust/atlas/user/ludovica/CraigMacro`

using the `AIW.C` macro.

## Run AIZ (ROOT + TLVUtils)

### 1) Activate ROOT environment

```bash
conda activate MAC
```

### 2) Go to the tutorial folder

```bash
cd /Users/dudu/Downloads/Low_muRun/WAi/predictions_for_W_Ai_13TeV/Phoweg_Sherpa_comparisonPLOTs/Tutorial_AI
```

### 3) Run in ROOT (recommended interpreted mode)

```bash
root -l
```

Then inside ROOT:

```cpp
.L TLVUtils.cxx
.L AIZ.C
test = true;      // optional: quick run
AIZ(false);       // pT binning
// AIZ(true);     // |y| binning
```

To focus the A4 sensitivity-loss PDF on the bulk of its distribution, set
`zoomA4SensitivityLoss = true` before calling `AIZ`. The zoomed plot is capped at
20 on the y-axis, so larger values are clipped; leave the flag false for the
full-range plot.

### Interpreting and comparing A4 sensitivity-loss plots

The plotted factor is the estimated statistical uncertainty on `A4` after the
fiducial cut divided by that before the cut:

`SensitivityLossFactor = sigma(A4)_after / sigma(A4)_before`.

In this estimate, the information per bin is approximated by
`N_eff * <P4^2>`, so
`SensitivityLossFactor = sqrt((N_eff_before * <P4^2>_before) / (N_eff_after * <P4^2>_after))`.
Fewer effective events or a smaller `P4` lever arm means less information and
therefore a larger uncertainty ratio.

The horizontal coordinate is `pT(Z)` when `AIZ(false)` is run and `|y(Z)|` when
`AIZ(true)` is run. A factor of 1 means unchanged estimated uncertainty; above 1
means the cut worsens precision, and below 1 means the accepted events retain
enough `A4` information to give a smaller estimated uncertainty despite the
reduced sample size. This is an uncertainty ratio, not the fraction of events
accepted. Empty or undefined bins are stored as zero; they do not mean zero
uncertainty. When the zoom flag is enabled, values above 20 are clipped.

For the currently available plots:

- `AI_Z_Truth_Fiducial_CFonly_testPowheg_Y_NormXsec_A4_sensitivity_loss.pdf`
  uses the CFonly cut: both leptons have `pT >= 25 GeV`, with one in the central
  region (`|eta| < 2.5`) and one outside it (`|eta| >= 2.5`). Most populated
  bins have factors above 1. The first populated bin (`|y| = 1.0–1.2`) is about
  46.7 with the corrected formula. It contains only one accepted event and is
  clipped by the plot's y-axis maximum of 20, so treat it as a low-statistics
  result.
- `AI_Z_Truth_Fiducial_CCAsym_Zai_finalbinningPowheg_Y_NormXsec_A4_sensitivity_loss.pdf`
  uses the asymmetric central-central cut: both leptons are central, with
  leading `pT >= 27 GeV` and subleading `pT >= 25 GeV`. Its corrected factors
  rise from about 1.56 at low `|y|` to about 95.2 at `|y| = 2.4–2.6`; the
  latter bin has very low acceptance and should be treated cautiously. Bins
  above `|y| = 2.6` are empty/undefined.

These PDFs use the same `|y|` binning, but different cuts and different named MC
productions (`testPowheg` versus `Zai_finalbinningPowheg`). Their difference
cannot be attributed to the cut alone; a controlled cut comparison should use
the same underlying events and production settings. Also note the different
y-axis ranges: the CFonly plot is zoomed and clipped at 20, while the CCAsym
plot shows its full range.

### 4) One-line batch command (no interactive ROOT prompt)

```bash
conda run -n MAC root -l -b -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'test=true;' -e 'AIZ(false);' -e '.q'
```

### Notes

- `AIZ.C` now includes `TLVUtils.h` and calls `TLVUtils::getCSFAngles` and `TLVUtils::getAiPolynoms`, so load `TLVUtils.cxx` before running `AIZ`.
- **Collins-Soper azimuth convention:** `phi_CS` is stored and plotted in the ATLAS range $$[0,2\pi)$$. An equivalent signed value in $[-\pi,\pi)$ is converted by adding $2\pi$ when it is negative. This convention applies only to `phi_CS`; ordinary laboratory lepton `phi` values retain ROOT's signed convention.
- In this environment, interpreted loading is reliable. ACLiC mode with `+` (`.L TLVUtils.cxx+`) may fail due a local shared-library loading issue on macOS/conda.

## Build/Run with Makefile

The repository now includes a simple `Makefile` with these targets:

- `libs`: builds `TLVUtils.o` and `libTLVUtils.so`
- `run`: runs `AIZ(false)`
- `run-y`: runs `AIZ(true)`
- `run-test`: runs a quick test (`test=true; AIZ(false);`)
- `clean`: removes built objects/libraries

Examples:

```bash
conda run -n MAC make libs
conda run -n MAC make run-test
conda run -n MAC make run
```

## Full Workflow for New Y-Slice Leading/Subleading Plots

This is the complete workflow to produce the input ROOT file with the new histograms
and then generate the comparison PDFs, including:

- `compare_projection_overlay_pT_lead_vs_sublead_in_yZ_slices.pdf`
- `compare_projection_overlay_eta_lead_vs_sublead_in_yZ_slices.pdf`

### 1) Activate environment and enter the repo

Why: ROOT and local macros must be available from this directory.

```bash
conda activate MAC
cd /data/dust/user/ludovica/CraigMacro/Tutorial_AI/Ai_tutorial
```

### 2) Run a quick validation production (test mode)

Why: This is a fast check that `AIZ.C` writes the new histograms before running full statistics.

```bash
root -l -q \
	-e '.L TLVUtils.cxx' \
	-e '.L AIZ.C' \
	-e 'override=true;' \
	-e 'test=true;' \
	-e 'ifTrueOnly=true;' \
	-e 'FiducialCut=false;' \
	-e 'FiducialCutEtaonly=false;' \
	-e 'normXS=true;' \
	-e 'AIZ(true);'
```

Expected output file (quick check):

- `AI_Z_Truth_testPowheg_Y_NormXsec.root`

### 3) Verify the new histograms exist in the ROOT file

Why: `CompareAIZProjections.C` can only make the new overlays if these keys exist.

```bash
root -l -q -e 'TFile f("AI_Z_Truth_testPowheg_Y_NormXsec.root"); f.GetListOfKeys()->Print();' \
	| grep -E 'zY_vs_pt_leading|zY_vs_pt_subleading|zY_vs_eta_leading|zY_vs_eta_subleading'
```

You should see all 4 keys listed:

- `zY_vs_pt_leading`
- `zY_vs_pt_subleading`
- `zY_vs_eta_leading`
- `zY_vs_eta_subleading`

### 4) Run the comparison macro on the produced file

Why: This generates all projection/overlay PDFs, including the two new Y-slice leading/subleading overlays.

```bash
root -l -q 'CompareAIZProjections.C("AI_Z_Truth_testPowheg_Y_NormXsec.root", true)'
```

### 5) Confirm the two new PDF outputs are created

Why: Final validation that the workflow succeeded.

```bash
ls -1 compare_projection_overlay_*_in_yZ_slices.pdf
```

You should see at least:

- `compare_projection_overlay_pT_lead_vs_sublead_in_yZ_slices.pdf`
- `compare_projection_overlay_eta_lead_vs_sublead_in_yZ_slices.pdf`

### 6) Run full-statistics production (final file)

Why: After the quick test succeeds, regenerate your final analysis file with full event statistics.

Inside ROOT:

```cpp
.L TLVUtils.cxx
.L AIZ.C
override = true;
test = false;            // full statistics
ifTrueOnly = true;
FiducialCut = true;
FiducialCutEtaonly = false;
FiducialCutCCCF = false;
FiducialCutCFonly = true;
normXS = true;
AIZ(true);
```

Then run:

```bash
root -l -b -q 'CompareAIZProjections.C("AI_Z_Truth_Fiducial_CF_Zai_finalbinningPowheg_Y_NormXsec.root", true)'
```

### Troubleshooting

- If you do not get the new PDFs, first check Step 3 again.
- If ROOT cannot find a macro, make sure you are in:
	- `/data/dust/user/ludovica/CraigMacro/Tutorial_AI/Ai_tutorial`
- If output files already exist and are not updated, ensure `override=true;` before `AIZ(true);`.
