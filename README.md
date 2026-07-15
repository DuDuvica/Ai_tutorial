# Ai_tutorial

Quick plotting utilities for `Ai`.

## Quick Start (3 commands)

```bash
conda activate MAC && cd /data/dust/user/ludovica/CraigMacro/Tutorial_AI/Ai_tutorial
root -l -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'override=true;' -e 'test=true;' -e 'ifTrueOnly=true;' -e 'FiducialCut=false;' -e 'FiducialCutEtaonly=false;' -e 'normXS=true;' -e 'AIZ(true);'
root -l -q 'CompareAIZProjections.C("AI_Z_Truth_testPowheg_Y_NormXsec.root", true)'
```

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

### 4) One-line batch command (no interactive ROOT prompt)

```bash
conda run -n MAC root -l -b -q -e '.L TLVUtils.cxx' -e '.L AIZ.C' -e 'test=true;' -e 'AIZ(false);' -e '.q'
```

### Notes

- `AIZ.C` now includes `TLVUtils.h` and calls `TLVUtils::getCSFAngles` and `TLVUtils::getAiPolynoms`, so load `TLVUtils.cxx` before running `AIZ`.
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
FiducialCut = false;
FiducialCutEtaonly = false;
FiducialCutCF = true;
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
