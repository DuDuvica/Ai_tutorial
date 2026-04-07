# Ai_tutorial

Quick plotting utilities for `Ai`.

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
