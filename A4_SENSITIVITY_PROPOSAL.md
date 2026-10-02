# A4 Sensitivity Loss from Fiducial Cuts — Proposal

Scope: `A4` (the Collins–Soper `P4` angular coefficient) only. No other `A_i` are affected.

## 1. Sensitivity-metric definition

`A4` is extracted as a weighted moment:

```
A4 = 4 * <P4(cosθ_CS, φ_CS)>   (per kinematic bin, weight = mcEventWeight)
```

The statistical precision of this moment estimator is driven by two
independent ingredients, both computable from quantities already produced by
the macro:

1. **Effective sample size** (Kish effective N), which accounts for MC event
   weights:

   ```
   N_eff = (Σw)^2 / Σ(w^2)
   ```

   `Σw` and `Σ(w^2)` are exactly the quantities already tracked internally by
   the `AIZP6::Moment` helper (`sumW` and the bin error of `sumW`, since
   `sumW->Sumw2()` stores `Σ(w^2)` as the squared bin error).

2. **Lever arm** (how spread the `P4` values are among surviving events):

   ```
   <P4^2> = Σ(w * P4^2) / Σw
   ```

   computed per bin the same way as `<P4>`, just filling the square of the
   polynomial value instead of the value itself.

For this sensitivity proxy, the information on the A4 coefficient is taken to
scale with the effective event count times the mean squared P4 lever arm:

```
I(A4) ∝ N_eff * <P4^2>
σ(A4) ∝ 1 / sqrt(I(A4))
```

**Sensitivity loss factor**, comparing the cut ("after") to the uncut
("before") sample in the same kinematic bin:

```
SensitivityLossFactor = σ(A4)_after / σ(A4)_before
                       = sqrt( (N_eff,before * <P4^2>_before)
                             / (N_eff,after  * <P4^2>_after ) )
```

Interpretation:
- `= 1`: the fiducial cut removes events without degrading A4 precision
  (e.g. the retained sample's information compensates for the reduced event
  count).
- `> 1`: the cut inflates the statistical uncertainty on A4 beyond what is
  explained by lost acceptance alone (remaining events carry less
  angular information on average).
- `< 1`: the surviving sample is, bin-by-bin, more informative per event
  enough to compensate for the event loss.

This is a single-moment proxy (no correlations with other `A_i`, no
background), chosen because it is transparent, cheap, and reuses existing
accumulators (`AIZP6::Moment`) rather than requiring a full fit.

## 2. Output format

### 2.1 Histogram (written to the output ROOT file)

- `A4_sensitivity_loss`: 1D histogram, same binning as `A4` (i.e. respecting
  the `isY` switch: `|y(Z)|` bins or `p_T(Z)` bins, using the existing
  `bins[]`/`Nbins` arrays).
  - Bin content: `SensitivityLossFactor` per kinematic bin.
  - Bin error: propagated from the statistical uncertainties of `N_eff` and
    `<P4^2>` in each bin (standard error propagation through the ratio/sqrt).
  - Axis titles follow the existing `A4` convention
    (`;|y(Z)| or p_{T}(Z) [GeV];#sigma(A_{4})_{after} / #sigma(A_{4})_{before}`).

- Supporting (non-essential, diagnostic) histograms, also written, mirroring
  the `before`/`after` naming convention already used for `eff_*` histograms:
  - `A4_Neff_before`, `A4_Neff_after`
  - `A4_P4sq_before`, `A4_P4sq_after` (i.e. `<P4^2>` per bin)

### 2.2 Plot (PDF)

- One canvas, e.g. `<outputPrefix>_A4_sensitivity_loss.pdf`, drawing
  `A4_sensitivity_loss` with error bars (`"E1"` draw option), a horizontal
  reference line at `y = 1`, following the same `AIZP6::Style` conventions
  used elsewhere in the macro (centered axis titles, no stats box).
- Optional second panel (or separate PDF) showing `N_eff_after/N_eff_before`
  and `<P4^2>_after / <P4^2>_before` separately, so the two contributing
  factors to the loss can be visually disentangled.

### 2.3 Printout (stdout)

A per-bin summary table, in the same style as existing efficiency printouts,
e.g.:

```
======================================================================
A4 SENSITIVITY LOSS FROM FIDUCIAL CUT
======================================================================
 Bin low   Bin high   N_eff(before)   N_eff(after)   <P4^2>(before)  <P4^2>(after)  LossFactor
  ...
----------------------------------------------------------------------
Inclusive (full range): LossFactor = X.XXX
======================================================================
```

plus a single inclusive (non-differential) number computed the same way
over the full kinematic range, for a quick headline comparison.


The plotted quantity is **uncertainty after / uncertainty before**, not
sensitivity itself. Values above 1 indicate worse precision after the cut;
values below 1 indicate that the retained sample has greater estimated
information despite having fewer events. It is not an event-acceptance ratio.

For the default central-central cut, both leptons must have
`pT >= 25 GeV` and `|eta| <= 2.5`. The factor uses the product
`N_eff * <P4^2>` as a simple information proxy. Both the event-count change and
the change in the P4 lever arm matter: a decrease in either one reduces the
proxy information and increases the uncertainty estimate. This is a
single-moment approximation, not a replacement for a full correlated fit.


Ran terminal command: pdftotext -layout AI_Z_Truth_Fiducial_CCAsym_testPowheg_Y_NormXsec_A4_sensitivity_loss.pdf - 2>&1

Ran terminal command: grep -n "A4_SensitivityLossFactor\|GetYaxis\|GetXaxis\|SetTitle" AIZ.C | grep -i "sensitiv" 

Read , lines 1155 to 1168
Completed with input: {
  "path": "/data/dust/user/ludovica/CraigMacro/Tutorial_AI/Ai_tutorial/AIZ.C",
  "view_range": [
    1155,
    1168
  ]
}

Ran terminal command: grep -n "Double_t bins\[\|Int_t Nbins\|double bins\[\|bins\[\] =" AIZ.C | head -10

Ran terminal command: sed -n '1100,1162p' AIZ.C

Ran terminal command: grep -n "Nbins\s*=\|Int_t Nbins\|const int Nbins\|Double_t bins\[Nbins\]\|TH1D \*A4 =\|TH1D\* A4 =" AIZ.C | head -20

Ran terminal command: sed -n '388,414p' AIZ.C

### PLOT explanation ::


**Y‑axis units:** The y-axis is **dimensionless** — it's a ratio of two uncertainties, `σ(A₄)_after / σ(A₄)_before`, so it has no physical unit (just like any ratio of two quantities with the same units). A value of 1 means no change in precision; above 1 means the cut increases the uncertainty on A₄ (loss of sensitivity); below 1 means the cut decreases the uncertainty on A₄ (gain in sensitivity).

Here's a short explanatory paragraph you can drop into a note/slide for this specific plot:

> **Figure: `AI_Z_Truth_Fiducial_CCAsym_testPowheg_Y_NormXsec_A4_sensitivity_loss.pdf`**
> This plot shows the **A₄ sensitivity-loss factor** as a function of `|y(Z)|` (x-axis, "bin coordinate", in bins of width 0.2 from 0 to ~4.6), for the asymmetric central-central (CCAsym) fiducial cut applied to truth-level Powheg Z→ee events. The y-axis, `σ(A₄)_after/σ(A₄)_before`, is a **dimensionless ratio** of the estimated statistical uncertainty on the A₄ angular coefficient *after* the fiducial cut to the uncertainty *before* the cut, computed per `|y(Z)|` bin as
> `sqrt[(N_eff,before·⟨P₄²⟩_before) / (N_eff,after·⟨P₄²⟩_after)]`,
> where `N_eff` is the effective (Kish) event count and `⟨P₄²⟩` is the mean squared P₄ Legendre-moment weight ("lever arm") in that bin. A horizontal dashed line at 1 marks "no change." Values above 1 (most bins here, rising steeply above `|y(Z)|≈1.4`) indicate the CCAsym cut degrades A₄ sensitivity in that bin, increasingly so at high rapidity where the forward/central acceptance becomes very restrictive and few events survive (e.g. the last populated bin reaches ~96, i.e. a ~96× larger statistical uncertainty than before the cut).

Let me know if you'd like this trimmed down further, or tailored for a specific audience (internal note vs. slide caption).