# AIZ vs TLVUtils: Differences and Changes to Apply

## Scope
Comparison of Collins-Soper angle and A0-A7 implementations between:
- `AIZ.C`
- `TLVUtils.cxx` (`getCSFAngles`, `getAiPolynoms`)

## Findings (differences)

### 1) `cos(theta_CS)` numerator sign is opposite

**AIZ (`AIZ.C`, lines ~397-401):**
```cpp
double lpe  = em.E()+em.Pz();
double lme  = em.E()-em.Pz();
double lpeplus = ep.E()+ep.Pz();
double lmeplus = ep.E()-ep.Pz();

double costheta = (1./z.M())*pow((z.M2()+pow(z.Pt(),2)),-0.5)*(lpeplus*lme - lmeplus*lpe);
costheta *= z.Rapidity()>=0. ? 1:-1;
```

**TLVUtils (`TLVUtils.cxx`, lines ~158-166):**
```cpp
costh = (Lplus * Pminus - Lminus * Pplus);
costh *= TMath::Abs(boson.Pz());
costh /= (boson.Mag() * boson.Pz());
costh /= TMath::Sqrt(boson.Mag2() + boson.Pt() * boson.Pt());
```

With `L` from negative lepton and `P` from positive lepton, AIZ uses
$$P^+L^- - P^-L^+ = -(L^+P^- - L^-P^+)$$
so the numerator sign is flipped relative to TLVUtils.

Impact: AIZ `costheta` may be globally sign-flipped relative to TLVUtils.

---

### 2) Sign application method differs (rapidity sign vs `Abs(Pz)/Pz`)

AIZ multiplies by `sign(y_Z)`.
TLVUtils multiplies by `Abs(Pz)/Pz`, i.e. `sign(Pz)`.

These are usually equivalent for massive bosons, but TLVUtils has a potential singular behavior when `Pz -> 0` because of division by `boson.Pz()`.

Impact: potential instability at/near `Pz=0` and non-identical edge behavior.

---

### 3) `phi_CS` is computed with a different construction

**AIZ** uses a transverse-plane construction with:
- `rt = pa.Cross(pz)` and then `rt.SetZ(0)`
- `qttrans = z.Vect(); qttrans.SetZ(0)`
- `phi = atan2(y, x)` with rapidity-sign factor in `y`

**TLVUtils** computes `phi` in the boson rest frame using boosted beam axes (`xAxis`, `yAxis`) and boosted negative lepton.

Impact: not guaranteed to be numerically identical event-by-event. Can induce sign/phase convention differences.

---

### 4) Event handling difference for `pT(Z)=0`

AIZ skips events when `qttrans.Mag()==0 || rt.Mag()==0`.
TLVUtils does not impose this event-level skip in `getCSFAngles`/`getAiPolynoms`.

Impact: AIZ drops events entirely (including contributions to A0/A4) where TLVUtils-style workflow may still keep them.

---

### 5) Angular polynomial basis differs in convention usage

AIZ fills **scaled** basis directly:
- A0: `20/3*(0.5-1.5 c^2)+2/3`
- A1: `5*(2 c s cosphi)`
- A2: `10*(s^2 cos2phi)`
- A3: `4*(s cosphi)`
- A4: `4*c`
- A5: `5*(s^2 sin2phi)`
- A6: `4*(sin2theta sinphi)`
- A7: `4*(s sinphi)`

TLVUtils `getAiPolynoms` returns **raw** polynomials without those prefactors.

Impact: histogram contents are not directly comparable unless converted to the same convention.

---

### 6) Safety differences

- TLVUtils clamps `costh` to `[-1, 1]`.
- AIZ does not clamp before `sqrt(1-costheta*costheta)`.

Impact: rare floating precision overshoot can make `sinth` NaN in AIZ.

## Changes to Apply (to make AIZ match TLVUtils)

### Minimal correctness alignment

1. Replace AIZ `costheta` numerator with TLVUtils sign convention:
```cpp
// current AIZ numerator term (opposite sign)
(lpeplus*lme - lmeplus*lpe)

// TLVUtils-compatible numerator sign
(lpe*lmeplus - lme*lpeplus)
```

2. Use TLVUtils-style sign factor and avoid singularity:
```cpp
const double pzBos = z.Pz();
if (fabs(pzBos) < 1e-12) continue; // or define a stable fallback
costheta *= fabs(pzBos) / pzBos;
```

3. Clamp `costheta` before trig:
```cpp
costheta = std::max(-1.0, std::min(1.0, costheta));
```

4. Remove unused variable:
```cpp
// remove: double tanphi = ...
```

### Full definition alignment (recommended)

5. Compute `costheta` and `phi` through `TLVUtils::getCSFAngles(...)` in AIZ, using:
- `lep1 = em` (negative lepton)
- `charge1 = -1`
- `lep2 = ep`
- same `ebeam` used in analysis

6. Compute raw polynomials through `TLVUtils::getAiPolynoms(costheta, phi, aipols)`.

7. Decide and document one convention for stored histograms:
- Option A (TLVUtils raw): store `aipols[i]` directly.
- Option B (current AIZ scaled): apply the scaling factors after `getAiPolynoms` and keep existing plot conventions.

8. Align `pT(Z)=0` behavior with TLVUtils policy (do not skip whole event unless explicitly intended).

## Practical patch order

1. Fix `costheta` sign and clamp.
2. Replace local `phi` construction with `TLVUtils::getCSFAngles`.
3. Replace manual polynomial block with `TLVUtils::getAiPolynoms`.
4. Keep or remove scaling factors depending on chosen output convention.
5. Validate with event-by-event comparison plots for `(costheta_AIZ - costheta_TLV)` and wrapped `(phi_AIZ - phi_TLV)`.

## Expected result after alignment

After steps 1-3, AIZ and TLVUtils should match at the level of angular definitions. Remaining differences should then only come from:
- chosen polynomial normalization convention (raw vs scaled), and
- event selection policy (especially handling of `pT(Z)=0` / `Pz(Z)=0`).
