# FINAL MODEL — exact drive and readout channel definitions (2026-08-09)

Two-face dressing model. One Hamiltonian, all-static Fe–Tm vertices:
`W1_yy = W1_xz = 0.01` (one exchange tensor, one scale), `W3_yz = 0.01`
+ thermal (2000 / 0.05 / 0.00333), `d_z^eff` drive entry 3.0, `kappaB = 0`
(exactly inert in both geometries), `e1 = 2.067834 meV`, `e2 = 4.9628 meV`.
Pulse: the measured THz waveform (`experimental_pulse_codeunits.dat`) in
every run. Units: energies meV; f[THz] = E[meV]/4.135667696;
1 code time unit = 0.6582119569 ps.

Mode frequencies (THz): qFM 0.38, E12 0.50, E23 0.70, qAFM 0.90, E13 1.20.

## Component index convention (all datasets)
- **SU2** (Fe, per-site sum): columns 0,1,2 = S_x, S_y, S_z (crystallographic
  a, b, c). S_x = F_x (the canted moment / order-parameter direction),
  S_z = F_z (≡ 0 in geometry A by symmetry, verified to 1e-15).
- **SU3** (Tm Gell-Mann): columns 0..7 = λ1..λ8.
  λ1,λ2 = Re/Im of the 1↔2 coherence (E12); λ3 = population imbalance;
  λ4,λ5 = 1↔3 (E13); λ6,λ7 = 2↔3 (E23); λ8 = diagonal.
  Global and local frames coincide for TmFeO3 (the four transported Tm ions
  are exactly equivalent; the staggered combination is bitwise zero).
- **Product channels** (supplied explicitly — see the warning below):
  - `QNL = (λ1·λ2)_NL` — the on-site hyperpolarisability coordinate, radiates
    at 2·E12 = 1.00 THz (SHG).
  - `DNL = (F_x·λ2)_NL` — the fluctuating part of the order-induced dipole,
    radiates at qAFM ± E12 = 1.40 / 0.40 THz (SFG / DFG).

> **Why the products are stored, and cannot be derived.** A product of
> nonlinear signals is *not* a function of the stored M_NL:
> `(XY)_NL = X_01·Y_01 − X_1·Y_1 − X_0·Y_0`, which needs the two-pulse (M01),
> probe-only (M1) and pump-only (M0) trajectories separately. They are
> therefore computed at pack-build time and shipped as their own datasets.

## THE TWO OPERATORS
There are only **two** detection operators in the whole atlas, one per
analysed *polarisation state*; each serves excitation and detection alike,
so the four channels are the 2×2 map O_in → O_out.

```
O2 ≡ (E∥a, H∥c) = F_z + 5.264 λ2 (m_z, M1) + 5.264 λ1 (d_x, E1) + β (λ1λ2),  β = 66
O1 ≡ (E∥c, H∥a) = F_x + 0.9128 λ7 (m_x, M1) + 4.4 λ6 (d_z, E1)
                      + g_d(⟨F_x⟩ + δF_x) λ2   (order-induced dipole)
```
with `μ_x13 = 0` (the 1↔3 transition is magnetically dark in *every*
polarisation — measured, not assumed) and `d_z(1→3) ≤ 0.006` (bounded by the
weak 1.15 THz wing; dropped from the minimal operator).
No population generators (λ3, λ8) appear in either operator.

Channel map: **A-cross = O1→O2, A-same = O1→O1, B-cross = O2→O1, B-same = O2→O2.**

## GEOMETRY A (one experiment, two analysers) — runs/geomA_gs{1,2,3}
**DRIVE** (single pulse, E∥c / H∥a):
- Fe: Zeeman B∥a·f(t), amplitude 0.12.
- Tm (su3 vector, manual): E∥c d_z^eff → 3.0·λ2 ; B∥a μ_x → 0.9128·λ7
  (dark-μ13: λ5 component = 0). Amplitude 0.021950.

**A-cross** — analyser (E∥a, H∥c) out — is O2 with `F_z ≡ 0`:
```
S(t) = 5.264·λ2 + 5.264·λ1 + 66·QNL
```
Datasets: `linear/l2`, `linear/l1`, `products/QNL`.
Census: (qAFM,E12) = 1.00 anchor, (E12,2E12) = 0.96 SHG, (−qAFM,E12) = 0.50.

**A-same-pol** — analyser (E∥c, H∥a) out — is O1 with the propagation filter:
```
S(ω) = 4.4·λ6·T13(ω_t)·T13(ω_τ)  +  0.30·c_E12·λ2  +  0.42·(c_E12/⟨F_x⟩)·DNL·T13(ω_t)
c_E12 = 1.488e-4,   c_Fe = 0,   ⟨F_x⟩ = dataset attribute `Fx_mean`
```
Datasets: `linear/l6`, `linear/l2`, `products/DNL`, attr `Fx_mean`.
The three g_d avatars (drive 3.0; condensed emission 0.30·c_E12; dynamical
emission 0.42·c_E12/⟨F_x⟩) are **one physical constant** in three unit systems.
Census: (E12,E23) = 1.00 anchor, SFG (qAFM,1.41) = 0.31, DFG (qAFM,0.40) = 0.43,
(E12,1.41) = 0.36, (qAFM,E12) = 0.39; magnon diagonal 0.13 (sub-threshold —
its non-observation is what bounds c_Fe).

## GEOMETRY B (one experiment, two analysers) — runs/geomB_gs1
**DRIVE** (single pulse, E∥a / H∥c):
- Fe: Zeeman B∥c·f(t), amplitude 0.10.
- Tm (su3 vector, manual, `auto_su3_pump = false`):
  - B∥c μ_z → **+0.030711**·λ2 (magnetic quadrature)
  - E∥a d_x → **−0.016500**·λ1 (electric quadrature; the RELATIVE MINUS SIGN
    is required — it places the transfer row at ω_τ = 0.52 as measured)
  - `pump_direction_su3 = -0.016500, 0.030711, 0,0,0,0,0,0`, amplitude 0.034863
  - ⇒ w_E1 = 0.0165/0.0307 = 0.54 in μ units: the E1/M1 balance, fixed by the
    measured (E12,qAFM)/(qFM,qAFM) = 0.43–0.46.

**B-cross** — analyser (E∥c, H∥a) out — is O1, whose Tm part vanishes
identically: the c-directed drive touches only λ1, λ2, so the {λ1,λ2,λ3}
subalgebra closes and λ4…λ7 ≡ 0 (verified exactly). The channel is therefore
pure Fe:
```
S(t) = m_x = Fe S_x
```
Dataset: `linear/Fx`. Census: (qFM,qAFM) = 1.00, (E12,qAFM) = 0.42 transfer.

> **Cross-geometry scale caveat.** The cross-family scales (w_E1, the g_d
> avatars, c_Fe) are calibrated *within* geometry A and connect its two
> channels exactly. They do **not** transfer literally to geometry B, which
> would additionally need the absolute internal field and the
> polarisation-resolved transmission at the entry face. Empirically the
> measured B-cross map is magnon-dominated (no Tm feature above ~0.15), which
> bounds the geometry-B Tm-emission normalisation ≥5× below the naive
> geometry-A transfer. Below that bound the un-truncated moment predicts
> B-side SFG/DFG tones at (qAFM, 0.40) and (qAFM, 1.41) — a low-level
> signature worth searching for.

**B-same-pol** — analyser (E∥a, H∥c) out — is O2, the *same operator* as
A-cross, so it is a zero-parameter prediction:
```
S(t) = F_z + 5.264·λ2 + 5.264·λ1 + 66·QNL
```
Datasets: `linear/Fz`, `linear/l2`, `linear/l1`, `products/QNL`. The linear
readouts carry the directly driven coherence, which the M_NL subtraction
removes at first order, so the channel collapses onto the quadratic term:
**a single ridge at (E12, 2E12) = 1.00, everything else ≤ 0.07.**

## Thermal ensembles
A panels at T = 10 K, B panels at T = 0. The thermal map is the Boltzmann sum
over the three level-seeded trajectories:
```
M(T) = Σ_n w_n M_n ,  w_n ∝ exp(−E_n / (0.086173·T[K])) ,  E = (0, 2.067834, 4.9628) meV
```
Exact, because the classical λ dynamics *is* the von Neumann equation of the
driven damped qutrit.

## Propagation filter (E1 emission only; M1 passes freely)
```
T(ω) = exp[ −d / (1 + ((ω − ω0)/γ)²) ]
d = 6 at E13 = 1.20 THz (γ = 0.05) — applied to the emission column AND the
                                     delay-axis FID (kills the 2Q row exactly)
d = 1 at E12 = 0.50 THz (γ = 0.04)
```

## Spectra recipe (identical in every channel unless noted)
1. Detection window `t ≥ 3` code units (excludes pulse overlap).
2. Apodisation — cross channels: Hann on both axes. Same-pol: cosine gate over
   the overlap region `|τ| < 6` plus Gaussian t-apodisation to 0.03.
3. `|FFT2|`; physical delay axis `ω_τ = −fftfreq(τ)` (+ = non-rephasing).
4. Gaussian smoothing σ = (1, 0.5) bins; **linear** amplitude display.
5. The ω_τ = 0 pump–probe row is **retained** in the geometry-A panels
   (auto-scaled to the measured relative row amplitude) and **excluded** from
   every census (`|ω_τ| < 0.18 THz`).
6. Blind census: 2D local maxima, amplitude ≥ 0.05 of the map max, prominence
   ≥ 1.5 over a 0.34-THz median background, dedup radius 0.09 THz.
