# A-same-pol iteration campaign (2026-08-08) — ledger and verdict

Target: the digitised experimental map (paper/data/experiment_geomA_samepol_digitized.json).
Analyzer: canonical pipeline validated against the trusted baseline table; same-pol
convention REQUIRES the tau-gate (|tau|<6 cosine ramp) + samepol t-apodization.

## Probes run (all geomA gs1, tau step 0.2)
- prA_d6p/d6m: lambda6 electric drive quadrature, both signs -> NULL (pathway
  ratios are structure-constant-fixed, quadrature-blind; unlike geometry B there
  is no interference lever).
- prA_l7m: lambda7 drive sign flip -> exact NULL.
- prA_kE5y/kE3/kE5: kappaE_5y (ESJ, the only symmetry-allowed magnon->E13 route,
  never tested before) = 0.003/0.005/0.01 -> real (qAFM,E13) lever BUT floods the
  12-block: at 0.003 already (qAFM,E12) 0.32->0.57 while the reabsorbed marked
  peak stays 0.05. NOT ADOPTED.

## What the iteration ESTABLISHED (adopted)
1. tau-gate is mandatory for this channel (documented convention; halves 2Q).
2. E13 SELF-ABSORPTION (f13=8.7, the largest oscillator strength): Lorentzian
   optical depth ~6 on BOTH axes (emission column + delay-row FID).
   - kills the (E13,E23) two-quantum phantom EXACTLY (0.33 -> 0.00)
   - explains the experimental E13 feature positions 1.15/1.25 = self-reversal
     WINGS (+-1 absorption linewidth around 1.20) — positional evidence.
3. The rectified (±E12, 0) pair (strongest measured feature): symmetric ±E12
   tau-label reproduced in-model but orders too weak under coherent pumping ->
   direct optical evidence of the INCOHERENT population channel (same reservoir
   as the geometry-B post-pulse buildup); calibration target for c_heat*W3.

## FINAL A-same scorecard (baseline + gate + self-absorption, published weights)
  (E12,E23)  exp 1.00  model 1.00     (E12,qFM)  0.45 / 0.44
  (E12,qAFM) 0.40 / ~0.45             (qAFM,E12) 0.30 / 0.32
  (E13,E23)2Q  0 / 0.00               (E12,1.15 wing) 0.15 / ~0.07
  three NAMED residuals:
  R1 (-E12,E23) model ~0.56 vs exp <0.15: exact-qutrit point-ion floor
     (0.47-0.58, verified vs standalone integrator). Extrinsic candidate:
     E12-line propagation asymmetry; needs sample thickness/absorption data.
  R2 (qAFM,1.25) marked peak 0.30: NO in-model route survives all constraints
     (kappaE floods 12-block; E1-emitted version self-absorbed; M1-transparent
     version re-injects the (E12,E13) cascade 8:1). THE open item, now sharply
     characterized: requires the 1->2->3 cascade to be ~15x weaker than ideal
     three-level dynamics, i.e. physics of the real level 3 beyond the qutrit.
  R3 (qAFM,qAFM) diagonal model 0.25 vs exp unlisted (digitisation threshold
     ~0.15) — borderline, not necessarily a conflict.

## What experiment can settle
- sample thickness / 1D absorption depth at E13 and E12 -> pins d13 and tests R1;
- T-dependence of the (qAFM,1.25) peak (CEF-population-fed vs magnon);
- a lineout across omega_t = 1.1-1.3 at omega_tau = qAFM: self-reversal predicts
  a DOUBLET (wings) around a dark center at 1.20.
Figures: asame_gap.png (gap quantification), asame_final.png (final comparison).

## 2026-08-08: the (0.90, >1.2) hunt — dynamical completion of d_z^eff

Goal: user wants the measured qAFM-row peak above 1.2 THz (digitized at
(0.90,1.25)=0.30) IN the A same-pol model, after the R2 exhaustion.

Sweep: 33 product emission channels (all Fe-Fe F/G bilinears, all Fe-Tm
F/G x lambda1/lambda2 cross moments) evaluated from the production
trajectories (M01*M01 - M1*M1 - ref*ref, geomA_gs1-3, T=10K mix).

Verdict table (own-max normalized target vs (0.50,col) partner):
- Fe-Fe striction (GxGy, FyGy): target at 1.28 (right position!), ratio 2.0,
  BUT unobserved delay rows at |w_tau|=1.05/1.49 at FULL strength -> rejected.
- Fe-Tm cross moments Fxl2/Fxl1/Gyl2/Gyl1: target (0.90,1.39)=1.00 = the
  map's TOP feature (first channel ever with the qAFM row dominant);
  partner ratio 1.4-1.6.
- Fxl2 is NOT a new coupling: d_z^eff = ESJ E_z(F_x lambda2) with the spin
  leg condensed on <F_x>. The full operator is g_d(<F_x>+dF_x)lambda2;
  we had truncated the fluctuation. dF_x*lambda2 radiates at
  qAFM±E12 = 1.40/0.40 with qAFM delay. E1 -> T13 column filter (mild, 0.68
  at 1.39).
- Zero-parameter check (weight = dF_x/<F_x> from trajectories, <F_x>=0.0084):
  OVERSHOOTS x2: (0.90,1.39)=0.66, (0.90,0.40)=0.95, diag 0.35.
- Compromise fit (s=0.42 on dyn, static lambda2 avatar x0.30, cFe=0):
  (0.90,1.39)=0.29 [target 0.30] | (0.90,0.50)=0.38 | (0.90,0.40)=0.42
  (companion, blends with 0.50 into one broad blob) | diag=0.13 (stays
  sub-threshold) | (0.50,1.39)=0.41 (LIABILITY: unobserved partner; base
  model already had 0.22 there) | anchors all unchanged.
Figure: asame_dynamical_dzeff.png (current vs +dyn side by side).

DECISION POINTS (experimental):
1. POSITION: this mechanism puts the peak at qAFM+E12 = 1.39-1.44, NOT at
   1.25. If the measured peak is truly at 1.25 (E13 upper self-reversal
   wing), the candidate fails on position and R2 stands (real level-3).
2. Is there anything measured at (E12, ~1.4)? Model predicts 0.41 there.
3. Does the measured (0.90, 0.5) blob extend down to 0.40 (the difference
   tone qAFM-E12)? Model predicts a merged 0.40+0.50 double structure.
Falsifiables if adopted: whole family ~ g_d -> collapses at spin
reorientation together with (qAFM,E12); peak tracks qAFM+E12 in T; CEF
linewidth; difference partner at 0.40 mandatory.
NOT yet wired into the atlas/tex - awaiting position confirmation.
