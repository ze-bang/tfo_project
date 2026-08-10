# The complete operator inventory (final κB-free scenario, 2026-08-08)

## 1. Qutrit coordinate dictionary (λ¹..λ⁸)

| λ  | content            | freq (THz) | T-parity | mirror m_z | multipole & owner                                   |
|----|--------------------|-----------|----------|------------|-----------------------------------------------------|
| λ¹ | Re ρ₁₂             | 0.50      | even     | even       | electric d_x — owned by E∥a; weight w_E1·μ = 2.84   |
| λ² | Im ρ₁₂             | 0.50      | odd      | even       | magnetic μ_z = 5.264 — owned by H∥c; PLUS the order-induced electric d_z^eff (E∥c, ∝⟨F_x⟩, strength g_d) |
| λ³ | pop. n₁−n₂         | 0 (rect.) | even     | even       | d_x diagonal symmetry-allowed but EXCLUDED from all detection (permanent moments don't radiate); dynamical coordinate only |
| λ⁴ | Re ρ₁₃             | 1.20      | even     | odd        | electric d_z(1↔3) = 0.006 — owned by E∥c (nearly dark) |
| λ⁵ | Im ρ₁₃             | 1.20      | odd      | odd        | magnetic μ_x(1↔3) = 0 — DARK (CEF value 2.39 rejected by data: the fully dark 1↔3 rule) |
| λ⁶ | Re ρ₂₃             | 0.70      | even     | odd        | electric d_z(2↔3) = 4.4 — owned by E∥c              |
| λ⁷ | Im ρ₂₃             | 0.70      | odd      | odd        | magnetic μ_x(2↔3) = 0.9128 — owned by H∥a           |
| λ⁸ | diagonal           | 0 (rect.) | even     | even       | excluded from detection (as λ³)                     |

Fe (per-site sum): S_x = F_x (M1, radiates into H∥a analyser; F_x ≡ order
parameter direction), S_y = the pumped magnon coordinate (radiates m_y —
seen by NEITHER analyser), S_z = F_z (M1, H∥c analyser; ≡ 0 in geometry A).

## 2. The two polarisation operators (drive content = detection content)

O₂ ≡ (E∥a, H∥c):
    S_z(Fe)  +  5.264 λ²  (μ_z, M1)
             +  2.84 λ¹   (d_x, E1; w_E1 = 0.54 μ, MEASURED; sign(d_x μ_z) < 0)
             +  β λ¹λ²    (hyperpolarisability, emission only; β from A-cross parity ≈ 50–66)
    — no λ³/λ⁸. Drives geometry B; detects A-cross and B-same.
    Implemented drive vector (geometry B): (−0.0165, 0.030711, 0,…), amp 0.034864.

O₁ ≡ (E∥c, H∥a)  [MINIMAL form, necessity audit 2026-08-08]:
    S_x(Fe)  +  0.9128 λ⁷ (μ_x(2↔3), M1)
             +  0 λ⁵      (μ_x(1↔3): dark)
             +  4.4 λ⁶    (d_z(2↔3), E1)   [λ⁴ DROPPED: d_z(1↔3) ≤ 0.006
                           kept as a bound from the 1.15 wing — audit: marginal]
             +  g_d (⟨F_x⟩ + δF_x) λ²
                          (d_z^eff UN-TRUNCATED, adopted 2026-08-08: the ESJ term
                           E_z·(F_x λ²).  Condensed piece ⟨F_x⟩λ² = the old static
                           avatars; dynamical piece δF_x·λ² is a PRODUCT emitter
                           radiating at qAFM±E12 = 1.40/0.40 with qAFM delay →
                           the measured (0.90,1.41)=0.30 sum tone, the 0.40
                           difference tone, and the (E12,1.41) partner.  The ONLY
                           channel in the ~40-entry catalogue whose map is led by
                           the qAFM-delay row (both factors first-order-written).
                           Whole family ∝ ⟨F_x⟩ ⇒ collapses at the reorientation)
    — no λ³/λ⁸ (mirror-odd analyser has no diagonals anyway).
    Drives geometry A; detects A-same and B-cross.
    Implemented drive vector (geometry A, production = O₁ × incoming transmission):
    (0, 3.0, 0, 0, 0, 0, 0.9128, 0), amp 0.021950; strict-form (bare O₁):
    (0, 3.0, 0, 0.006, 0, 4.4, 0.9128, 0), amp 0.037822 — pathway-ratio-equivalent.

Channel map:  A-cross = O₁→O₂;  A-same = O₁→O₁;  B-cross = O₂→O₁;  B-same = O₂→O₂.

## 3. Hamiltonian operators (the κB-free vertex set)

| term | operator | value | role |
|------|----------|-------|------|
| W (ONE tensor, ONE scale) | S^T W S λ¹: components W₁ʸʸ = W₁ˣᶻ = 0.01 | 0.01 | the exchange-modulated crystal field. NOT two mechanisms (audit 2026-08-08): the two Γ₂-allowed components share every conversion task — A-cross anchor (yy-led), forced (qAFM,qAFM) sideband (yy-dominant 0.43/0.64 split), B transfer (multi-path; xz removal costs 25%) |
| W₃ʸᶻ | S_y S_z λ³ (+ thermal feedback c_heat=2000, cap 0.05) | 0.01 | ATLAS-INERT (audit: all coherent features ≤5% moved without it). Enters ONLY the incoherent observables: rectified (±E₁₂,0) line magnitude and the buildup reservoir |
| m_z quad. | μλ² + βλ¹λ²   | β≈50–66 | hyperpolarisability (Tm 4c lacks inversion) → 2E₁₂ harmonic |
| κᴮ (BSλ), κᴱ (ESλ) | —   | 0    | catalogued; superseded / rejected (κE_5y floods the 1↔2 block) |
| other quadratics | λ¹λ⁶, λ²λ⁷ (d_z E1); λ¹λ⁷, λ²λ⁶ (m_x M1) | 0 | all four symmetry-allowed; radiate at E₁₂+E₂₃=1.20; tested for the (qAFM,>1.2) peak — each carries a 10–12× (E₁₂,1.2) partner (opposite of the measured ratio) → rejected. The peak is instead the dynamical d_z^eff sum tone at qAFM+E12=1.40 (see O₁) — R2 RESOLVED 2026-08-08 |
| Fe–Fe striction | d ∝ G_xG_y / F_yG_y pair dipoles | 0 | right position (qAFM+qFM=1.28) but full-strength unobserved combination delay rows at \|ω_τ\|=1.05/1.49 → rejected |
| Fe sector | J₁=4.74/5.15, J₂=0.15/0.30, Ka=−0.0153, Kc=−0.0187, D₁=0.049, D₂=0.014 | | magnons; (qFM,qAFM) magnon–magnon peak |
| Tm onsite | −(Δ₁₂/2)λ³ + ((Δ₁₂−2Δ₁₃)/2√3)λ⁸, e₁=2.0678, e₂=4.9628 | | qutrit spectrum |

## 4. Propagation + calibration layer (shared by drive and detection)

- Transmission filter, E1 components only (M1 transparent):
  E₁₂ line depth ≈ 1 (E∥c d_z^eff emission), E₁₃ line depth ≈ 6 (both the
  emission column and the delay-row FID; self-reversal wings at 1.15/1.25).
- Three cross-family scales (one number each, used consistently):
  w_E1 = 0.54 (E1:M1 — measured by the B transfer ratio);
  g_d (d_z^eff = the ESJ vertex E_z·(F_x λ²), one parameter seen through THREE
  avatars: drive avatar = su3 vector entry 3.0 (condensed, calibrated by the
  E12-excitation amplitude); static emission avatar = 0.30 × c_E12 with
  c_E12 = 1.488e-4 per raw λ² map (condensed, rebalanced 2026-08-08 when the
  dynamical piece was adopted; the (qAFM,E12) blob is now fed jointly);
  dynamical emission avatar = s·c_E12/⟨F_x⟩ on the (F_x·λ²) product map with
  s = 0.42 (condensation-ratio correction, calibrated by the measured
  (qAFM,1.41) = 0.30 sum tone; the naive trajectory value δF_x/⟨F_x⟩,
  ⟨F_x⟩ = 0.0084, overshoots ×2 — the model's canting-angle uncertainty));
  Fe:Tm emission scale (the g-knob: the atlas now uses c_Fe = 0;
  bounded from ABOVE by the experimental non-observation of the (qAFM,qAFM)
  magnon diagonal in A same-pol (diagonal 0.13 at c_Fe=0, below the
  digitisation threshold; the λ⁶-carried (E₁₂,qAFM) transfer is insensitive).

Strictly reciprocal: operator content, within-family dipole ratios, and the
filter. Calibrated: the three cross-family scales. Everything else is zero.

## Cross-geometry scale caveat (2026-08-09, user consistency check)
The cross-family scales (w_E1, g_d avatars, c_Fe) are calibrated WITHIN
geometry A and connect its two channels exactly. They do NOT transfer
literally to geometry B: a strict full-O1 B-cross composite with the
A-calibrated numbers becomes Tm-dominated (λ²/SFG-DFG at 0.6–1.0),
whereas the measured B-cross is magnon-dominated (no Tm feature > ~0.15).
=> per-geometry absolute emission normalisation (internal field +
polarisation-resolved entry-face transmission) is a separate constant;
the measured magnon dominance of B-cross is a NULL-BOUND placing the
B-side Tm-emission scale ≥5× below the naive A-transfer. Prediction at
whatever scale the bound allows: B-cross SFG/DFG tones at (0.90, 0.40)
and (0.90, 1.41). The atlas B-cross panel displays the Fe part (the
measured dominance); operator CONTENT is O1 everywhere.
