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

O₁ ≡ (E∥c, H∥a):
    S_x(Fe)  +  0.9128 λ⁷ (μ_x(2↔3), M1)
             +  0 λ⁵      (μ_x(1↔3): dark)
             +  4.4 λ⁶ + 0.006 λ⁴  (d_z, E1)
             +  g_d λ²    (d_z^eff: order-induced ESJ term, E-field × condensed
                           T-odd spin leg → lands legally on the T-odd quadrature;
                           ∝ ⟨F_x⟩ ⇒ collapses at the spin reorientation)
    — no λ³/λ⁸ (mirror-odd analyser has no diagonals anyway).
    Drives geometry A; detects A-same and B-cross.
    Implemented drive vector (geometry A, production = O₁ × incoming transmission):
    (0, 3.0, 0, 0, 0, 0, 0.9128, 0), amp 0.021950; strict-form (bare O₁):
    (0, 3.0, 0, 0.006, 0, 4.4, 0.9128, 0), amp 0.037822 — pathway-ratio-equivalent.

Channel map:  A-cross = O₁→O₂;  A-same = O₁→O₁;  B-cross = O₂→O₁;  B-same = O₂→O₂.

## 3. Hamiltonian operators (the κB-free vertex set)

| term | operator | value | role |
|------|----------|-------|------|
| W₁ʸʸ | S_y² λ¹            | 0.01 | two-magnon displacive vertex → geometry-A (qAFM, E₁₂) |
| W₁ˣᶻ | S_x S_z λ¹         | 0.01 | static converter: electrically written λ¹ ↔ qAFM magnon → geometry-B (E₁₂, qAFM); also the forced (qAFM,qAFM) sideband |
| W₃ʸᶻ | S_y S_z λ³ (+ thermal feedback c_heat=2000, cap 0.05) | 0.01 | population channel: reorientation physics; 2DCS tail only |
| m_z quad. | μλ² + βλ¹λ²   | β≈50–66 | hyperpolarisability (Tm 4c lacks inversion) → 2E₁₂ harmonic |
| κᴮ (BSλ), κᴱ (ESλ) | —   | 0    | catalogued; superseded / rejected (κE_5y floods the 1↔2 block) |
| other quadratics | λ¹λ⁶, λ²λ⁷ (d_z E1); λ¹λ⁷, λ²λ⁶ (m_x M1) | 0 | all four symmetry-allowed; radiate at E₁₂+E₂₃=1.20; tested for the (qAFM,1.2) peak — each carries a 10–12× (E₁₂,1.2) partner (opposite of the measured ratio) → rejected; R2 stays open |
| Fe sector | J₁=4.74/5.15, J₂=0.15/0.30, Ka=−0.0153, Kc=−0.0187, D₁=0.049, D₂=0.014 | | magnons; (qFM,qAFM) magnon–magnon peak |
| Tm onsite | −(Δ₁₂/2)λ³ + ((Δ₁₂−2Δ₁₃)/2√3)λ⁸, e₁=2.0678, e₂=4.9628 | | qutrit spectrum |

## 4. Propagation + calibration layer (shared by drive and detection)

- Transmission filter, E1 components only (M1 transparent):
  E₁₂ line depth ≈ 1 (E∥c d_z^eff emission), E₁₃ line depth ≈ 6 (both the
  emission column and the delay-row FID; self-reversal wings at 1.15/1.25).
- Three cross-family scales (one number each, used consistently):
  w_E1 = 0.54 (E1:M1 — measured by the B transfer ratio);
  g_d (d_z^eff = the ESJ vertex E_z·S·λ² with the spin leg condensed on
  ⟨F_x⟩: an order-induced electric dipole of the 1↔2 transition, g_d ∝ ⟨F_x⟩,
  hence collapses at the spin reorientation. ONE parameter seen through two
  unit systems: drive avatar = the su3 vector entry 3.0 (g_d × internal field,
  calibrated by the E12-excitation amplitude); detection avatar =
  c_E12 = 1.488e-4 per raw λ² map (g_d × emission conversion, calibrated by
  the measured (qAFM,E12) = 0.30). Connecting the two would need the absolute
  field inside the sample — reciprocity holds but is untestable on this term
  by ratios alone); Fe:Tm emission scale (the g-knob,
  c_Fe ≤ 2.0e−5 per raw S_x map vs 4.4 per raw λ⁶ map — detector/mode overlap;
  bounded from ABOVE by the experimental non-observation of the (qAFM,qAFM)
  magnon diagonal in A same-pol (atlas uses 2.0e−5 → diagonal 0.14, below the
  digitisation threshold; the λ⁶-carried (E₁₂,qAFM) transfer is insensitive,
  0.46→0.43).

Strictly reciprocal: operator content, within-family dipole ratios, and the
filter. Calibrated: the three cross-family scales. Everything else is zero.
