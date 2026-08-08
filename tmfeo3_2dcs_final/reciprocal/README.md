# The RECIPROCAL atlas (2026-08-08): two operators, drive = detection

O2 (E||a,H||c) = S_z + 5.264λ2 + 2.84λ1 (w_E1=0.54μ) + βλ1λ2  [β recal 50.1]
O1 (E||c,H||a) = 1.65 S_x + 3.0λ2(d_z^eff) + 0.006λ4 + 4.4λ6 + 0.9128λ7, μx13=0
Drives ARE the operators: geometry B already ran O2 exactly; new geometry-A runs
(geomA_gs1-3) carry the complete O1 Tm vector (0,3.0,0,0.006,0,4.4,0.9128,0).
Propagation: E13 self-absorption depth 6 (established) + E12 depth 1 on E1 members.

VERDICT
+ Both headline cross peaks are reciprocity-proof: A-cross (0.91,0.49)=1.00 with
  harmonic 0.96/rephasing 0.50 unchanged under the completed O1 drive; B-cross
  (0.38,0.90)=1.00, (0.53,0.90)=0.45 unchanged under the O1 detection.
+ The reciprocity price in B-cross (the legally-required d_z^eff λ2 readout)
  lands WITHIN the measured residual bound once E12 self-absorption d≈1 is
  included: strays 0.36 → 0.17-0.18, and they split into a self-reversal
  DOUBLET at ω_t = 0.45/0.54 — a sharp new experimental fingerprint.
− A-same under the FULL reciprocal drive over-rotates: the map becomes
  qAFM-delay-led ((0.90,0.89)=1.00) instead of the measured E12-delay-led
  hierarchy. Interpretation: the incoming E||c field is itself attenuated at
  the CEF lines in a thick sample (same propagation physics, drive side), so
  the effective electric drive components are smaller than their bare dipoles;
  the adopted production A-drive (weaker electric writing) remains the better
  effective model until the sample's optical depths are measured.
=> Single unknown funnel: EVERY remaining discrepancy (R1 rephasing floor,
   A-same drive strength, B-cross stray magnitude, E13 wing amplitudes) now
   maps onto the sample absorption depths at E12/E23/E13 — one measurable set.
Files: atlas_reciprocal.py/.png, census_atlas_reciprocal.json, runs geomA_gs1-3.
