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

## FINAL reciprocity accounting (2026-08-08, after the drive-attenuation scan)
1. DRIVE-COMPLETION NULL: scanning the lambda6 electric drive over t23 in
   {1.0, 0.5, 0.2, ~0} changes the A-same pathway RATIOS negligibly (variants
   near-identical under fixed detection). The earlier "over-rotation" was a
   detection-UNIT error (drive-vector entries used as map weights), not drive
   physics. Strict reciprocity of the drive is therefore FREE: the complete O1
   vector and the production vector give the same physics up to the anchor
   normalization (direct lambda6 writing inflates the (E12,E23) anchor ~2x,
   deflating the magnon-emission family relatively; production drive matches
   the data best and equals the reciprocal drive filtered by the sample's
   incoming-field transmission at the strongly absorbed E1 lines).
2. WHAT IS STRICT: operator CONTENT per polarisation (selection rules, all
   verified); within-family dipole ratios shared by drive and detection;
   the propagation filter (E12 d~1, E13 d~6) applied to fields both ways.
3. WHAT IS CALIBRATED (one number each, used consistently): Fe:Tm emission
   scale (the g-knob, ~300x vs bare moments — mode overlap/detector), the
   d_z^eff strength (enters drive and detection as the same unknown g_d;
   ratios cannot overdetermine it without absolute field calibration), and
   w_E1 = 0.54 (E1:M1, MEASURED by the geometry-B transfer).
4. The canonical experiment-facing atlas remains: production A-drive
   (= transmission-filtered O1) + adopted A-same composite (asame_final.png);
   the strict-form atlas (bare O1 drive, fitted scales) is atlas_reciprocal.png
   -- identical headline peaks, anchor-rescaled magnon family (x~0.5).
