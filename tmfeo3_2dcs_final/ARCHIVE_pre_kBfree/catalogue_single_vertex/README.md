# Single-vertex 2DCS catalogue (one interaction at a time, small sweeps, both geometries)

Baseline (all Fe-Tm couplings OFF, drives ON): Fe J/K/D + CEF (e1,e2) + gamma(l4-7)=0.15,
pump amp 0.1 / width 0.5, g_tm=0.058333; A: pump 1,0,0 + su3 dir (0,d,0,0,2.3915,0,0.9128,0), d=3.0
unless swept; B: pump 0,0,1 + auto mu-projection. Catalogue resolution: tau [-60,0] step 0.2,
t [-70,40] step 0.02, corrected recipe (honest M1, no chunking).

Runs: {A,B}_baseline; {A,B}_W1yy_{0.01,0.03,0.05}; {A,B}_W1xz_{0.01,0.03,0.09};
{A,B}_kB1y_{0.005,0.01,0.03}; {A,B}_Km2y_{0.001,0.002,0.005}; {A,B}_kE5x_{0.02,0.1};
A_d{0.0,1.6} (drive sweep; d=3.0 is baseline).

Per run (h5 deleted after distillation):
- census.json : all local maxima >=7% per channel (FeSx, l2, l4, l5, l6, l7), labeled
- spectra.npz : windowed |FFT| maps (float32; keys=channels, wt, wta axes)
- panels.png  : 6-channel overview (each channel self-normalized)
- reproduce with the .param recipe above (exact params in ../../..//archive params or regenerate
  via scratchpad catalogue_gen.py)

Quick findings visible in the censuses:
- A_kB1y_* == A_baseline exactly: the B_z gate in action (kappaB silent for H||a).
- kE5x inert in both geometries at both strengths.
- B_kB1y sweep: (E12,qAFM) 0.42 -> 0.88 -> dominant (1.00 at 0.03, hierarchy flips).
- A_W1yy sweep: (qAFM,E12)=1.00 in l2 at all values; (E12,E12) grows with W (0.13->0.23).
- A_W1xz sweep: Fe echo (E12,qAFM) grows; at 0.09 dominates Fe channel (leak).
- A_Km2y sweep: l2 acquires (qAFM,E12)-family + growing contamination by 0.005.

## Full-resolution regeneration (on demand)
Any entry rebuilds at full resolution (tau [-120,0] step 0.1, t [-130,40]) with:
    cd <repo>/build
    python3 ../tfo_project/tmfeo3_2dcs_final/catalogue_single_vertex/regen_fullres.py <RUN> [--keep-h5] [--reader]
Output -> catalogue_single_vertex/fullres/<RUN>/sample_0/{census.json, spectra.npz, panels.png}.
The ~5 GB h5 is deleted unless --keep-h5; --reader also produces the reader_TmFeO3
_components PDFs (--tau-gate 6 --dc-remove --norm linear). Runtime ~1-2 min/entry.
Demo included: fullres/B_kB1y_0.01 (census: (qFM,qAFM)=1.00, (E12,qAFM)=0.68).
