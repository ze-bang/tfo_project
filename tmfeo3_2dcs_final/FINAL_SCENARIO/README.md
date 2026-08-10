# FINAL MODEL — TmFeO3 2DCS reproduction package (2026-08-09)

Everything needed to reproduce the published atlas, in **both the time and
the frequency domain**, without the solver or this repository.

**Start here:** [`READOUT_CHANNELS.md`](READOUT_CHANNELS.md) — the exact drive
and detection operator of each of the four channels, with dataset names.
Then [`OPERATORS.md`](OPERATORS.md) — the complete operator inventory
(what is in the model, what is bounded, what is excluded and why).

## The model in one paragraph
The Tm³⁺ crystal-field levels are non-Kramers singlets, so time reversal
blocks the Fe exchange field from coupling to the CEF coherences at first
order: the magnon–CEF interaction is intrinsically nonlinear and symmetry
admits exactly **two** second-order objects — an *energy face*
`W·SS·λ¹` (one exchange tensor at one scale, `W1_yy = W1_xz = 0.01`), which
converts coherences but never radiates, and a *moment face*
`g_d·E_z(F_x λ²)`, the order-induced electric dipole, which drives, radiates
linearly through its condensed part `⟨F_x⟩`, and radiates internal
magnon–CEF sum/difference-frequency tones through its fluctuation `δF_x`.
Plus the bare subsystems, the on-site hyperpolarisability `β`, and one
self-absorption filter. `κB = 0` (exactly inert); `W3` moves no coherent
feature (it feeds only the rectified line and the slow buildup).

## THE PACKAGE — `trajectories/atlas_construction_pack.h5` (~840 MB)
Self-contained. Rebuild all four panels and their blind censuses with:

```bash
python3 ../make_atlas_from_pack.py trajectories/atlas_construction_pack.h5
```

That reader touches **only** the pack — no solver, no run directories — and
reproduces the published census exactly (max deviation 0.001 in normalised
amplitude, from float32 storage). It writes `atlas_from_pack.pdf/.png` and
`census_from_pack.json`.

### Layout
```
/                      root attrs: model, operators, drive, propagation filter,
                       Boltzmann rule, and the complete spectra recipe
/axes/geom{A,B}/       t, t_ps, t_dec, t_dec_ps, tau, tau_ps
/geomA/gs{1,2,3}/      geometry A, the three level seeds (Boltzmann-mix for T>0)
/geomB/gs1/            geometry B (evaluated at T = 0)
   linear/{Fx,Fz,l1,l2,l4,l6,l7}   M_NL(tau, t), FULL time resolution
   products/{QNL,DNL}              M_NL of the two BILINEAR emitters (see below)
   components/{MNL_SU2,MNL_SU3}    all coordinates (S_xyz, lambda1..8), t/5
   reference/{M_global_SU2,M_global_SU3,M_local_SU3}   pump-only M0(t)
   attrs: Fx_mean = <|F_x|>, the condensation scale of d_z^eff
/spectra/              FREQUENCY DOMAIN: the four published |FFT2| panels
   A_cross, A_same, B_cross, B_same + omega_t_THz, omega_tau_THz (+ *_B axes)
/mixed/                HYBRID: |FFT along tau| vs detection time t
   A_l6, A_l2, A_DNL, A_l1, B_Fx + omega_tau_THz, use /axes/*/t_dec_ps
```

### The two product channels — why they are shipped
Two emitters in the final model are **bilinear** in the dynamical coordinates:

| dataset | product | radiates at | feature |
|---|---|---|---|
| `products/QNL` | λ¹·λ² | 2·E12 = 1.00 THz | on-site SHG (the β term) |
| `products/DNL` | F_x·λ² | qAFM ± E12 = 1.40 / 0.40 THz | magnon–CEF SFG / DFG |

A product of nonlinear signals is **not** a function of the stored `M_NL`:
`(XY)_NL = X₀₁Y₀₁ − X₁Y₁ − X₀Y₀` needs the two-pulse, probe-only and
pump-only trajectories separately. They are therefore computed when the pack
is built and shipped as their own datasets — without them the SFG/DFG triplet
of the same-polarised channel cannot be reconstructed.

## Conventional formats (no HDF5/Python needed)
Generated from the pack by `export_conventional.py` (single source of truth):

- `trajectories/matlab/` — `<geom>_<seed>.mat` (time domain, full resolution:
  `linear_*`, `products_*`, `components_*`, `reference_*`, axes in both code
  units and ps), plus `spectra.mat` (the four frequency-domain panels + axes)
  and `mixed_domain.mat`. Every file carries a `readme` string with the
  operators, units and recipe — `load('geomA_gs1.mat')` and go.
- `trajectories/csv/` — gzipped CSV of every time-domain array
  (rows = τ, columns = t decimated ×5), of the four spectra and of the mixed
  maps (first row = the ω_t or t axis, first column = ω_τ), plus plain-text
  axes and full-resolution M0 reference CSVs. For Origin / Igor / Excel.

## Regenerating from scratch
```bash
# 1. the runs themselves (~40 s each, from the repo root)
mpirun -np 16 build/spin_solver FINAL_SCENARIO/params/geomA_gs1.param   # etc.
# 2. the pack, then the conventional exports
python3 FINAL_SCENARIO/build_atlas_pack.py
python3 FINAL_SCENARIO/export_conventional.py
# 3. the published figure (reads runs/ directly)
python3 FINAL_SCENARIO/figures/atlas_kBfree.py
```
`params/` holds the exact solver inputs of every production run.
Data files (`*.h5`, `*.mat`, `*.csv.gz`) are gitignored — share directly.

## What the atlas contains (blind censuses, normalised per channel)
| channel | leading peaks | reading |
|---|---|---|
| **A cross** O1→O2 | 1.00 (qAFM, E12); 0.96 (E12, 2E12); 0.50 (−qAFM, E12) | energy-face anchor; on-site SHG; rephasing partner |
| **B cross** O2→O1 | 1.00 (qFM, qAFM); 0.42 (E12, qAFM); 0.15 (qAFM, qAFM) | magnon–magnon (Fe only); electric-write → W transfer; Fe M1 diagonal |
| **A same** O1→O1 | 1.00 (E12, E23); 0.43 (qAFM, 0.40); 0.36 (E12, 1.41); 0.31 (qAFM, 1.41) | population-route anchor; **DFG**; SFG partner; **SFG** |
| **B same** O2→O2 | 1.00 (E12, 2E12), all else ≤ 0.07 | zero-parameter β prediction |

Structural absences that act as measurements: no ω_τ = E13 row (the 1↔3
transition is dark in every polarisation), no (qAFM, qAFM) diagonal in A-same
(bounds the Fe emission weight, `c_Fe = 0`), no strong (E12, 1.20) feature
(excludes the four remaining quadratic moments of the 4c site).

## Provenance
κB-free scenario established 2026-08-06/07 (exhaustive vertex exclusion);
dynamical `d_z^eff` adopted 2026-08-08 after a 33-channel product sweep
identified `δF_x·λ²` as the only emitter whose map is led by the magnon-delay
row; every component then certified by a leave-one-out necessity audit, and
the SU(2)⊗SU(3) implementation verified at code level. Full audit trails in
`../../tmfeo3_foundation.tex` and `../probes_asame/README.md`.
