# TmFeO3 2DCS — reproduction package

Everything needed to reproduce every published spectrum, in **both the time
and the frequency domain**, with no solver and no repository.

```bash
python3 make_atlas_from_pack.py          # rebuilds all four atlas panels + censuses
```
That reads **only** `atlas_construction_pack.h5` and writes
`atlas_from_pack.pdf/.png` and `census_from_pack.json`. It reproduces the
published census exactly (largest deviation 0.001 in normalised amplitude,
from float32 storage).

## Contents
| file / dir | what it is |
|---|---|
| `atlas_construction_pack.h5` (800 MB) | **the package** — time domain, frequency domain and hybrid domain, fully self-documented in its own HDF5 attributes |
| `make_atlas_from_pack.py` | the reader above; reads only the pack |
| `READOUT_CHANNELS.md` | **read this first** — exact drive and detection operator per channel, with dataset names |
| `OPERATORS.md` | complete operator inventory: what is in the model, what is bounded, what is excluded and why |
| `matlab/` (786 MB) | `<geom>_<seed>.mat` time domain, `spectra.mat` frequency domain, `mixed_domain.mat`; each carries a `readme` string |
| `csv/` (245 MB) | gzipped CSV of every array (time, frequency and hybrid) plus plain-text axes — for Origin / Igor / Excel |
| `params/` | the exact solver inputs of the five production runs |
| `census_atlas_kBfree.json` | the published blind censuses |
| `experiment_geomA_samepol_digitized.json` | the digitised experimental map used for comparison |

## Pack layout
```
/                      root attrs: model, the two operators, drive, propagation
                       filter, Boltzmann rule, full spectra recipe
/axes/geom{A,B}/       t, t_ps, t_dec, t_dec_ps, tau, tau_ps
/geomA/gs{1,2,3}/      geometry A, three level seeds (Boltzmann-mix for T > 0)
/geomB/gs1/            geometry B (evaluated at T = 0)
   linear/{Fx,Fz,l1,l2,l4,l6,l7}   M_NL(tau, t) per coordinate, FULL resolution
   products/{QNL,DNL}              the two BILINEAR emitters  <-- see below
   components/{MNL_SU2,MNL_SU3}    all coordinates (S_xyz, lambda1..8), t/5
   reference/                      pump-only M0(t)
   attrs: Fx_mean = <|F_x|>
/spectra/              the four published |FFT2| panels + THz axes
/mixed/                |FFT along tau| vs real detection time t
```

### The two product channels — why they are shipped
| dataset | product | radiates at | feature |
|---|---|---|---|
| `products/QNL` | λ¹·λ² | 2·E12 = 1.00 THz | on-site SHG (the β term) |
| `products/DNL` | F_x·λ² | qAFM ± E12 = 1.40 / 0.40 THz | magnon–CEF SFG / DFG |

A product of nonlinear signals is **not** a function of the stored `M_NL`:
`(XY)_NL = X₀₁Y₀₁ − X₁Y₁ − X₀Y₀` needs the two-pulse, probe-only and
pump-only trajectories separately. They are therefore computed when the pack
is built and shipped as their own datasets — without them the SFG/DFG triplet
of the same-polarised channel cannot be reconstructed.

## The model in one paragraph
Tm³⁺ crystal-field levels are non-Kramers singlets, so time reversal blocks
the Fe exchange field from coupling to the CEF coherences at first order: the
magnon–CEF interaction is intrinsically nonlinear, and symmetry admits exactly
**two** second-order objects — an *energy face* `W·SS·λ¹` (one exchange
tensor at one scale, `W1_yy = W1_xz = 0.01`) which converts coherences but
never radiates, and a *moment face* `g_d·E_z(F_x λ²)` (the order-induced
electric dipole) which drives, radiates linearly through its condensed part
`⟨F_x⟩`, and radiates the SFG/DFG tones through its fluctuation `δF_x`.
Plus the bare subsystems, the on-site hyperpolarisability `β`, and one
self-absorption filter. `κB = 0` (exactly inert in these geometries); `W3`
moves no coherent feature (it feeds only the rectified line and the buildup).

Detection is **two operators**, one per analysed polarisation state, each
serving excitation and readout alike:
```
O2 = (E∥a, H∥c) = F_z + 5.264 λ2 + 5.264 λ1 + 66 (λ1λ2)
O1 = (E∥c, H∥a) = F_x + 0.9128 λ7 + 4.4 λ6 + g_d(⟨F_x⟩ + δF_x) λ2
```
Channels: **A-cross = O1→O2, A-same = O1→O1, B-cross = O2→O1, B-same = O2→O2.**

## Regenerating the package from the runs
```bash
mpirun -np 16 build/spin_solver <params/*.param>          # ~40 s each
python3 ../../tmfeo3_2dcs_final/FINAL_SCENARIO/build_atlas_pack.py
python3 ../../tmfeo3_2dcs_final/FINAL_SCENARIO/export_conventional.py
```
The data files (`*.h5`, `*.mat`, `*.csv.gz`) are gitignored — share directly.
