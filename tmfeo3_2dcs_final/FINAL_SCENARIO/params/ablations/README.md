# Leave-one-out ablation runs (the necessity audit)

Each file is a production param with exactly ONE term removed, used for the
term-by-term verification table in the paper and for the ablation curves of
the buildup figure. Drive-leg ablations rescale the remaining vector so all
other effective fields are unchanged.

| param | term removed |
|---|---|
| `A_noZeeman` | the Fe Zeeman drive leg (geometry A) |
| `A_noL2drive` | the d_z^eff drive leg (3.0 λ²) |
| `A_noL7drive` | the μ_x drive leg (0.9128 λ⁷) |
| `A_noW1yy` | `W1_yy` (the S_y² λ¹ component) |
| `A_noW1xz` | `W1_xz` (the S_xS_z λ¹ component) |
| `A_noKappaB` | `kappaB_1y` — verified bitwise inert |
| `A_Wxz{0.005,0.02,0.03}` | W1_xz sweep (forced-dipole linearity panel) |
| `B_noL1drive` | the d_x electric writing (geometry B) |
| `B_noW1xz` | `W1_xz` (geometry B) |
| `B_noW3` | the `W3` population channel + thermal feedback |
| `B_kappaB` | REINSTATES `kappaB_1y = W = 0.010` (not a removal): (R1,qAFM) rises 0.42 → map max, inverting the excitation-row hierarchy — the intensity-based exclusion of a field-assisted vertex |
| `B_flu0.7` | production drive × 0.7 (field-scaling check: exponents 2.02 / 1.85) |
| `B_flu0.5` | production drive × 0.5 (field-scaling check: exponents 2.02 / 1.88) |
| `B_FeOnly` | W1 = 0 AND Tm drive = 0 (pure Fe sector): (0.38,0.90) survives at 1.18× production — Fe-only origin of qFM→qAFM; qAFM-row excitation profile monotonic through 0.52 — no R1 cross peak |
| `B_noTmDrive` | Tm drive = 0, W1 kept: (0.52,0.90) local max vanishes (profile monotonic) — the R1 cross peak requires the R1 writing, not just the vertex |

Run as `mpirun -np 16 build/spin_solver <file>` (~40 s each); the figure
scripts read the outputs from the paths set inside each param.
