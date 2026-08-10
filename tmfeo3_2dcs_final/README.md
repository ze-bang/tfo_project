# tmfeo3_2dcs_final

- `FINAL_SCENARIO/` — the working model. `runs/` (production simulation
  output), `params/` (exact solver inputs), `figures/` (the atlas generator
  and the publication figure scripts), `build_atlas_pack.py` +
  `export_conventional.py` (the reproduction-package builders).
  `OPERATORS.md` and `READOUT_CHANNELS.md` document the model.

The collaborator-facing reproduction package — the construction pack, the
MATLAB/CSV exports, the reader and the documentation — lives in
`../paper/data/`. Start there.

Legacy campaign directories (pre-final scenario probes, the reciprocity
campaign, the gate-era archive) were removed 2026-08-10; their verdicts are
recorded in `../tmfeo3_foundation.tex` (necessity audit, code-level
verification) and `../paper/supp.tex`.
