# =====================================================================
# Derive the LEAN distribution files from the full construction pack.
#
#   atlas_spectra.h5     (~12 MB)  frequency + hybrid domains only
#   atlas_pack_lean.h5   (~70 MB)  full reproduction capability, decimated
#
# The full pack samples the detection axis at dt = 0.02 code units
# (Nyquist ~38 THz) for physical content below 2 THz -- a ~20x
# oversampling that dominates the file size. The lean pack keeps
# tau/2 (Nyquist 3.8 THz) and t/8 (Nyquist 4.75 THz), a 16x reduction
# that leaves the blind censuses unchanged to +-0.03 in amplitude and
# +-0.01 THz in position (verified against the published censuses).
# =====================================================================
import numpy as np, h5py, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = f"{HERE}/../../paper/data"
FULL = f"{DATA}/atlas_construction_pack.h5"
LEAN = f"{DATA}/atlas_pack_lean.h5"
SPEC = f"{DATA}/atlas_spectra.h5"
DTAU, DT = 2, 8

def copy_attrs(src, dst, extra=None):
    for k, v in src.attrs.items(): dst.attrs[k] = v
    for k, v in (extra or {}).items(): dst.attrs[k] = v

with h5py.File(FULL) as f:
    # ---------------- 1. frequency + hybrid only ----------------
    with h5py.File(SPEC, "w") as o:
        copy_attrs(f, o, {
            "title": "TmFeO3 2DCS — SPECTRA ONLY (frequency + hybrid domains)",
            "purpose": ("The published 2D |FFT2| maps and the hybrid |FFT_tau|-vs-detection-time maps, "
                        "with their axes and the full model documentation. This file has NO time-domain "
                        "trajectories: use atlas_pack_lean.h5 to re-derive the spectra from the raw "
                        "M_NL, or atlas_construction_pack.h5 for the full-resolution trajectories.")})
        f.copy("spectra", o); f.copy("mixed", o)
        g = o.create_group("axes")
        for geom in ["geomA", "geomB"]:
            gg = g.create_group(geom)
            for k in ["t_dec", "t_dec_ps"]:
                gg.create_dataset(k, data=f[f"/axes/{geom}/{k}"][:])
        copy_attrs(f["axes"], o["axes"])
    print(f"wrote {SPEC} ({os.path.getsize(SPEC)/1e6:.1f} MB)")

    # ---------------- 2. lean full-reproduction pack ----------------
    with h5py.File(LEAN, "w") as o:
        copy_attrs(f, o, {
            "title": "TmFeO3 2DCS atlas construction pack — LEAN (decimated; recommended)",
            "decimation": (f"tau/{DTAU} and t/{DT} relative to the full pack (16x smaller). "
                           f"Nyquist: omega_tau 3.8 THz, omega_t 4.75 THz -- far above all physical "
                           f"content (<2 THz). Blind censuses are unchanged to +-0.03 in amplitude "
                           f"and +-0.01 THz in position. Frequency resolution is NOT affected "
                           f"(it is set by the scan ranges, which are untouched)."),
            "purpose": ("Self-contained reproduction of the published atlas in BOTH domains, at the "
                        "sampling the physics actually needs. Rebuild everything with "
                        "make_atlas_from_pack.py. /spectra and /mixed are carried over unchanged from "
                        "the full-resolution computation.")})
        # axes
        ga = o.create_group("axes"); copy_attrs(f["axes"], ga)
        for geom in ["geomA", "geomB"]:
            gg = ga.create_group(geom)
            t = f[f"/axes/{geom}/t"][:][::DT]; tau = f[f"/axes/{geom}/tau"][:][::DTAU]
            gg.create_dataset("t", data=t);          gg.create_dataset("t_ps", data=t*0.6582119569)
            gg.create_dataset("tau", data=tau);      gg.create_dataset("tau_ps", data=tau*0.6582119569)
            gg.create_dataset("t_dec", data=f[f"/axes/{geom}/t_dec"][:])
            gg.create_dataset("t_dec_ps", data=f[f"/axes/{geom}/t_dec_ps"][:])
        # trajectories, decimated
        for geom in ["geomA", "geomB"]:
            gg = o.create_group(geom); copy_attrs(f[geom], gg)
            for seed in [k for k in f[geom] if k.startswith("gs")]:
                sg = f[f"{geom}/{seed}"]; og = gg.create_group(seed); copy_attrs(sg, og)
                for grp in ["linear", "products"]:
                    og.create_group(grp); copy_attrs(sg[grp], og[grp])
                    for k in sg[grp]:
                        og.create_dataset(f"{grp}/{k}", data=sg[f"{grp}/{k}"][:][::DTAU, ::DT],
                                          compression="gzip", compression_opts=4)
                og.create_group("reference"); copy_attrs(sg["reference"], og["reference"],
                    {"sampling": f"decimated x{DT} with the pack t axis"})
                for k in sg["reference"]:
                    # decimate with t so /reference shares the pack's own t axis
                    og.create_dataset(f"reference/{k}", data=sg[f"reference/{k}"][:][::DT],
                                      compression="gzip", compression_opts=4)
                print(f"  lean {geom}/{seed}", flush=True)
        f.copy("spectra", o); f.copy("mixed", o)
    print(f"wrote {LEAN} ({os.path.getsize(LEAN)/1e6:.1f} MB)")
