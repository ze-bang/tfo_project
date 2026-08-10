# =====================================================================
# Convert atlas_construction_pack.h5 into conventional formats, in BOTH
# domains, for collaborators who do not want to touch HDF5/Python:
#
#   paper/data/matlab/<geom>_<seed>.mat      time domain, FULL resolution
#   paper/data/matlab/spectra.mat            frequency domain (4 panels + axes)
#   paper/data/csv/<geom>_<seed>_<chan>_MNL.csv.gz   time domain (t decimated)
#   paper/data/csv/spectrum_<panel>.csv.gz          frequency domain maps
#   paper/data/csv/mixed_<name>.csv.gz              |FFT_tau| vs detection time
#   paper/data/csv/*_axes*.csv, *_reference_M0.csv
#
# Single source of truth: the pack. Run after build_atlas_pack.py.
# =====================================================================
import numpy as np, h5py, os, gzip, sys
from scipy.io import savemat

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = f"{HERE}/../../paper/data"
# Default source is the LEAN pack (16x smaller, censuses unchanged). Pass the
# full pack explicitly for full-resolution exports.
PACK = sys.argv[1] if len(sys.argv) > 1 else f"{DATA}/atlas_pack_lean.h5"
MDIR = f"{DATA}/matlab"; CDIR = f"{DATA}/csv"
os.makedirs(MDIR, exist_ok=True); os.makedirs(CDIR, exist_ok=True)
PS = 0.6582119569
# the lean pack is already decimated; do not decimate again
DEC = 1 if "lean" in os.path.basename(PACK) else 5

README = (
 "TmFeO3 2DCS reproduction package - FINAL model (2026-08-09). "
 "Units: time in code units, 1 unit = 0.6582119569 ps (t_ps columns provided); "
 "frequencies in THz. MNL = M01 - M1 - M0 (two-pulse minus the two single-pulse responses). "
 "Matrix orientation: ROWS = tau (pump-probe delay), COLUMNS = t (detection time). "
 "TIME DOMAIN: linear_* are the M_NL of individual coordinates (Fx,Fz = Fe S_x,S_z; "
 "l1,l2,l4,l6,l7 = Gell-Mann coordinates); products_QNL = (lambda1*lambda2)_NL feeds the on-site "
 "SHG at 2*E12; products_DNL = (F_x*lambda2)_NL feeds the magnon-CEF SFG/DFG tones at "
 "qAFM+-E12 = 1.40/0.40 THz. The products are supplied explicitly because a product of nonlinear "
 "signals cannot be reconstructed from MNL alone. "
 "DETECTION (two operators, one per analysed polarisation state): "
 "O2 (E||a,H||c) = Fz + 5.264*l2 + 5.264*l1 + 66*(l1*l2)  -> channels A_cross (Fz=0 there), B_same; "
 "O1 (E||c,H||a) = Fx + 0.9128*l7 + 4.4*l6 + g_d(<Fx>+dFx)*l2, mu_x13 = 0 (dark 1->3) "
 "-> channels A_same, B_cross. A_same composite = 4.4*l6*T13 + 0.30*cE12*l2 "
 "+ 0.42*(cE12/<Fx>)*DNL*T13 with cE12 = 1.488e-4, c_Fe = 0. "
 "PROPAGATION: T(w) = exp(-d/(1+((w-w0)/g)^2)) on E1 emission only, d = 6 at E13 = 1.20 THz "
 "(g = 0.05, on the emission column AND the delay-axis FID), d = 1 at E12 = 0.50 (g = 0.04). "
 "THERMAL: sum over seeds gs1,gs2,gs3 with weights ~ exp(-E_n/(0.086173*T[K])), "
 "E = (0, 2.067834, 4.9628) meV; A panels at T = 10 K, B panels at T = 0. "
 "SPECTRA: t >= 3 window; cross channels Hann on both axes, same-pol cosine gate |tau| < 6 plus "
 "Gaussian t-apodisation to 0.03; |FFT2|; physical omega_tau = -fftfreq(tau) (+ = non-rephasing); "
 "Gaussian smoothing (1, 0.5) bins.")

with h5py.File(PACK) as f:
    # ---------------- time domain ----------------
    for geom in ["geomA", "geomB"]:
        t   = f[f"/axes/{geom}/t"][:]; tau = f[f"/axes/{geom}/tau"][:]
        td  = f[f"/axes/{geom}/t_dec"][:]
        for nm, arr, hdr in [("axis_t_full", t, "t_code_units"), ("axis_t_decimated", td, "t_code_units"),
                             ("axis_tau", tau, "tau_code_units")]:
            np.savetxt(f"{CDIR}/{geom}_{nm}.csv",
                       np.column_stack([arr, arr*PS]), delimiter=",",
                       header=f"{hdr},ps", comments="")
        for seed in [k for k in f[geom] if k.startswith("gs")]:
            sg = f[f"{geom}/{seed}"]
            mat = {"readme": README, "t": t, "t_ps": t*PS, "tau": tau, "tau_ps": tau*PS,
                   "t_decimated": td, "Fx_mean": float(sg.attrs["Fx_mean"])}
            for k in sg["linear"]:
                arr = sg[f"linear/{k}"][:]
                mat[f"linear_{k}"] = arr
                with gzip.open(f"{CDIR}/{geom}_{seed}_linear_{k}_MNL.csv.gz", "wt") as fh:
                    fh.write(f"# {README}\n# rows: tau ({len(tau)}), columns: t decimated x{DEC} ({len(td)})\n")
                    np.savetxt(fh, arr[:, ::DEC], delimiter=",", fmt="%.6e")
            for k in sg["products"]:
                arr = sg[f"products/{k}"][:]
                mat[f"products_{k}"] = arr
                with gzip.open(f"{CDIR}/{geom}_{seed}_products_{k}.csv.gz", "wt") as fh:
                    fh.write(f"# {README}\n# rows: tau ({len(tau)}), columns: t decimated x{DEC} ({len(td)})\n")
                    np.savetxt(fh, arr[:, ::DEC], delimiter=",", fmt="%.6e")
            if "components" in sg:
                mat["components_MNL_SU2"] = sg["components/MNL_SU2"][:]
                mat["components_MNL_SU3"] = sg["components/MNL_SU3"][:]
                mat["components_note"] = ("(tau, t decimated, component); SU2 cols S_x,S_y,S_z; "
                                          "SU3 cols lambda1..lambda8 -- the complete coordinate set")
            for name in ["M_global_SU2", "M_global_SU3", "M_local_SU3"]:
                mat[f"reference_{name}"] = sg[f"reference/{name}"][:]
            savemat(f"{MDIR}/{geom}_{seed}.mat", mat, do_compression=True)
            ref = np.column_stack([t, t*PS, sg["reference/M_global_SU2"][:], sg["reference/M_global_SU3"][:]])
            np.savetxt(f"{CDIR}/{geom}_{seed}_reference_M0.csv", ref, delimiter=",", fmt="%.6e",
                       header="t_code,t_ps,Sx,Sy,Sz,l1,l2,l3,l4,l5,l6,l7,l8", comments="")
            print("exported", geom, seed, flush=True)

    # ---------------- frequency domain ----------------
    sp = f["spectra"]
    smat = {"readme": README,
            "omega_t_THz": sp["omega_t_THz"][:], "omega_tau_THz": sp["omega_tau_THz"][:],
            "omega_t_THz_B": sp["omega_t_THz_B"][:], "omega_tau_THz_B": sp["omega_tau_THz_B"][:],
            "beta_fitted": float(sp.attrs["beta_fitted"]),
            "note": ("2D |FFT| maps as published; rows = omega_tau (physical sign, + = non-rephasing), "
                     "columns = omega_t. Geometry-B panels use the *_B axes.")}
    for panel in ["A_cross", "A_same", "B_cross", "B_same"]:
        arr = sp[panel][:]; smat[f"spectrum_{panel}"] = arr
        wt = sp["omega_t_THz_B" if panel.startswith("B") else "omega_t_THz"][:]
        wT = sp["omega_tau_THz_B" if panel.startswith("B") else "omega_tau_THz"][:]
        with gzip.open(f"{CDIR}/spectrum_{panel}.csv.gz", "wt") as fh:
            fh.write(f"# {README}\n# {sp[panel].attrs['content']}\n")
            fh.write("# first row = omega_t (THz); first column = omega_tau (THz, physical sign)\n")
            top = np.concatenate([[np.nan], wt])
            np.savetxt(fh, np.vstack([top, np.column_stack([wT, arr])]), delimiter=",", fmt="%.6e")
        print("exported spectrum", panel, flush=True)
    savemat(f"{MDIR}/spectra.mat", smat, do_compression=True)

    # ---------------- mixed domain ----------------
    mg = f["mixed"]
    mmat = {"readme": README, "omega_tau_THz": mg["omega_tau_THz"][:],
            "omega_tau_THz_B": mg["omega_tau_THz_B"][:],
            "t_dec_ps_A": f["/axes/geomA/t_dec_ps"][:], "t_dec_ps_B": f["/axes/geomB/t_dec_ps"][:],
            "note": mg.attrs["purpose"]}
    for nm in mg:
        if nm.startswith("omega"): continue
        arr = mg[nm][:]; mmat[f"mixed_{nm}"] = arr
        wT = mg["omega_tau_THz_B" if nm.startswith("B") else "omega_tau_THz"][:]
        tt = f["/axes/geomB/t_dec_ps"][:] if nm.startswith("B") else f["/axes/geomA/t_dec_ps"][:]
        with gzip.open(f"{CDIR}/mixed_{nm}.csv.gz", "wt") as fh:
            fh.write(f"# {README}\n# {mg.attrs['purpose']}\n")
            fh.write("# first row = detection time t (ps); first column = omega_tau (THz)\n")
            top = np.concatenate([[np.nan], tt])
            np.savetxt(fh, np.vstack([top, np.column_stack([wT, arr])]), delimiter=",", fmt="%.6e")
        print("exported mixed", nm, flush=True)
    savemat(f"{MDIR}/mixed_domain.mat", mmat, do_compression=True)

tot = sum(os.path.getsize(os.path.join(d, x)) for d in [MDIR, CDIR] for x in os.listdir(d))
print(f"done: matlab + csv = {tot/1e6:.0f} MB")
