# Convert atlas_construction_pack.h5 into conventional formats:
#   trajectories/matlab/<geom>_<seed>.mat   (MATLAB v5, compressed; FULL time resolution)
#   trajectories/csv/<geom>_<seed>_<chan>_MNL.csv.gz  (plain text; t decimated x5)
#   trajectories/csv/<geom>_axes*.csv, *_reference.csv
# Single source of truth: the pack. Run after build_atlas_pack.py.
import numpy as np, h5py, os, gzip
from scipy.io import savemat

HERE = os.path.dirname(os.path.abspath(__file__))
PACK = f"{HERE}/trajectories/atlas_construction_pack.h5"
MDIR = f"{HERE}/trajectories/matlab"; CDIR = f"{HERE}/trajectories/csv"
os.makedirs(MDIR, exist_ok=True); os.makedirs(CDIR, exist_ok=True)
DEC = 5

README = (
 "TmFeO3 2DCS atlas construction pack (kappaB-free final scenario, 2026-08-07). "
 "Units: time in code units, 1 unit = 0.6582119569 ps; f[THz] = E[meV]/4.135667696. "
 "MNL = M01 - M1 - M0 (two-pulse minus single-pulse responses). "
 "Channels (detected signals): P2 = Fz + 5.264*lambda2 + 5.264*lambda1 + 66*(l1*l2)_local "
 "[analyser (E||a,H||c): A_cross uses geomA (Fz=0 there), B_same uses geomB; NO lambda3/8]; "
 "mx = Fx + 2.3915*lambda5 + 0.9128*lambda7 [analyser (E||c,H||a): B_cross (pure Fe there), A_same full P1]; "
 "cef_dz = 0.006*lambda4 + 4.4*lambda6 [A_same atlas panel]. "
 "Thermal maps: sum over seeds gs1,gs2,gs3 with weights ~ exp(-E_n/(0.086173*T[K])), "
 "E = (0, 2.067834, 4.9628) meV. "
 "Atlas recipe: keep t>=3; subtract per-t tau-mean; half-Hann apodization on tau and t; "
 "|FFT2|; physical omega_T = -fftfreq(tau); census on local maxima (amp>=0.05, prominence>=1.5). "
 "Matrix orientation: rows = tau (delay), columns = t (detection time).")

with h5py.File(PACK) as f:
    for geom in ["geomA", "geomB"]:
        t   = f[f"{geom}/t"][:]; tau = f[f"{geom}/tau"][:]
        td  = t[::DEC]
        np.savetxt(f"{CDIR}/{geom}_axis_t_full.csv", t, delimiter=",", header="t_code_units", comments="")
        np.savetxt(f"{CDIR}/{geom}_axis_t_decimated.csv", td, delimiter=",", header="t_code_units", comments="")
        np.savetxt(f"{CDIR}/{geom}_axis_tau.csv", tau, delimiter=",", header="tau_code_units", comments="")
        for seed in f[geom]:
            if not seed.startswith("gs"):
                continue
            sg = f[f"{geom}/{seed}"]
            mat = {"readme": README, "t": t, "tau": tau, "t_decimated": td}
            for chan in ["P2", "mx", "cef_dz"]:
                arr = sg[f"channels/{chan}"][:]
                mat[f"MNL_{chan}"] = arr                       # full resolution
                # CSV: decimated in t to keep files manageable
                out = f"{CDIR}/{geom}_{seed}_{chan}_MNL.csv.gz"
                with gzip.open(out, "wt") as fh:
                    fh.write(f"# {README}\n# rows: tau ({len(tau)}), columns: t decimated x{DEC} ({len(td)})\n")
                    np.savetxt(fh, arr[:, ::DEC], delimiter=",", fmt="%.6e")
                print("csv ", out, flush=True)
            mat["components_MNL_SU2"] = sg["components/MNL_SU2"][:]
            mat["components_MNL_SU3"] = sg["components/MNL_SU3"][:]
            mat["components_QNL"]     = sg["components/QNL"][:]
            mat["components_note"]    = f"components on t decimated x{DEC}; SU2 cols S_x,S_y,S_z; SU3 cols lambda1..8"
            for name in ["M_global_SU2", "M_global_SU3", "M_local_SU3"]:
                mat[f"reference_{name}"] = sg[f"reference/{name}"][:]
            savemat(f"{MDIR}/{geom}_{seed}.mat", mat, do_compression=True)
            print("mat  ", f"{MDIR}/{geom}_{seed}.mat", flush=True)
            # reference CSV (small, full resolution)
            ref = np.column_stack([t, sg["reference/M_global_SU2"][:], sg["reference/M_global_SU3"][:]])
            np.savetxt(f"{CDIR}/{geom}_{seed}_reference_M0.csv", ref, delimiter=",", fmt="%.6e",
                       header="t,Sx,Sy,Sz,l1,l2,l3,l4,l5,l6,l7,l8", comments="")
with open(f"{CDIR}/README.txt", "w") as fh:
    fh.write(README + "\n")
print("DONE")
for d in (MDIR, CDIR):
    total = sum(os.path.getsize(os.path.join(d, x)) for x in os.listdir(d))
    print(d, f"{total/1e6:.0f} MB")
