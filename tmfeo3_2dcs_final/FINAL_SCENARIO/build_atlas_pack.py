# Build the ATLAS CONSTRUCTION PACK: one self-contained HDF5 from which the
# complete 2DCS atlas (all four channels, any temperature <~20 K) can be
# reconstructed with no access to the solver or the run directories.
# Companion reader: make_atlas_from_pack.py (reads ONLY the pack).
import numpy as np, h5py, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = f"{HERE}/trajectories/atlas_construction_pack.h5"
BETA, MU, WE1_DET = 66.0, 5.264, 5.264       # detection weights (atlas convention)
E_LEVELS = [0.0, 2.067834, 4.9628]           # meV; Boltzmann weights w_n ~ exp(-E_n/kT)
DEC = 5                                      # t decimation for the per-component group

def channels_from_run(f):
    """Return dict of full-resolution detected-channel M_NL(tau, t) composites."""
    t   = f['/reference/times'][:]
    tau = f['/tau_scan/tau_values'][:]
    R2  = f['/reference/M_global_SU2'][:]
    R3  = f['/reference/M_global_SU3'][:]
    Rl  = f['/reference/M_local_SU3'][:]
    q0  = Rl[:, 0]*Rl[:, 1]
    n_tau, n_t = len(tau), len(t)
    P2   = np.zeros((n_tau, n_t), np.float32)   # (E||a,H||c) analyser
    MX   = np.zeros((n_tau, n_t), np.float32)   # (E||c,H||a) analyser, m_x part
    CEF  = np.zeros((n_tau, n_t), np.float32)   # (E||c,H||a) d_z CEF composite
    for i in range(n_tau):
        g = f[f'/tau_scan/tau_{i}']
        d2 = g['M01_global_SU2'][:] - g['M1_global_SU2'][:] - R2
        d3 = g['M01_global_SU3'][:] - g['M1_global_SU3'][:] - R3
        A_ = g['M01_local_SU3'][:]; B_ = g['M1_local_SU3'][:]
        qn = A_[:, 0]*A_[:, 1] - B_[:, 0]*B_[:, 1] - q0
        # P2 = F_z + mu*l2 + w_E1*l1 + beta*(l1 l2)_local   [NO l3/l8]
        P2[i]  = d2[:, 2] + MU*d3[:, 1] + WE1_DET*d3[:, 0] + BETA*qn
        # m_x  = F_x + 2.3915*l5 + 0.9128*l7
        MX[i]  = d2[:, 0] + 2.3915*d3[:, 4] + 0.9128*d3[:, 6]
        # d_z CEF composite (atlas A-same panel): 0.006*l4 + 4.4*l6
        CEF[i] = 0.006*d3[:, 3] + 4.4*d3[:, 5]
    return t, tau, P2, MX, CEF

def components_from_run(f, td_idx):
    t = f['/reference/times'][:]; tau = f['/tau_scan/tau_values'][:]
    R2 = f['/reference/M_global_SU2'][:]; R3 = f['/reference/M_global_SU3'][:]
    Rl = f['/reference/M_local_SU3'][:]; q0 = Rl[:, 0]*Rl[:, 1]
    M2 = np.zeros((len(tau), len(td_idx), 3), np.float32)
    M3 = np.zeros((len(tau), len(td_idx), 8), np.float32)
    Q  = np.zeros((len(tau), len(td_idx)),    np.float32)
    for i in range(len(tau)):
        g = f[f'/tau_scan/tau_{i}']
        M2[i] = (g['M01_global_SU2'][:] - g['M1_global_SU2'][:] - R2)[td_idx]
        M3[i] = (g['M01_global_SU3'][:] - g['M1_global_SU3'][:] - R3)[td_idx]
        A_ = g['M01_local_SU3'][:]; B_ = g['M1_local_SU3'][:]
        Q[i] = (A_[:, 0]*A_[:, 1] - B_[:, 0]*B_[:, 1] - q0)[td_idx]
    return M2, M3, Q

with h5py.File(OUT, "w") as out:
    out.attrs["title"] = "TmFeO3 2DCS atlas construction pack — kappaB-free final scenario (2026-08-07)"
    out.attrs["purpose"] = ("Self-contained: reconstruct all four atlas channels at any T<~20K from "
                            "this file alone with make_atlas_from_pack.py. All time evolution included.")
    out.attrs["units"] = ("time in code units (1 unit = 0.6582119569 ps); energies meV; "
                          "frequency conversion f[THz] = E[meV]/4.135667696")
    out.attrs["MNL_definition"] = "M_NL = M01 - M1 - M0; single-pulse references in /reference"
    out.attrs["detection_operators"] = (
        "P2 (E||a,H||c) = F_z + 5.264*lambda2 + 5.264*lambda1 + 66*(lambda1*lambda2)_local, NO lambda3/8. "
        "Used by: A_cross (F_z=0 there) and B_same. "
        "P1 (E||c,H||a): m_x = F_x + 2.3915*lambda5 + 0.9128*lambda7 (B_cross; Tm part =0 by closure) "
        "and d_z CEF composite 0.006*lambda4 + 4.4*lambda6 (A_same atlas panel).")
    out.attrs["drive"] = (
        "Geometry A: Zeeman B||a f(t) amp 0.12 on Fe; Tm su3 (0,3.0,0,0,0,0,0.9128,0) amp 0.02195 "
        "(E||c d_z^eff 3.0 lambda2 + B||a mu_x 0.9128 lambda7; dark-mu13). "
        "Geometry B: Zeeman B||c amp 0.10 on Fe; Tm su3 (-0.0165,0.030711,0,...) amp 0.034864 "
        "(E||a d_x lambda1 with w_E1=0.54 and NEGATIVE sign + B||c mu_z lambda2). kappaB = 0.")
    out.attrs["boltzmann"] = ("Thermal map at T: sum_n w_n * M_n over seeds n=gs1,gs2,gs3, "
                              "w_n ~ exp(-E_n/(0.086173*T[K])), E = (0, 2.067834, 4.9628) meV. "
                              "Exact because the qutrit dynamics is the exact von Neumann equation.")
    out.attrs["atlas_processing"] = (
        "Per channel: keep t>=3; subtract per-t tau-mean (removes omega_tau=0 row); "
        "half-Hann apodization on tau (max at tau=0) and on t; |FFT2|; physical omega_T = -fftfreq(tau); "
        "display: Gaussian smoothing sigma=(2,1.5) bins, mask |omega_T|<0.18 THz; "
        "census: local maxima, amp>=0.05 of max, prominence>=1.5 vs 0.34-THz median background, "
        "dedup radius 0.09 THz. A panels at T=10 K, B panels at T=0.")

    runs = {"geomA": ["gs1", "gs2", "gs3"], "geomB": ["gs1"]}
    for geom, seeds in runs.items():
        gg = out.create_group(geom)
        for seed in seeds:
            path = f"{HERE}/runs/{geom}_{seed}/sample_0/pump_probe_spectroscopy.h5"
            with h5py.File(path) as f:
                t, tau, P2, MX, CEF = channels_from_run(f)
                sg = gg.create_group(seed)
                if "t" not in gg:
                    gg.create_dataset("t", data=t)
                    gg.create_dataset("tau", data=tau)
                    gg.create_dataset("t_dec", data=t[np.arange(0, len(t), DEC)])
                for name, arr in [("P2", P2), ("mx", MX), ("cef_dz", CEF)]:
                    sg.create_dataset(f"channels/{name}", data=arr,
                                      compression="gzip", compression_opts=4,
                                      chunks=(min(64, arr.shape[0]), min(4096, arr.shape[1])))
                sg["channels"].attrs["axes"] = "(tau, t) FULL time resolution, float32"
                sg["channels"].attrs["map"] = ("A_cross = P2 of geomA; A_same = cef_dz (atlas) or mx (full P1) of geomA; "
                                               "B_cross = mx of geomB (equals pure F_x there); B_same = P2 of geomB")
                # single-pulse references, full resolution
                for name in ["M_global_SU2", "M_global_SU3", "M_local_SU3"]:
                    sg.create_dataset(f"reference/{name}", data=f[f'/reference/{name}'][:],
                                      compression="gzip", compression_opts=4)
                # raw components, t-decimated
                M2, M3, Q = components_from_run(f, np.arange(0, len(t), DEC))
                sg.create_dataset("components/MNL_SU2", data=M2, compression="gzip", compression_opts=4)
                sg.create_dataset("components/MNL_SU3", data=M3, compression="gzip", compression_opts=4)
                sg.create_dataset("components/QNL",     data=Q,  compression="gzip", compression_opts=4)
                sg["components"].attrs["axes"] = f"(tau, t decimated x{DEC}, component); SU2=S_x,S_y,S_z; SU3=lambda1..8; QNL=(l1*l2)_local"
            print("packed", geom, seed, flush=True)
print("wrote", OUT, f"({os.path.getsize(OUT)/1e6:.0f} MB)")
