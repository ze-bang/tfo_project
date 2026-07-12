#!/usr/bin/env python3
"""Regenerate any single-vertex catalogue entry at FULL resolution.

Usage (from the repo root, so ../tfo_project paths resolve like the solver expects):
    cd <repo>/build
    python3 ../tfo_project/tmfeo3_2dcs_final/catalogue_single_vertex/regen_fullres.py A_W1yy_0.03 [--keep-h5] [--reader]

Full resolution = tau in [-120,0] step 0.1, t in [-130,40] step 0.02
(vs the catalogue's tau [-60,0] step 0.2).  Output goes to
catalogue_single_vertex/fullres/<RUN>/ with census.json, spectra.npz,
panels.png; the ~5 GB h5 is deleted unless --keep-h5.
"""
import os, sys, math, json, glob, subprocess

RUN = sys.argv[1]
KEEP = "--keep-h5" in sys.argv
READER = "--reader" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(f"{HERE}/../../..")
OUT = f"{HERE}/fullres/{RUN}"

# ---- run table (must match the catalogue generator) ----
def overrides(run):
    g, rest = run.split("_", 1)
    d = 3.0; kw = {}
    if rest == "baseline": pass
    elif rest.startswith("W1yy"): kw["W1_yy"] = float(rest.split("_")[1])
    elif rest.startswith("W1xz"): kw["W1_xz"] = float(rest.split("_")[1])
    elif rest.startswith("kB1y"): kw["kappaB_1y"] = float(rest.split("_")[1])
    elif rest.startswith("Km2y"): kw["Kminus_2y"] = float(rest.split("_")[1])
    elif rest.startswith("kE5x"): kw["kappaE_5x"] = float(rest.split("_")[1])
    elif rest.startswith("d"): d = float(rest[1:])
    else: raise SystemExit(f"unknown run tag: {run}")
    return g, d, kw

def param(g, d, kw):
    if g == "A":
        v = math.sqrt(d*d + 2.3915**2 + 0.9128**2)
        drive = (f"pump_direction   = 1,0,0\nauto_su3_pump      = false\n"
                 f"pump_direction_su3 = 0,{d},0,0,2.3915,0,0.9128,0\n"
                 f"pump_amplitude_su3 = {0.1*0.058333*v:.6f}")
    else:
        drive = ("pump_direction   = 0,0,1\nauto_su3_pump      = true\n"
                 "pump_direction_su3 = 0,0,0,0,0,0,0,0\npump_amplitude_su3 = 0.0")
    body = f"""system = tmfeo3
simulation_mode = 2dcs
lattice_size = 1,1,1
num_trials = 1
J1ab = 4.74
J1c  = 5.15
J2ab = 0.15
J2c  = 0.30
Ka   = -0.0153
Kb   = 0.0
Kc   = -0.0187
D1   = 0.049
D2   = 0.0
e1 = 2.067834
e2 = 4.9628
initial_spin_config = ../tfo_project/tmfeo3_gamma2_1x1x1_seed
T_start = 1e-4
T_end   = 1e-5
annealing_steps  = 0
n_deterministics = 10000
T_zero           = true
spin_length     = 2.5
spin_length_su3 = 1.0
pump_amplitude   = 0.1
pump_width       = 0.5
pump_frequency   = 0.0
pump_time        = 0.0
pump_width_su3     = 0.5
pump_frequency_su3 = 0.0
g_ratio_tm = 0.058333
alpha_gilbert = 0.0
gamma_su3 = 0.0
gamma_su3_lambda4 = 0.15
gamma_su3_lambda5 = 0.15
gamma_su3_lambda6 = 0.15
gamma_su3_lambda7 = 0.15
tau_start    = -120.0
tau_end      =    0.0
tau_step     =    0.1
parallel_tau = true
md_time_start  = -130.0
md_time_end    =  40.0
md_timestep    =   0.02
md_integrator  = rk4
use_gpu        = false
save_observables = false
reuse_m0_for_m1 = false
pulse_window_chunking = false
{drive}
"""
    for k, vv in kw.items(): body += f"{k} = {vv}\n"
    body += f"output_dir       = {OUT}\n"
    return body

g, d, kw = overrides(RUN)
os.makedirs(OUT, exist_ok=True)
pfile = f"{OUT}/fullres.param"
open(pfile, "w").write(param(g, d, kw))
print(f"[regen] running {RUN} at full resolution ...")
r = subprocess.run([f"{REPO}/build/spin_solver", pfile], capture_output=True, text=True, cwd=f"{REPO}/build")
if "completed successfully" not in r.stdout: sys.exit(f"solver failed:\n{r.stdout[-2000:]}\n{r.stderr[-500:]}")

# ---- distill (adaptive tau tapers) ----
import numpy as np, h5py
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SCALE = 2*np.pi/4.135667696
MODES = {"qFM":0.38, "E12":0.50, "E23":0.70, "qAFM":0.90, "E13":1.20}
def name(x):
    b = min(MODES.items(), key=lambda kv: abs(kv[1]-abs(x)))
    return b[0] if abs(b[1]-abs(x)) < 0.07 else f"{x:+.2f}"
CH = [("SU2",0,"FeSx"), ("SU3",1,"l2"), ("SU3",3,"l4"), ("SU3",4,"l5"), ("SU3",5,"l6"), ("SU3",6,"l7")]
dd = f"{OUT}/sample_0"
with h5py.File(f"{dd}/pump_probe_spectroscopy.h5") as f:
    t = f['/reference/times'][:]; tau = f['/tau_scan/tau_values'][:]
    tmin = tau.min()
    out = {}; cen = {}
    for sp, l, nm in CH:
        M0 = f[f'/reference/M_global_{sp}'][:, l]
        M = np.zeros((len(tau), len(t)))
        for i in range(len(tau)):
            gg = f[f'/tau_scan/tau_{i}']
            M[i] = gg[f'M01_global_{sp}'][:, l] - gg[f'M1_global_{sp}'][:, l] - M0
        tm = t >= 3.0; tk = t[tm]
        M = M[:, tm] - M[:, tm].mean(axis=1, keepdims=True)
        apod = np.exp(np.log(0.03)*((tk-tk[0])/(tk[-1]-tk[0]))**2)
        w = np.ones_like(tau)
        e1 = tau > -10; w[e1] = 0.5*(1-np.cos(np.pi*(-tau[e1])/10))
        e2 = tau < tmin+10; w[e2] = 0.5*(1-np.cos(np.pi*(-tmin+tau[e2])/10))
        A = np.abs(np.fft.fftshift(np.fft.fft2(M*w[:, None]*apod[None, :])))
        wt = np.fft.fftshift(np.fft.fftfreq(tm.sum(), t[1]-t[0]))*SCALE
        wta = np.fft.fftshift(np.fft.fftfreq(len(tau), tau[1]-tau[0]))*SCALE
        mx = (wt >= 0) & (wt <= 1.7); my = np.abs(wta) <= 1.7
        A = A[np.ix_(my, mx)].astype(np.float32); wtc = wt[mx]; wtac = wta[my]
        out[nm] = A
        B = A[np.ix_((np.abs(wtac) > 0.10) & (np.abs(wtac) < 1.6), (wtc > 0.15) & (wtc < 1.6))]
        wtb = wtc[(wtc > 0.15) & (wtc < 1.6)]; wtab = wtac[(np.abs(wtac) > 0.10) & (np.abs(wtac) < 1.6)]
        n = float(B.max()); peaks = []; C = B.copy()
        while n > 1e-8:
            iy, ix = np.unravel_index(np.argmax(C), C.shape)
            v = float(C[iy, ix]/n)
            if v < 0.07 or len(peaks) >= 10: break
            peaks.append({"rel": round(v, 3), "wt": round(float(wtb[ix]), 3),
                          "wtau": round(float(wtab[iy]), 3),
                          "label": f"({name(wtab[iy])},{name(wtb[ix])})"})
            dy = int(round(0.07/abs(wtab[1]-wtab[0]))); dx = int(round(0.07/abs(wtb[1]-wtb[0])))
            C[max(0, iy-dy):iy+dy+1, max(0, ix-dx):ix+dx+1] = 0
        cen[nm] = {"max": n, "peaks": peaks}
np.savez_compressed(f"{dd}/spectra.npz", wt=wtc, wta=wtac, **out)
json.dump(cen, open(f"{dd}/census.json", "w"), indent=1)
fig, axs = plt.subplots(1, 6, figsize=(24, 4))
for ax, (sp, l, nm) in zip(axs, CH):
    A = out[nm]; n = cen[nm]["max"] or 1
    ax.pcolormesh(wtc, wtac, A/n, shading="auto", cmap="inferno", vmin=0, vmax=1, rasterized=True)
    ax.set_xlim(0.15, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_title(f"{nm} [{n:.1e}]", fontsize=9)
fig.suptitle(f"{RUN} (FULL RES)"); fig.tight_layout(); fig.savefig(f"{dd}/panels.png", dpi=90)
if READER:
    subprocess.run(["python3", f"{REPO}/util/readers_new/reader_TmFeO3.py", dd, "2dcs",
                    "--e1", "2.067834", "--e2", "4.9628", "--kc", "3.72", "--qfm", "1.571",
                    "--tau-gate", "6", "--dc-remove", "--norm", "linear"], cwd=REPO)
if not KEEP:
    os.remove(f"{dd}/pump_probe_spectroscopy.h5")
    for p in glob.glob(f"{dd}/*.npy"): os.remove(p)
print(f"[regen] done -> {dd}  (h5 {'kept' if KEEP else 'deleted'})")
for nm in ("FeSx", "l2"):
    print(f"  {nm}: " + " ".join(f"{p['label']}={p['rel']}" for p in cen[nm]["peaks"][:4]))
