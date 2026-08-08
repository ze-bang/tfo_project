# The RECIPROCAL atlas: exactly TWO operators, each serving excitation AND detection.
#   O2 (E||a,H||c) = S_z + 5.264*l2 + 2.84*l1 (w_E1=0.54 mu) + beta*l1*l2
#       drive:  geometry B su3 vector (-2.84, 5.264)-normalized  [ALREADY the production run]
#       detect: A-cross and B-same
#   O1 (E||c,H||a) = 1.65*S_x + 3.0*l2(d_z^eff) + 0.006*l4 + 4.4*l6 + 0.9128*l7, mu_x13 = 0 (dark)
#       drive:  geometry A su3 vector (0,3.0,0,0.006,0,4.4,0.9128,0)  [runs reciprocal/geomA_gs*]
#       detect: A-same and B-cross
#   1.65 = calibrated Fe/Tm emission scale (the one empirical detection number).
#   Same-pol A processed with the tau-gate; E13 self-absorption (depth 6, both axes)
#   applied to the E1 members of O1 (l4, l6) per the adopted A-same scenario.
import numpy as np, h5py, json, os
from scipy.ndimage import gaussian_filter, maximum_filter, median_filter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SCALE = 2*np.pi/4.135667696
ROOT = "/home/pc_linux/ClassicalSpin_Cpp/tfo_project/tmfeo3_2dcs_final"
HERE = f"{ROOT}/reciprocal"
E_LEV = [0.0, 2.067834, 4.9628]
BETA_GUESS = 66.0

def weights(T):
    if T <= 0: return {"gs1": 1.0}
    w = np.array([np.exp(-e/(0.086173*T)) for e in E_LEV]); w /= w.sum()
    return {f"gs{i+1}": w[i] for i in range(3) if w[i] > 1e-6}

def load_all(path):
    with h5py.File(path) as f:
        t = f['/reference/times'][:]; tau = f['/tau_scan/tau_values'][:]
        R2 = f['/reference/M_global_SU2'][:]; R3 = f['/reference/M_global_SU3'][:]
        Rl = f['/reference/M_local_SU3'][:]; q0 = Rl[:, 0]*Rl[:, 1]
        n_tau, n_t = len(tau), len(t)
        D2 = np.zeros((n_tau, n_t, 3), np.float32); D3 = np.zeros((n_tau, n_t, 8), np.float32)
        Q = np.zeros((n_tau, n_t), np.float32)
        for i in range(n_tau):
            g = f[f'/tau_scan/tau_{i}']
            D2[i] = g['M01_global_SU2'][:] - g['M1_global_SU2'][:] - R2
            D3[i] = g['M01_global_SU3'][:] - g['M1_global_SU3'][:] - R3
            A_ = g['M01_local_SU3'][:]; B_ = g['M1_local_SU3'][:]
            Q[i] = A_[:, 0]*A_[:, 1] - B_[:, 0]*B_[:, 1] - q0
    return t, tau, D2, D3, Q

def spec(t, tau, M, samepol=False):
    tm = t >= 3.0
    if samepol:
        tk = t[tm]; at = np.exp(np.log(0.03)*((tk-tk[0])/(tk[-1]-tk[0]))**2)
        atau = np.ones_like(tau)
        m1 = tau > -6; atau[m1] = 0.5*(1-np.cos(np.pi*(-tau[m1])/6))
        m2 = tau < -110; atau[m2] = 0.5*(1-np.cos(np.pi*(120+tau[m2])/10))
    else:
        at = np.hanning(2*tm.sum())[tm.sum():]; atau = np.hanning(2*len(tau))[:len(tau)]
    Md = (M[:, tm]-M[:, tm].mean(axis=0, keepdims=True))*atau[:, None]*at[None, :]
    wt = np.fft.fftshift(np.fft.fftfreq(tm.sum(), t[1]-t[0]))*SCALE
    wta = np.fft.fftshift(np.fft.fftfreq(len(tau), tau[1]-tau[0]))*SCALE
    mx = (wt > 0.02) & (wt < 1.6); my = (np.abs(wta) < 1.6)
    return wt[mx], -wta[my], gaussian_filter(np.abs(np.fft.fftshift(np.fft.fft2(Md)))[np.ix_(my, mx)], sigma=(2, 1.5))

def blind(wtb, wTb, A, amp_min=0.05):
    B = A.copy(); B[np.abs(wTb) < 0.18, :] = 0
    dyx = abs(wTb[1]-wTb[0]); dxx = abs(wtb[1]-wtb[0])
    n = B.max(); fp = (int(0.10/dyx)|1, int(0.10/dxx)|1)
    ismax = (B == maximum_filter(B, size=fp)); bg = median_filter(B, size=(int(0.34/dyx)|1, int(0.34/dxx)|1))
    pk = []
    for y, x in zip(*np.where(ismax)):
        a = B[y, x]/n; p = B[y, x]/max(bg[y, x], 1e-30)
        if a >= amp_min and p >= 1.5:
            pk.append((round(float(a),3), round(float(p),1), round(float(wTb[y]),2), round(float(wtb[x]),2)))
    pk.sort(reverse=True); out = []
    for q in pk:
        if all(abs(q[2]-r[2]) > 0.09 or abs(q[3]-r[3]) > 0.09 for r in out): out.append(q)
    return out[:8]

# ---- load runs ----
A_runs = {s: load_all(f"{HERE}/geomA_{s}/sample_0/pump_probe_spectroscopy.h5") for s in ["gs1", "gs2", "gs3"]}
B_run  = load_all(f"{ROOT}/FINAL_SCENARIO/runs/geomB_gs1/sample_0/pump_probe_spectroscopy.h5")

def mixA(idx3=None, idx2=None, quad=False, T=10):
    tot = None
    for s, w in weights(T).items():
        t, tau, D2, D3, Q = A_runs[s]
        M = Q if quad else (D2[:, :, idx2] if idx2 is not None else D3[:, :, idx3])
        tot = w*M if tot is None else tot + w*M
    return t, tau, tot

# O2 detection pieces for A-cross (T=10K)
t, tau, l2 = mixA(idx3=1); _, _, l1 = mixA(idx3=0); _, _, Sz = mixA(idx2=2); _, _, Qd = mixA(quad=True)
wtb, wTb, A2m = spec(t, tau, l2); _, _, A1m = spec(t, tau, l1)
_, _, AZ = spec(t, tau, Sz);      _, _, AQ = spec(t, tau, Qd)
lin = 5.264*A2m + 2.84*A1m + AZ
def at2(A, wT, wt_): return A[np.argmin(np.abs(wTb-wT))][np.argmin(np.abs(wtb-wt_))]
beta = at2(lin, 0.90, 0.50)/max(at2(AQ, 0.49, 1.00), 1e-30)   # A-cross parity calibration
Across = lin + beta*AQ
print(f"beta recalibrated from A-cross parity: {beta:.1f}")

# O1 detection for A-same (tau-gated) + E13 self-absorption on E1 members
t, tau, l2g = mixA(idx3=1); _, _, l4g = mixA(idx3=3); _, _, l6g = mixA(idx3=5)
_, _, l7g = mixA(idx3=6);   _, _, Sxg = mixA(idx2=0)
wtg, wTg, G2 = spec(t, tau, l2g, True); _, _, G4 = spec(t, tau, l4g, True)
_, _, G6 = spec(t, tau, l6g, True); _, _, G7 = spec(t, tau, l7g, True); _, _, GF = spec(t, tau, Sxg, True)
T13c = np.exp(-6.0/(1+((wtg-1.20)/0.05)**2))[None, :]
T13r = np.exp(-6.0/(1+((np.abs(wTg)-1.20)/0.05)**2))[:, None]
T12c = np.exp(-1.0/(1+((wtg-0.50)/0.04)**2))[None, :]      # E12-line self-absorption (E||c E1), depth 1
Asame = (0.006*G4 + 4.4*G6)*T13c*T13r + 8.205e-5*GF + 1.488e-4*G2*T12c   # fitted Fe:Tm and d_z^eff scales (g-knob), operator content strict

# O1 detection for B-cross; O2 for B-same (T=0, geometry B run)
tB, tauB, D2B, D3B, QB = B_run
wtB, wTB, B2 = spec(tB, tauB, D3B[:, :, 1]); _, _, B1 = spec(tB, tauB, D3B[:, :, 0])
_, _, B4 = spec(tB, tauB, D3B[:, :, 3]);     _, _, B6 = spec(tB, tauB, D3B[:, :, 5])
_, _, B7 = spec(tB, tauB, D3B[:, :, 6]);     _, _, BFx = spec(tB, tauB, D2B[:, :, 0])
_, _, BFz = spec(tB, tauB, D2B[:, :, 2]);    _, _, BQ = spec(tB, tauB, QB)
T13cB = np.exp(-6.0/(1+((wtB-1.20)/0.05)**2))[None, :]
T12cB = np.exp(-1.0/(1+((wtB-0.50)/0.04)**2))[None, :]
Bcross = (0.006*B4 + 4.4*B6)*T13cB + 8.205e-5*BFx + 1.488e-4*B2*T12cB
Bsame  = BFz + 5.264*B2 + 2.84*B1 + beta*BQ

panels = [
    ("A cross  DRIVE O1 (H\\,$\\parallel$a Zeeman + Tm[3.0$\\lambda^2$,0.006$\\lambda^4$,4.4$\\lambda^6$,0.913$\\lambda^7$])\n"
     "DETECT O2 = S$_z$+5.264$\\lambda^2$+2.84$\\lambda^1$+%.0f$\\lambda^1\\lambda^2$,  $T$=10 K" % beta,
     wtb, wTb, Across, "A_cross"),
    ("B cross  DRIVE O2 (H$\\parallel$c Zeeman + Tm[$-$2.84$\\lambda^1$,5.264$\\lambda^2$])\n"
     "DETECT O1 (fitted cross-family scales),  $T$=0",
     wtB, wTB, Bcross, "B_cross"),
    ("A same-pol  DRIVE O1, DETECT O1 (same operator; $\\tau$-gate;\nE$_{13}$ self-absorption d=6 on the E1 members),  $T$=10 K",
     wtg, wTg, Asame, "A_same"),
    ("B same-pol  DRIVE O2, DETECT O2 (same operator),  $T$=0",
     wtB, wTB, Bsame, "B_same"),
]
census = {"operators": {
    "O2": "S_z + 5.264 l2 + 2.84 l1 (w_E1=0.54 mu) + beta l1 l2; beta recal = %.1f" % beta,
    "O1": "content: S_x + d_z^eff l2 + d_z(l4,l6) + mu_x l7; mu_x13=0. Map weights: 4.4 l6 + 0.006 l4 + 8.205e-5 S_x + 1.488e-4 l2 (fitted cross-family scales, used identically in drive+detection where testable)"},
    "drives": {"A": "0,3.0,0,0.006,0,4.4,0.9128,0 amp 0.037822 (= O1 Tm part)",
               "B": "-0.0165,0.030711 amp 0.034864 (= O2 Tm part, w_E1 ratio)"}}
fig, axs = plt.subplots(2, 2, figsize=(12.6, 9.6))
for ax, (ttl, wx, wy, A, key) in zip(axs.ravel(), panels):
    pk = blind(wx, wy, A); census[key] = pk
    print(f"\n=== {key} ===")
    for a, p, yy, xx in pk: print(f"   {a:5.2f} x{p:4.1f}  ({yy:+.2f}, {xx:.2f})")
    Z = A.copy(); Z[np.abs(wy) < 0.18, :] = 0
    ax.pcolormesh(wx, wy, Z/Z.max(), shading="auto", cmap="inferno", vmin=0, vmax=1, rasterized=True)
    for nm, v in [("qFM",0.38),("E12",0.50),("E23",0.70),("qAFM",0.90),("E13",1.20)]:
        ax.axvline(v, color="w", ls="--", lw=0.6, alpha=0.3)
        ax.axhline(v, color="w", ls="--", lw=0.6, alpha=0.3); ax.axhline(-v, color="w", ls="--", lw=0.6, alpha=0.3)
    for a, p, yy, xx in pk[:6]:
        ax.plot(xx, yy, "x", color="cyan", ms=8, mew=1.6)
        ax.annotate(f"({yy:+.2f},{xx:.2f})", (xx, yy), textcoords="offset points", xytext=(4, 5), color="cyan", fontsize=6.5)
    ax.set_xlim(0.15, 1.55); ax.set_ylim(-1.4, 1.4)
    ax.set_xlabel("$\\omega_t$ (THz)", fontsize=9); ax.set_ylabel("$\\omega_\\tau$ (THz)", fontsize=9)
    ax.set_title(ttl, fontsize=8.3)
fig.suptitle("The RECIPROCAL atlas — two operators, each both excitation and detection\n"
             "O2 $\\equiv$ (E$\\parallel$a, H$\\parallel$c);  O1 $\\equiv$ (E$\\parallel$c, H$\\parallel$a);  $\\kappa^B$-free Hamiltonian",
             fontsize=11.5)
fig.tight_layout()
fig.savefig(f"{HERE}/atlas_reciprocal.png", dpi=125)
json.dump(census, open(f"{HERE}/census_atlas_reciprocal.json", "w"), indent=1)
print("\nsaved reciprocal/atlas_reciprocal.png + census_atlas_reciprocal.json")
