# Reconstruct the full TmFeO3 2DCS atlas from atlas_construction_pack.h5 ALONE.
# Usage: python3 make_atlas_from_pack.py [pack.h5] [T_A_kelvin]
# No solver, no run directories, no other files needed.
import sys, numpy as np, h5py
from scipy.ndimage import gaussian_filter, maximum_filter, median_filter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PACK = sys.argv[1] if len(sys.argv) > 1 else "trajectories/atlas_construction_pack.h5"
T_A  = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
SCALE = 2*np.pi/4.135667696
E_LEV = [0.0, 2.067834, 4.9628]

def weights(T):
    if T <= 0: return {"gs1": 1.0}
    w = np.array([np.exp(-e/(0.086173*T)) for e in E_LEV]); w /= w.sum()
    return {f"gs{i+1}": w[i] for i in range(3) if w[i] > 1e-6}

def mixed_channel(f, geom, chan, T):
    tot = None
    for seed, w in weights(T).items():
        if seed not in f[geom]: continue
        M = f[f"{geom}/{seed}/channels/{chan}"][:]
        tot = w*M if tot is None else tot + w*M
    return f[f"{geom}/t"][:], f[f"{geom}/tau"][:], tot

def spec(t, tau, M):
    tm = t >= 3.0
    at = np.hanning(2*tm.sum())[tm.sum():]; atau = np.hanning(2*len(tau))[:len(tau)]
    Md = (M[:, tm] - M[:, tm].mean(axis=0, keepdims=True))*atau[:, None]*at[None, :]
    wt  = np.fft.fftshift(np.fft.fftfreq(tm.sum(), t[1]-t[0]))*SCALE
    wta = np.fft.fftshift(np.fft.fftfreq(len(tau), tau[1]-tau[0]))*SCALE
    mx = (wt > 0.12) & (wt < 1.6); my = (np.abs(wta) < 1.6)
    return wt[mx], -wta[my], gaussian_filter(np.abs(np.fft.fftshift(np.fft.fft2(Md)))[np.ix_(my, mx)], sigma=(2, 1.5))

def blind(wtb, wTb, A, amp_min=0.05, prom=1.5):
    B = A.copy(); B[np.abs(wTb) < 0.18, :] = 0
    dyx = abs(wTb[1]-wTb[0]); dxx = abs(wtb[1]-wtb[0])
    n = B.max(); fp = (int(0.10/dyx)|1, int(0.10/dxx)|1)
    ismax = (B == maximum_filter(B, size=fp))
    bg = median_filter(B, size=(int(0.34/dyx)|1, int(0.34/dxx)|1))
    pk = []
    for y, x in zip(*np.where(ismax)):
        a = B[y, x]/n; p = B[y, x]/max(bg[y, x], 1e-30)
        if a >= amp_min and p >= prom:
            pk.append((round(float(a), 3), round(float(wTb[y]), 2), round(float(wtb[x]), 2)))
    pk.sort(reverse=True); out = []
    for q in pk:
        if all(abs(q[1]-r[1]) > 0.09 or abs(q[2]-r[2]) > 0.09 for r in out): out.append(q)
    return out[:8]

with h5py.File(PACK) as f:
    panels = [
        (f"A cross ($T$={T_A:.0f} K): P2 operator", *mixed_channel(f, "geomA", "P2", T_A)),
        ("B cross ($T$=0): $m_x$ = Fe $S_x$",       *mixed_channel(f, "geomB", "mx", 0)),
        (f"A same-pol ($T$={T_A:.0f} K): d$_z$ CEF", *mixed_channel(f, "geomA", "cef_dz", T_A)),
        ("B same-pol ($T$=0): P2 operator",          *mixed_channel(f, "geomB", "P2", 0)),
    ]
fig, axs = plt.subplots(2, 2, figsize=(12.4, 9.4))
for ax, (ttl, t, tau, M) in zip(axs.ravel(), panels):
    wtb, wTb, A = spec(t, tau, M)
    pk = blind(wtb, wTb, A)
    print(f"\n=== {ttl}")
    for a, yy, xx in pk: print(f"   {a:5.2f}  ({yy:+.2f}, {xx:.2f})")
    Z = A.copy(); Z[np.abs(wTb) < 0.18, :] = 0
    ax.pcolormesh(wtb, wTb, Z/Z.max(), shading="auto", cmap="inferno", vmin=0, vmax=1, rasterized=True)
    for v in [0.38, 0.50, 0.70, 0.90, 1.20]:
        ax.axvline(v, color="w", ls="--", lw=0.6, alpha=0.3)
        ax.axhline(v, color="w", ls="--", lw=0.6, alpha=0.3); ax.axhline(-v, color="w", ls="--", lw=0.6, alpha=0.3)
    for a, yy, xx in pk[:6]:
        ax.plot(xx, yy, "x", color="cyan", ms=9, mew=1.8)
        ax.annotate(f"({yy:+.2f},{xx:.2f})", (xx, yy), textcoords="offset points", xytext=(5, 6), color="cyan", fontsize=7)
    ax.set_xlim(0.15, 1.6); ax.set_ylim(-1.4, 1.4)
    ax.set_xlabel("$\\omega_t$ (THz)"); ax.set_ylabel("$\\omega_\\tau$ (THz)"); ax.set_title(ttl, fontsize=9.5)
fig.suptitle("TmFeO$_3$ 2DCS atlas — reconstructed from atlas_construction_pack.h5 alone", fontsize=12)
fig.tight_layout(); fig.savefig("atlas_from_pack.png", dpi=120)
print("\nsaved atlas_from_pack.png")
