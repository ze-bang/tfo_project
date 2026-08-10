# =====================================================================
# Rebuild the complete published atlas from the CONSTRUCTION PACK ALONE.
# No solver, no run directories, no other file is read.
#
#   python3 make_atlas_from_pack.py [path/to/atlas_construction_pack.h5]
#
# Prints the blind census of each channel (compare with
# figures/census_atlas_kBfree.json) and writes atlas_from_pack.pdf/.png.
# =====================================================================
import sys, os, json
import numpy as np, h5py
from scipy.ndimage import gaussian_filter, maximum_filter, median_filter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# default: the LEAN pack (16x smaller, censuses unchanged); fall back to the full one
_lean, _full = f"{HERE}/atlas_pack_lean.h5", f"{HERE}/atlas_construction_pack.h5"
PACK = sys.argv[1] if len(sys.argv) > 1 else (_lean if os.path.exists(_lean) else _full)
SCALE = 2*np.pi/4.135667696
E_LEVELS = [0.0, 2.067834, 4.9628]
MODES = [("qFM", 0.38), ("E12", 0.50), ("E23", 0.70), ("qAFM", 0.90), ("E13", 1.20)]

# ---- FINAL detection weights (see pack .attrs['detection_operators']) ----
MU, WE1, D_Z23   = 5.264, 5.264, 4.4
C_E12, S_ST, S_DY = 1.488e-4, 0.30, 0.42
C_FE, D13        = 0.0, 6.0

def ws(T):
    w = np.array([1.,0,0]) if T <= 0 else np.array([np.exp(-e/(0.086173*T)) for e in E_LEVELS])
    return w/w.sum()

def spec(t, tau, M, samepol, subtract=True):
    tm = t >= 3.0
    if samepol:
        tk = t[tm]; at = np.exp(np.log(0.03)*((tk-tk[0])/(tk[-1]-tk[0]))**2)
        atau = np.ones_like(tau)
        m1 = tau > -6;   atau[m1] = 0.5*(1-np.cos(np.pi*(-tau[m1])/6))
        m2 = tau < -110; atau[m2] = 0.5*(1-np.cos(np.pi*(120+tau[m2])/10))
    else:
        at = np.hanning(2*tm.sum())[tm.sum():]; atau = np.hanning(2*len(tau))[:len(tau)]
    Md = ((M[:, tm] - M[:, tm].mean(axis=0, keepdims=True)) if subtract else M[:, tm])
    Md = Md*atau[:, None]*at[None, :]
    wt  = np.fft.fftshift(np.fft.fftfreq(tm.sum(), t[1]-t[0]))*SCALE
    wta = np.fft.fftshift(np.fft.fftfreq(len(tau), tau[1]-tau[0]))*SCALE
    mx = (wt > 0.12) & (wt < 1.9); my = np.abs(wta) < 1.6
    return wt[mx], -wta[my], gaussian_filter(
        np.abs(np.fft.fftshift(np.fft.fft2(Md)))[np.ix_(my, mx)], sigma=(1, 0.5))

def blind(wtb, wTb, A, amp_min=0.05):
    B = A.copy(); B[np.abs(wTb) < 0.18, :] = 0
    dy = abs(wTb[1]-wTb[0]); dx = abs(wtb[1]-wtb[0]); n = B.max()
    ismax = (B == maximum_filter(B, size=(int(0.10/dy) | 1, int(0.10/dx) | 1)))
    bg = median_filter(B, size=(int(0.34/dy) | 1, int(0.34/dx) | 1))
    ys, xs = np.where(ismax); pk = []
    for y, x in zip(ys, xs):
        a = B[y, x]/n; p = B[y, x]/max(bg[y, x], 1e-30)
        if a >= amp_min and p >= 1.5:
            pk.append((round(float(a), 3), round(float(p), 1),
                       round(float(wTb[y]), 2), round(float(wtb[x]), 2)))
    pk.sort(reverse=True); out = []
    for q in pk:
        if all(abs(q[2]-r[2]) > 0.09 or abs(q[3]-r[3]) > 0.09 for r in out): out.append(q)
    return out[:8]

with h5py.File(PACK) as f:
    tA   = f["/axes/geomA/t"][:];   tauA = f["/axes/geomA/tau"][:]
    tB   = f["/axes/geomB/t"][:];   tauB = f["/axes/geomB/tau"][:]
    seedsA = ["gs1", "gs2", "gs3"]
    FXM = float(f["/geomA/gs1"].attrs["Fx_mean"])
    def mixA(path, T):
        tot = None
        for wi, s in zip(ws(T), seedsA):
            if wi < 1e-6: continue
            arr = f[f"/geomA/{s}/{path}"][:]
            tot = wi*arr if tot is None else tot + wi*arr
        return tot
    l1A, l2A, l6A = mixA("linear/l1", 10), mixA("linear/l2", 10), mixA("linear/l6", 10)
    FzA, QA, DA   = mixA("linear/Fz", 10), mixA("products/QNL", 10), mixA("products/DNL", 10)
    FxB = f["/geomB/gs1/linear/Fx"][:]; l1B = f["/geomB/gs1/linear/l1"][:]
    l2B = f["/geomB/gs1/linear/l2"][:]; FzB = f["/geomB/gs1/linear/Fz"][:]
    QB  = f["/geomB/gs1/products/QNL"][:]

# --- A cross (O1 -> O2) ---
wt, wT, a2 = spec(tA, tauA, l2A, False, False); _, _, a1 = spec(tA, tauA, l1A, False, False)
_, _, az = spec(tA, tauA, FzA, False, False);   _, _, aq = spec(tA, tauA, QA, False)
lin = MU*a2 + az + WE1*a1
beta = lin[np.argmin(abs(wT-0.90))][np.argmin(abs(wt-0.50))] / \
       aq[np.argmin(abs(wT-0.49))][np.argmin(abs(wt-1.00))]
A_cross = lin + beta*aq
# --- A same (O1 -> O1) ---
_, _, s6 = spec(tA, tauA, l6A, True, False); _, _, s2 = spec(tA, tauA, l2A, True, False)
_, _, sd = spec(tA, tauA, DA,  True, False)
T13c = np.exp(-D13/(1+((wt-1.20)/0.05)**2))[None, :]
T13r = np.exp(-D13/(1+((np.abs(wT)-1.20)/0.05)**2))[:, None]
A_same = D_Z23*s6*T13c*T13r + S_ST*C_E12*s2 + S_DY*(C_E12/FXM)*sd*T13c
# --- B panels ---
wtB, wTB, B_cross = spec(tB, tauB, FxB, False)
_, _, b2 = spec(tB, tauB, l2B, False); _, _, b1 = spec(tB, tauB, l1B, False)
_, _, bz = spec(tB, tauB, FzB, False); _, _, bq = spec(tB, tauB, QB,  False)
B_same = bz + MU*b2 + WE1*b1 + beta*bq

panels = [("(a) A cross   $O_1\\!\\to\\!O_2$", wt, wT, A_cross, "A_cross"),
          ("(b) B cross   $O_2\\!\\to\\!O_1$", wtB, wTB, B_cross, "B_cross"),
          ("(c) A same-pol   $O_1\\!\\to\\!O_1$", wt, wT, A_same, "A_same"),
          ("(d) B same-pol   $O_2\\!\\to\\!O_2$ (prediction)", wtB, wTB, B_same, "B_same_pred")]

census = {}
fig, axs = plt.subplots(2, 2, figsize=(12.4, 9.4))
for ax, (ttl, x, y, Araw, key) in zip(axs.ravel(), panels):
    A = Araw*(1-np.exp(-(((y[:, None]/0.25)**2)+((x[None, :]/0.25)**2))))
    off = A[np.ix_(np.abs(y) > 0.18, x > 0.25)].max()
    rowmax = A[np.ix_(np.abs(y) < 0.10, x > 0.25)].max()
    if rowmax > 0.45*off:
        k = 0.45*off/rowmax
        A = A*(1-(1-k)*np.exp(-((y[:, None]/0.08)**2)))
    pk = blind(x, y, A); census[key] = pk
    print(f"\n=== {key} ===")
    for a_, p, yy, xx in pk: print(f"   {a_:5.2f} x{p:4.1f}  ({yy:+.2f}, {xx:.2f})")
    zref = A[np.ix_(np.abs(y) > 0.18, x > 0.25)].max()
    ax.pcolormesh(x, y, np.minimum(A/zref, 1.0), shading="auto", cmap="inferno",
                  vmin=0, vmax=1, rasterized=True)
    for nm, v in MODES:
        ax.axvline(v, color="w", ls="--", lw=0.5, alpha=0.3)
        ax.axhline(v, color="w", ls="--", lw=0.45, alpha=0.25)
        ax.axhline(-v, color="w", ls="--", lw=0.45, alpha=0.25)
        ax.text(v, 1.37, nm, ha="center", va="top", color="w", alpha=0.7, fontsize=6)
    for a_, p, yy, xx in pk[:6]:
        ax.plot(xx, yy, "x", color="cyan", ms=7, mew=1.4)
        ax.annotate(f"({yy:+.2f},{xx:.2f})", (xx, yy), textcoords="offset points",
                    xytext=(4, -9 if yy > 1.12 else 4), color="cyan", fontsize=6)
    ax.set_xlim(0.15, 1.85); ax.set_ylim(-1.4, 1.4)
    ax.set_xlabel(r"$\omega_t/2\pi$ (THz)"); ax.set_ylabel(r"$\omega_\tau/2\pi$ (THz)")
    ax.set_title(ttl, fontsize=9, loc="left")
fig.suptitle("Atlas rebuilt from the construction pack alone (no solver, no run directories)", fontsize=11)
fig.tight_layout()
fig.savefig(f"{HERE}/atlas_from_pack.pdf", dpi=300)
fig.savefig(f"{HERE}/atlas_from_pack.png", dpi=120)
json.dump(census, open(f"{HERE}/census_from_pack.json", "w"), indent=1)
print("\nwrote atlas_from_pack.{pdf,png} + census_from_pack.json")
