# =====================================================================
# Build the ATLAS CONSTRUCTION PACK: ONE self-contained HDF5 from which
# the complete 2DCS atlas can be rebuilt with no access to the solver or
# the run directories -- in BOTH domains:
#   /geom*/  time domain      M_NL(tau, t)  + the two product channels
#   /spectra/ frequency domain |FFT2| maps, exactly as published
#   /mixed/   hybrid domain    |FFT_tau| vs detection time t
# Companion reader: make_atlas_from_pack.py (reads ONLY this file).
#
# IMPORTANT (why product channels are stored explicitly): two emitters in
# the final model are BILINEAR in the dynamical coordinates,
#     QNL = (lambda1 * lambda2)      -> on-site SHG at 2 E12
#     DNL = (F_x * lambda2)          -> magnon-CEF SFG/DFG at qAFM +- E12
# and a product of nonlinear signals is NOT a function of the stored
# M_NL: it needs M01, M1 and the reference separately, at build time.
# =====================================================================
import numpy as np, h5py, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = f"{HERE}/trajectories/atlas_construction_pack.h5"
PS_PER_UNIT = 0.6582119569
SCALE = 2*np.pi/4.135667696          # code angular freq -> THz
E_LEVELS = [0.0, 2.067834, 4.9628]   # meV
DEC = 5                              # t decimation for per-component arrays

# ---- FINAL detection weights (necessity-audited minimal operators) ----
MU, WE1_DET, BETA = 5.264, 5.264, 66.0     # O2: mu*l2 + w_E1*l1 + beta*(l1 l2)
D_Z23, D_Z13      = 4.4, 0.0               # O1: 4.4*l6 ; l4 dropped (bound <=0.006)
MU_X23, MU_X13    = 0.9128, 0.0            # O1: 0.913*l7 ; mu_x13 = 0 (dark)
C_E12             = 1.488e-4               # g_d emission conversion (per raw l2 map)
S_STATIC          = 0.30                   # condensed avatar share
S_DYN             = 0.42                   # condensation-ratio correction
C_FE              = 0.0                    # bounded above by the absent magnon diagonal
D13_DEPTH, D12_DEPTH = 6.0, 1.0            # self-absorption optical depths

def ws(T):
    w = np.array([1.,0,0]) if T <= 0 else np.array([np.exp(-e/(0.086173*T)) for e in E_LEVELS])
    return w/w.sum()

def read_run(path, dec_idx):
    """All quantities the atlas needs, in one pass over the delays."""
    with h5py.File(path) as f:
        t   = f['/reference/times'][:]
        tau = f['/tau_scan/tau_values'][:]
        R2  = f['/reference/M_global_SU2'][:]
        R3  = f['/reference/M_global_SU3'][:]
        Rl  = f['/reference/M_local_SU3'][:]
        q0  = Rl[:, 0]*Rl[:, 1]          # (l1 l2) reference
        d0  = R2[:, 0]*Rl[:, 1]          # (F_x l2) reference
        FXM = float(np.abs(R2[:, 0]).mean())
        n_tau, n_t, n_d = len(tau), len(t), len(dec_idx)
        # full-resolution linear coordinates needed by the composites
        L   = {k: np.zeros((n_tau, n_t), np.float32) for k in
               ["Fx", "Fz", "l1", "l2", "l4", "l6", "l7"]}
        QNL = np.zeros((n_tau, n_t), np.float32)     # (l1 l2)  -> SHG
        DNL = np.zeros((n_tau, n_t), np.float32)     # (F_x l2) -> SFG/DFG
        M2  = np.zeros((n_tau, n_d, 3), np.float32)
        M3  = np.zeros((n_tau, n_d, 8), np.float32)
        for i in range(n_tau):
            g  = f[f'/tau_scan/tau_{i}']
            A2 = g['M01_global_SU2'][:]; B2 = g['M1_global_SU2'][:]
            A3 = g['M01_global_SU3'][:]; B3 = g['M1_global_SU3'][:]
            Al = g['M01_local_SU3'][:];  Bl = g['M1_local_SU3'][:]
            d2 = A2 - B2 - R2
            d3 = A3 - B3 - R3
            L["Fx"][i] = d2[:, 0]; L["Fz"][i] = d2[:, 2]
            L["l1"][i] = d3[:, 0]; L["l2"][i] = d3[:, 1]
            L["l4"][i] = d3[:, 3]; L["l6"][i] = d3[:, 5]; L["l7"][i] = d3[:, 6]
            QNL[i] = Al[:, 0]*Al[:, 1] - Bl[:, 0]*Bl[:, 1] - q0
            DNL[i] = A2[:, 0]*Al[:, 1] - B2[:, 0]*Bl[:, 1] - d0
            M2[i]  = d2[dec_idx]; M3[i] = d3[dec_idx]
        ref = {k: f[f'/reference/{k}'][:] for k in
               ["M_global_SU2", "M_global_SU3", "M_local_SU3"]}
    return dict(t=t, tau=tau, L=L, QNL=QNL, DNL=DNL, M2=M2, M3=M3, ref=ref, FXM=FXM)

# ---------------- analysis (identical to figures/atlas_kBfree.py) -----------
def spec(t, tau, M, samepol, subtract=True):
    tm = t >= 3.0
    if samepol:
        tk = t[tm]; at = np.exp(np.log(0.03)*((tk-tk[0])/(tk[-1]-tk[0]))**2)
        atau = np.ones_like(tau)
        m1 = tau > -6;    atau[m1] = 0.5*(1-np.cos(np.pi*(-tau[m1])/6))
        m2 = tau < -110;  atau[m2] = 0.5*(1-np.cos(np.pi*(120+tau[m2])/10))
    else:
        at = np.hanning(2*tm.sum())[tm.sum():]; atau = np.hanning(2*len(tau))[:len(tau)]
    Md = ((M[:, tm] - M[:, tm].mean(axis=0, keepdims=True)) if subtract else M[:, tm])
    Md = Md*atau[:, None]*at[None, :]
    wt  = np.fft.fftshift(np.fft.fftfreq(tm.sum(), t[1]-t[0]))*SCALE
    wta = np.fft.fftshift(np.fft.fftfreq(len(tau), tau[1]-tau[0]))*SCALE
    mx = (wt > 0.12) & (wt < 1.9); my = np.abs(wta) < 1.6
    from scipy.ndimage import gaussian_filter
    A = gaussian_filter(np.abs(np.fft.fftshift(np.fft.fft2(Md)))[np.ix_(my, mx)], sigma=(1, 0.5))
    return wt[mx], -wta[my], A

def mixed_domain(t, tau, M, samepol, dec_idx):
    """|FFT along tau| retaining the detection time t (hybrid representation)."""
    if samepol:
        atau = np.ones_like(tau)
        m1 = tau > -6;    atau[m1] = 0.5*(1-np.cos(np.pi*(-tau[m1])/6))
        m2 = tau < -110;  atau[m2] = 0.5*(1-np.cos(np.pi*(120+tau[m2])/10))
    else:
        atau = np.hanning(2*len(tau))[:len(tau)]
    Md = (M - M.mean(axis=0, keepdims=True))*atau[:, None]
    F  = np.abs(np.fft.fftshift(np.fft.fft(Md, axis=0), axes=0))
    wta = np.fft.fftshift(np.fft.fftfreq(len(tau), tau[1]-tau[0]))*SCALE
    my = np.abs(wta) < 1.6
    return -wta[my], F[my][:, dec_idx].astype(np.float32)

def thermal(packs, key, T, seeds):
    w = ws(T); tot = None
    for wi, s in zip(w, seeds):
        if wi < 1e-6: continue
        arr = packs[s]["L"][key] if key in packs[s]["L"] else packs[s][key]
        tot = wi*arr if tot is None else tot + wi*arr
    return tot

# ---------------------------------------------------------------- build
print("reading runs ...", flush=True)
data = {}
for geom, seeds in [("geomA", ["gs1", "gs2", "gs3"]), ("geomB", ["gs1"])]:
    data[geom] = {}
    for s in seeds:
        p = f"{HERE}/runs/{geom}_{s}/sample_0/pump_probe_spectroscopy.h5"
        with h5py.File(p) as f:
            nt = len(f['/reference/times'][:])
        data[geom][s] = read_run(p, np.arange(0, nt, DEC))
        print(f"  {geom}/{s}", flush=True)

with h5py.File(OUT, "w") as out:
    out.attrs["title"] = ("TmFeO3 2DCS atlas construction pack - FINAL model "
                          "(two-face dressing; dynamical d_z^eff; 2026-08-09)")
    out.attrs["purpose"] = (
        "Self-contained reproduction of the published atlas in BOTH domains. "
        "/geom*: time domain M_NL(tau,t) incl. the two PRODUCT channels; "
        "/spectra: the published 2D |FFT| maps; /mixed: |FFT_tau| vs detection time t. "
        "Rebuild everything from this file alone with make_atlas_from_pack.py.")
    out.attrs["units"] = ("time in code units (1 unit = 0.6582119569 ps) and ps in /axes; "
                          "frequency axes in THz; energies meV, f[THz] = E[meV]/4.135667696")
    out.attrs["MNL_definition"] = ("M_NL = M01 - M1 - M0 per delay (honest three-run subtraction); "
                                   "single-pulse pump-only references in each seed's /reference")
    out.attrs["product_channels"] = (
        "QNL = (lambda1*lambda2)_NL and DNL = (F_x*lambda2)_NL are stored EXPLICITLY because a "
        "product of nonlinear signals cannot be reconstructed from M_NL alone (it needs M01, M1 and "
        "the reference separately). QNL feeds the on-site SHG at 2E12; DNL is the fluctuating part of "
        "the order-induced dipole g_d(<F_x>+dF_x)lambda2 and feeds the magnon-CEF SFG/DFG tones at "
        "qAFM+-E12 = 1.40/0.40 THz.")
    out.attrs["detection_operators"] = (
        "TWO operators, one per analysed polarisation state, serving drive and readout alike. "
        "O2 (E||a,H||c) = F_z + 5.264*l2 + 5.264*l1 + 66*(l1 l2); no l3/l8. Channels: A_cross (F_z=0 "
        "identically there), B_same. "
        "O1 (E||c,H||a) = F_x + 0.9128*l7 + 4.4*l6 + g_d(<F_x>+dF_x)*l2; mu_x13 = 0 (fully dark 1->3), "
        "l4 dropped (d_z(1->3) <= 0.006 bound). Channels: A_same, B_cross. "
        "A_same composite = 4.4*l6*T13 + 0.30*c_E12*l2 + 0.42*(c_E12/<F_x>)*DNL*T13, c_E12 = 1.488e-4, "
        "c_Fe = 0 (bounded above by the non-observation of the (qAFM,qAFM) magnon diagonal). "
        "B_cross is displayed as its Fe part: the measured channel is magnon-dominated, which bounds "
        "the geometry-B Tm-emission normalisation >=5x below the naive geometry-A transfer.")
    out.attrs["propagation_filter"] = (
        "Self-absorption of E1 light only (M1 passes freely): "
        "T(w) = exp(-d/(1+((w-w0)/gamma)^2)), d=6 at E13=1.20 THz (gamma=0.05) applied to the emission "
        "column AND the delay-axis FID; d=1 at E12=0.50 THz (gamma=0.04). Already folded into /spectra.")
    out.attrs["drive"] = (
        "Geometry A: Zeeman B||a f(t) amp 0.12 on Fe; Tm su3 (0,3.0,0,0,0,0,0.9128,0) amp 0.021950 "
        "(E||c d_z^eff 3.0*l2 + B||a mu_x 0.9128*l7). "
        "Geometry B: Zeeman B||c amp 0.10 on Fe; Tm su3 (-0.016500,0.030711,0,...) amp 0.034863 "
        "(E||a d_x*l1 with w_E1 = 0.54 and sign(d_x mu_z) < 0, + B||c mu_z*l2). "
        "kappaB = 0 (exactly inert in these geometries); W1_yy = W1_xz = 0.01; W3_yz = 0.01 + thermal.")
    out.attrs["boltzmann"] = ("Thermal map at T: sum_n w_n * M_n over seeds gs1,gs2,gs3 with "
                              "w_n ~ exp(-E_n/(0.086173*T[K])), E = (0, 2.067834, 4.9628) meV. Exact, "
                              "because the classical lambda dynamics IS the von Neumann equation.")
    out.attrs["atlas_processing"] = (
        "t >= 3 detection window (excludes pulse overlap). Cross channels: Hann on both axes. "
        "Same-pol channel: cosine gate over the overlap region |tau| < 6 plus Gaussian t-apodisation "
        "to 0.03. |FFT2|; physical omega_tau = -fftfreq(tau) (+ = non-rephasing); Gaussian smoothing "
        "sigma = (1, 0.5) bins; LINEAR amplitude display. The omega_tau = 0 pump-probe row is RETAINED "
        "in the geometry-A panels (auto-scaled to the measured relative row amplitude) and EXCLUDED "
        "from all censuses (|omega_tau| < 0.18 THz). Census: 2D local maxima, amp >= 0.05 of map max, "
        "prominence >= 1.5 over a 0.34-THz median background, dedup radius 0.09 THz. "
        "A panels at T = 10 K, B panels at T = 0.")

    # ---------------- axes ----------------
    ax = out.create_group("axes")
    for geom in data:
        s0 = list(data[geom])[0]
        t, tau = data[geom][s0]["t"], data[geom][s0]["tau"]
        td = t[np.arange(0, len(t), DEC)]
        g = ax.create_group(geom)
        g.create_dataset("t", data=t);              g["t"].attrs["units"] = "code units"
        g.create_dataset("t_ps", data=t*PS_PER_UNIT)
        g.create_dataset("t_dec", data=td);         g.create_dataset("t_dec_ps", data=td*PS_PER_UNIT)
        g.create_dataset("tau", data=tau);          g.create_dataset("tau_ps", data=tau*PS_PER_UNIT)
    ax.attrs["note"] = ("t = detection time after the probe; tau = pump-probe delay (negative = pump "
                        "first). *_dec are the DEC-fold decimated grids of the per-component arrays.")

    # ------------- time domain, per geometry and seed -------------
    for geom, seeds in [("geomA", ["gs1", "gs2", "gs3"]), ("geomB", ["gs1"])]:
        gg = out.create_group(geom)
        gg.attrs["seeds"] = ("gs1,gs2,gs3 = trajectories seeded in CEF levels 1,2,3; combine with the "
                             "Boltzmann weights (see /.attrs['boltzmann'])" if len(seeds) > 1
                             else "gs1 only (geometry B is evaluated at T = 0)")
        for s in seeds:
            D  = data[geom][s]
            sg = gg.create_group(s)
            sg.attrs["Fx_mean"] = D["FXM"]
            sg.attrs["Fx_mean_note"] = "<|F_x|> of the pump-only reference: the condensation scale of d_z^eff"
            for k, arr in D["L"].items():
                sg.create_dataset(f"linear/{k}", data=arr, compression="gzip", compression_opts=4,
                                  chunks=(min(64, arr.shape[0]), min(4096, arr.shape[1])))
            sg["linear"].attrs["axes"] = "(tau, t) at FULL time resolution, float32; M_NL of one coordinate"
            sg["linear"].attrs["content"] = ("Fx,Fz = Fe global S_x,S_z; l1,l2,l4,l6,l7 = Gell-Mann "
                                             "coordinates (global frame == local frame for TmFeO3: the four "
                                             "transported Tm ions are exactly equivalent)")
            for k, arr in [("QNL", D["QNL"]), ("DNL", D["DNL"])]:
                sg.create_dataset(f"products/{k}", data=arr, compression="gzip", compression_opts=4,
                                  chunks=(min(64, arr.shape[0]), min(4096, arr.shape[1])))
            sg["products"].attrs["axes"] = "(tau, t) FULL resolution, float32"
            sg["products"].attrs["definition"] = ("QNL = (l1*l2)_M01 - (l1*l2)_M1 - (l1*l2)_ref ; "
                                                  "DNL = (F_x*l2)_M01 - (F_x*l2)_M1 - (F_x*l2)_ref")
            for name, arr in D["ref"].items():
                sg.create_dataset(f"reference/{name}", data=arr, compression="gzip", compression_opts=4)
            sg["reference"].attrs["note"] = "pump-only single-pulse trajectory M0(t) (used in the M_NL subtraction)"
            sg.create_dataset("components/MNL_SU2", data=D["M2"], compression="gzip", compression_opts=4)
            sg.create_dataset("components/MNL_SU3", data=D["M3"], compression="gzip", compression_opts=4)
            sg["components"].attrs["axes"] = (f"(tau, t decimated x{DEC}, component); "
                                              "SU2 = S_x,S_y,S_z; SU3 = lambda1..lambda8 -- the complete "
                                              "coordinate set, for any composite the reader wishes to form")
            print("packed", geom, s, flush=True)

    # ------------- frequency domain: the four published panels -------------
    sp = out.create_group("spectra")
    A, B = data["geomA"], data["geomB"]
    tA, tauA = A["gs1"]["t"], A["gs1"]["tau"]
    tB, tauB = B["gs1"]["t"], B["gs1"]["tau"]
    seedsA = ["gs1", "gs2", "gs3"]
    def mixA(key, T, samepol, subtract):
        w = ws(T); tot = None
        for wi, s in zip(w, seedsA):
            if wi < 1e-6: continue
            arr = A[s]["L"][key] if key in A[s]["L"] else A[s][key]
            tot = wi*arr if tot is None else tot + wi*arr
        return spec(tA, tauA, tot, samepol, subtract)

    # --- A cross: O2 = F_z(=0) + mu l2 + w_E1 l1 + beta (l1 l2)
    wt, wT, a2 = mixA("l2", 10, False, False); _, _, a1 = mixA("l1", 10, False, False)
    _, _, az   = mixA("Fz", 10, False, False); _, _, aq = mixA("QNL", 10, False, True)
    lin  = MU*a2 + az + WE1_DET*a1
    beta = lin[np.argmin(abs(wT-0.90))][np.argmin(abs(wt-0.50))] / \
           aq[np.argmin(abs(wT-0.49))][np.argmin(abs(wt-1.00))]
    A_cross = lin + beta*aq
    # --- A same: minimal O1 with the un-truncated d_z^eff and the E13 filter
    _, _, s6 = mixA("l6", 10, True, False); _, _, s2 = mixA("l2", 10, True, False)
    _, _, sd = mixA("DNL", 10, True, False)
    FXM = A["gs1"]["FXM"]
    T13c = np.exp(-D13_DEPTH/(1+((wt-1.20)/0.05)**2))[None, :]
    T13r = np.exp(-D13_DEPTH/(1+((np.abs(wT)-1.20)/0.05)**2))[:, None]
    A_same = D_Z23*s6*T13c*T13r + S_STATIC*C_E12*s2 + S_DYN*(C_E12/FXM)*sd*T13c
    # --- B cross: Fe m_x (Tm identically zero by closure) ; B same: O2
    wtB, wTB, B_cross = spec(tB, tauB, B["gs1"]["L"]["Fx"], False)
    _, _, b2 = spec(tB, tauB, B["gs1"]["L"]["l2"], False)
    _, _, b1 = spec(tB, tauB, B["gs1"]["L"]["l1"], False)
    _, _, bz = spec(tB, tauB, B["gs1"]["L"]["Fz"], False)
    _, _, bq = spec(tB, tauB, B["gs1"]["QNL"],     False)
    B_same = bz + MU*b2 + WE1_DET*b1 + beta*bq

    sp.create_dataset("omega_t_THz", data=wt);   sp.create_dataset("omega_tau_THz", data=wT)
    sp.create_dataset("omega_t_THz_B", data=wtB); sp.create_dataset("omega_tau_THz_B", data=wTB)
    for name, arr, note in [
        ("A_cross", A_cross, "O1 -> O2, T=10 K; anchor (qAFM,E12)=1.00, SHG (E12,2E12)"),
        ("A_same",  A_same,  "O1 -> O1, T=10 K; anchor (E12,E23)=1.00, SFG/DFG triplet, E13 self-absorbed"),
        ("B_cross", B_cross, "O2 -> O1, T=0; Fe m_x part (magnon-dominated as measured)"),
        ("B_same",  B_same,  "O2 -> O2, T=0; PREDICTION: 2E12 harmonic alone")]:
        d = sp.create_dataset(name, data=arr.astype(np.float32), compression="gzip", compression_opts=4)
        d.attrs["axes"] = "(omega_tau, omega_t) in THz; use *_B axes for the geometry-B panels"
        d.attrs["content"] = note
        d.attrs["normalisation"] = "raw composite amplitude; normalise per panel for display"
    sp.attrs["beta_fitted"] = float(beta)
    sp.attrs["note"] = ("These are the published panels BEFORE the display-only steps (DC-corner taper, "
                        "pump-probe row auto-scale, per-panel normalisation) -- see .attrs['atlas_processing'].")

    # ------------- mixed domain: |FFT_tau| retaining detection time -------------
    mg = out.create_group("mixed")
    decA = np.arange(0, len(tA), DEC); decB = np.arange(0, len(tB), DEC)
    wTm, MA = mixed_domain(tA, tauA, sum(wi*A[s]["L"]["l6"] for wi, s in zip(ws(10), seedsA)), True, decA)
    _,   MB = mixed_domain(tA, tauA, sum(wi*A[s]["L"]["l2"] for wi, s in zip(ws(10), seedsA)), True, decA)
    _,   MD = mixed_domain(tA, tauA, sum(wi*A[s]["DNL"]     for wi, s in zip(ws(10), seedsA)), True, decA)
    _,   MC = mixed_domain(tA, tauA, sum(wi*A[s]["L"]["l1"] for wi, s in zip(ws(10), seedsA)), False, decA)
    wTmB, MFb = mixed_domain(tB, tauB, B["gs1"]["L"]["Fx"], False, decB)
    mg.create_dataset("omega_tau_THz", data=wTm); mg.create_dataset("omega_tau_THz_B", data=wTmB)
    for nm, arr in [("A_l6", MA), ("A_l2", MB), ("A_DNL", MD), ("A_l1", MC), ("B_Fx", MFb)]:
        d = mg.create_dataset(nm, data=arr, compression="gzip", compression_opts=4)
        d.attrs["axes"] = f"(omega_tau [THz], t decimated x{DEC}) -- see /axes/geom*/t_dec_ps"
    mg.attrs["purpose"] = ("Hybrid representation: the delay axis Fourier-transformed (which coherence was "
                           "stored) while the detection axis stays in real time (how the emission decays). "
                           "Directly comparable to a delay-scanned, time-resolved measurement.")

print("wrote", OUT, f"({os.path.getsize(OUT)/1e6:.0f} MB)")
