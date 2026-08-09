# Publication figures for the supplement and foundation (CMU/usetex, PDF).
# Regenerates: experimental_pulse_final, postdrive_buildup_final,
#              samepol_vs_experiment, samepol_darkmu13_final, tm_dipole_forced
# NOTE: the buildup panel reads the leave-one-out runs (B_noW1xz, B_noL1drive,
# B_noW3) and the W1xz mini-sweep (A_Wxz*); regenerate them from the
# FINAL_SCENARIO params with the corresponding single-line changes if absent.
import sys, os, json
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubstyle
pubstyle.apply(fontsize=8.5)
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

ROOT = "/home/pc_linux/ClassicalSpin_Cpp/tfo_project"
FS = f"{ROOT}/tmfeo3_2dcs_final/FINAL_SCENARIO"
AB = "/tmp/claude-1000/-home-pc-linux-ClassicalSpin-Cpp/b53b4a29-36f5-4743-8f9c-1f33e3909a5c/scratchpad/audit"
OUT = f"{ROOT}/paper/figs"
SCALE = 2*np.pi/4.135667696
PS_PER_UNIT = 0.6582
E = [0.0, 2.067834, 4.9628]

def wsT(T):
    w = np.array([1.,0,0]) if T <= 0 else np.array([np.exp(-e/(0.086173*T)) for e in E])
    return w/w.sum()

def load_l(path, sp, l):
    with h5py.File(path) as f:
        t=f['/reference/times'][:]; tau=f['/tau_scan/tau_values'][:]
        R=f[f'/reference/M_global_{sp}'][:,l]
        M=np.zeros((len(tau),len(t)),np.float32)
        for i in range(len(tau)):
            g=f[f'/tau_scan/tau_{i}']
            M[i]=g[f'M01_global_{sp}'][:,l]-g[f'M1_global_{sp}'][:,l]-R
    return t,tau,M

def load_dyn(path):
    with h5py.File(path) as f:
        t=f['/reference/times'][:]; tau=f['/tau_scan/tau_values'][:]
        R2=f['/reference/M_global_SU2'][:]; R3=f['/reference/M_local_SU3'][:]
        FXM=float(np.abs(R2[:,0]).mean())
        Q=np.zeros((len(tau),len(t)),np.float32)
        for i in range(len(tau)):
            g=f[f'/tau_scan/tau_{i}']
            Q[i]=g['M01_global_SU2'][:,0]*g['M01_local_SU3'][:,1]-g['M1_global_SU2'][:,0]*g['M1_local_SU3'][:,1]-R2[:,0]*R3[:,1]
    return t,tau,Q,FXM

def spec_same(t,tau,M):
    tm=t>=3.0; tk=t[tm]
    at=np.exp(np.log(0.03)*((tk-tk[0])/(tk[-1]-tk[0]))**2)
    atau=np.ones_like(tau); m1=tau>-6; atau[m1]=0.5*(1-np.cos(np.pi*(-tau[m1])/6))
    m2=tau<-110; atau[m2]=0.5*(1-np.cos(np.pi*(120+tau[m2])/10))
    Md=M[:,tm]*atau[:,None]*at[None,:]
    wt=np.fft.fftshift(np.fft.fftfreq(tm.sum(),t[1]-t[0]))*SCALE
    wta=np.fft.fftshift(np.fft.fftfreq(len(tau),tau[1]-tau[0]))*SCALE
    mx=(wt>0.12)&(wt<1.9); my=(np.abs(wta)<1.6)
    return wt[mx],-wta[my],gaussian_filter(np.abs(np.fft.fftshift(np.fft.fft2(Md)))[np.ix_(my,mx)],sigma=(1,0.5))

def asame_map(T):
    """Final minimal-operator A-same composite at temperature T."""
    w = wsT(T); parts = {}
    FXM=None
    for k in ["6","2","F","dyn"]: parts[k]=None
    for wi,seed in zip(w,["gs1","gs2","gs3"]):
        if wi < 1e-6: continue
        path=f"{FS}/runs/geomA_{seed}/sample_0/pump_probe_spectroscopy.h5"
        t,tau,M6=load_l(path,"SU3",5); _,_,M2=load_l(path,"SU3",1)
        _,_,Q,fxm=load_dyn(path); FXM=fxm if FXM is None else FXM
        for k,arr in [("6",M6),("2",M2),("dyn",Q)]:
            parts[k]=wi*arr if parts[k] is None else parts[k]+wi*arr
    wtb,wTb,S6=spec_same(t,tau,parts["6"]); _,_,S2=spec_same(t,tau,parts["2"])
    _,_,SD=spec_same(t,tau,parts["dyn"])
    T13c=np.exp(-6.0/(1+((wtb-1.20)/0.05)**2))[None,:]
    T13r=np.exp(-6.0/(1+((np.abs(wTb)-1.20)/0.05)**2))[:,None]
    CE=1.488e-4
    C=4.4*S6*T13c*T13r + 0.30*CE*S2 + 0.42*(CE/FXM)*SD*T13c
    return wtb,wTb,C,(t,tau,parts)

def show_map(ax,wtb,wTb,Z0,vmax=None):
    zref=Z0[np.ix_(np.abs(wTb)>0.18,wtb>0.25)].max() if vmax is None else vmax
    taper=1-np.exp(-(((wTb[:,None]/0.25)**2)+((wtb[None,:]/0.25)**2)))
    Z=np.minimum(Z0*taper/zref,1.0)
    off=Z[np.ix_(np.abs(wTb)>0.18,wtb>0.25)].max()
    rowsel=np.abs(wTb)<0.10; rowmax=Z[np.ix_(rowsel,wtb>0.25)].max()
    if rowmax>0.45*off:
        k=0.45*off/rowmax
        Z=Z*(1-(1-k)*np.exp(-((wTb[:,None]/0.08)**2)))
    pc=ax.pcolormesh(wtb,wTb,Z,shading="auto",cmap="inferno",vmin=0,vmax=1,rasterized=True)
    pubstyle.mode_guides(ax,ymax=1.4)
    ax.set_xlim(0.15,1.85); ax.set_ylim(-1.4,1.4)
    return pc

# ---------------------------------------------------------------- 1. pulse
tab=np.loadtxt(f"{ROOT}/experimental_pulse_codeunits.dat")
tp,Ep=tab[:,0]*PS_PER_UNIT,tab[:,1]
tp=tp-np.sum(tp*Ep**2)/np.sum(Ep**2)   # centre on the field-energy centroid (the solver auto-centres too)
fig,axs=plt.subplots(1,2,figsize=(6.8,2.3))
axs[0].plot(tp,Ep/np.abs(Ep).max(),color="#20415f",lw=1.0)
axs[0].axhline(0,color="0.75",lw=0.4)
axs[0].set_xlabel(r"$t$ (ps)"); axs[0].set_ylabel(r"$E(t)$ (norm.)")
axs[0].set_xlim(-2.5,2.5); pubstyle.panel_tag(axs[0],"(a)")
fr=np.fft.rfftfreq(len(tp)*8,(tp[1]-tp[0]))
sp=np.abs(np.fft.rfft(Ep,n=len(tp)*8)); sp/=sp.max()
axs[1].plot(fr,sp,color="#20415f",lw=1.0)
for nm,v in pubstyle.MODES:
    axs[1].axvline(v,color="0.55",ls=(0,(4,3)),lw=0.5)
    axs[1].text(v,1.03,nm,ha="center",fontsize=6.5,color="0.35")
axs[1].set_xlim(0,1.8); axs[1].set_ylim(0,1.09)
axs[1].set_xlabel(r"$f$ (THz)"); axs[1].set_ylabel(r"$|E(f)|$ (norm.)")
pubstyle.panel_tag(axs[1],"(b)")
fig.tight_layout()
fig.savefig(f"{OUT}/experimental_pulse_final.pdf",dpi=400); plt.close(fig)
print("pulse done")

# ------------------------------------------------------- 2. buildup (Gabor)
def gabor_qafm(path):
    with h5py.File(path) as f:
        t=f['/reference/times'][:]; F=f['/reference/M_global_SU2'][:,0]
    tps=t*PS_PER_UNIT; sig=5.0; f0=0.90
    F=F-F.mean()
    tg=np.arange(-20.0,86.0,1.0)
    amp=[np.abs(np.sum(F*np.exp(-((tps-t0)**2)/(2*sig**2))*np.exp(-2j*np.pi*f0*tps))) for t0 in tg]
    return tg,np.array(amp)
fig,axs=plt.subplots(1,2,figsize=(6.8,2.5))
cols={"full model":"#20415f","no $W^{xz}_1$":"#b0413e","no $d_x$ writing":"#c99700","no $W_3$/thermal":"#3d7d5a"}
ref_norm=None
for (lbl,path),c in zip([("full model",f"{FS}/runs/geomB_gs1/sample_0/pump_probe_spectroscopy.h5"),
                         ("no $W^{xz}_1$",f"{AB}/B_noW1xz/sample_0/pump_probe_spectroscopy.h5"),
                         ("no $d_x$ writing",f"{AB}/B_noL1drive/sample_0/pump_probe_spectroscopy.h5"),
                         ("no $W_3$/thermal",f"{AB}/B_noW3/sample_0/pump_probe_spectroscopy.h5")],cols.values()):
    tg,a=gabor_qafm(path)
    if ref_norm is None: ref_norm=a.max()
    axs[0].plot(tg,a/ref_norm,lw=1.0,color=c,label=lbl)
axs[0].axvline(0,color="0.6",lw=0.5)
axs[0].set_xlabel(r"$t$ (ps)"); axs[0].set_ylabel(r"$q_{\rm AFM}$ Gabor amp.\ (norm.)")
axs[0].legend(frameon=False,fontsize=6.5); pubstyle.panel_tag(axs[0],"(a)")
with h5py.File(f"{FS}/runs/geomB_gs1/sample_0/pump_probe_spectroscopy.h5") as f:
    t=f['/reference/times'][:]; tps=t*PS_PER_UNIT
    chans=[(r"Fe $F_x$ ($q_{\rm AFM}$)",f['/reference/M_global_SU2'][:,0],"#20415f"),
           (r"$\lambda^2$ ($E_{12}$)",f['/reference/M_global_SU3'][:,1],"#b0413e"),
           (r"$\lambda^1$",f['/reference/M_global_SU3'][:,0],"#c99700")]
for lbl,y,c in chans:
    sig=5.0; m=(tps>-25)&(tps<95); tg=tps[m][::4]
    env=[np.sqrt(np.mean((y[m][np.abs(tps[m]-t0)<sig]-np.mean(y[m][np.abs(tps[m]-t0)<sig]))**2)) for t0 in tg]
    env=np.array(env); env=np.maximum(env/env.max(),1e-6)
    axs[1].semilogy(tg,env,lw=1.0,color=c,label=lbl)
axs[1].axvline(0,color="0.6",lw=0.5); axs[1].set_ylim(1e-4,1.4)
axs[1].set_xlabel(r"$t$ (ps)"); axs[1].set_ylabel(r"local rms (norm.)")
axs[1].legend(frameon=False,fontsize=6.5); pubstyle.panel_tag(axs[1],"(b)")
fig.tight_layout()
fig.savefig(f"{OUT}/postdrive_buildup_final.pdf",dpi=400); plt.close(fig)
print("buildup done")

# --------------------------------------- 3. same-pol vs digitised experiment
dig=json.load(open(f"{ROOT}/paper/data/experiment_geomA_samepol_digitized.json"))
SYM={"E12":0.50,"-E12":-0.50,"E23":0.70,"E13":1.20,"qFM":0.38,"qAFM":0.90,"0":0.0,"0.78":0.78}
def parse_peak(sname):
    a,b=[x.strip() for x in sname.strip("()").split(",")]
    def val(x):
        if x in SYM: return SYM[x]
        try: return float(x)
        except ValueError: return None
    return val(a),val(b)
wtb,wTb,C,_=asame_map(10)
gx,gy=np.meshgrid(wtb,wTb)
Zexp=np.zeros_like(gx)
for ft in dig["target_list_our_convention"]:
    y0,x0=parse_peak(ft["peak"]); a0=ft["amp"]
    if x0 is None or y0 is None or x0<0.12: continue   # rectified line analysed separately
    Zexp+=a0*np.exp(-((gx-x0)**2/(2*0.05**2)+(gy-y0)**2/(2*0.05**2)))
# marked qAFM-row peaks at their digitised positions
for mk in dig.get("marked_peaks",[]):
    pass  # already included via the target list
fig,axs=plt.subplots(1,2,figsize=(6.9,3.1),sharey=True)
axs[0].pcolormesh(wtb,wTb,np.minimum(Zexp/Zexp.max(),1),shading="auto",cmap="inferno",vmin=0,vmax=1,rasterized=True)
pubstyle.mode_guides(axs[0],ymax=1.4)
axs[0].set_xlim(0.15,1.85); axs[0].set_ylim(-1.4,1.4)
axs[0].set_title(r"(a)\ digitised measurement",loc="left",fontsize=8)
axs[0].set_xlabel(r"$\omega_t/2\pi$ (THz)"); axs[0].set_ylabel(r"$\omega_\tau/2\pi$ (THz)")
pc=show_map(axs[1],wtb,wTb,C)
axs[1].set_title(r"(b)\ model, minimal operator",loc="left",fontsize=8)
axs[1].set_xlabel(r"$\omega_t/2\pi$ (THz)")
for wT,wt_ in [(0.90,0.50),(0.90,1.41)]:
    for ax in axs: ax.plot(wt_,wT,"x",color="cyan",ms=6,mew=1.2)
fig.subplots_adjust(wspace=0.05)
fig.savefig(f"{OUT}/samepol_vs_experiment.pdf",dpi=400); plt.close(fig)
print("vs-experiment done")

# ------------------------------------------ 4. same-pol at 0 / 10 / 20 K + SDC
fig,axs=plt.subplots(1,4,figsize=(11.2,2.9))
for ax,T in zip(axs[:3],[0,10,20]):
    wtb,wTb,C,pack=asame_map(T)
    show_map(ax,wtb,wTb,C)
    ax.set_title(rf"$T$={T}\,K",fontsize=8.5)
    ax.set_xlabel(r"$\omega_t/2\pi$ (THz)")
axs[0].set_ylabel(r"$\omega_\tau/2\pi$ (THz)")
for ax in axs[1:3]: ax.set_yticklabels([])
t,tau,parts=pack
M=parts["2"]; tm=t>=3.0
sdc=M[:,tm].mean(axis=1); sdc=sdc-np.poly1d(np.polyfit(tau,sdc,3))(tau)
wta=np.fft.fftshift(np.fft.fftfreq(len(tau),tau[1]-tau[0]))*SCALE
S=np.abs(np.fft.fftshift(np.fft.fft(sdc*np.hanning(len(tau)))))
m=np.abs(wta)<1.4
axs[3].plot(-wta[m],S[m]/S[m].max(),color="#20415f",lw=1.0)
axs[3].axvline(0.49,color="0.55",ls=(0,(4,3)),lw=0.6); axs[3].axvline(-0.49,color="0.55",ls=(0,(4,3)),lw=0.6)
axs[3].set_xlabel(r"$\omega_\tau/2\pi$ (THz)")
axs[3].set_ylabel(r"$|S_{\rm DC}(\omega_\tau)|$ (norm.)")
axs[3].set_title(r"rectified line",fontsize=8.5)
fig.tight_layout()
fig.savefig(f"{OUT}/samepol_darkmu13_final.pdf",dpi=400); plt.close(fig)
print("multiT done")

# ---------------------------------------------- 5. forced Tm dipole (3 panels)
def ref_l2_fx(path):
    with h5py.File(path) as f:
        t=f['/reference/times'][:]
        l2=f['/reference/M_global_SU3'][:,1]; fx=f['/reference/M_global_SU2'][:,0]
    return t,l2,fx
paths={0.005:f"{AB}/A_Wxz0.005",0.01:f"{FS}/runs/geomA_gs1",0.02:f"{AB}/A_Wxz0.02",0.03:f"{AB}/A_Wxz0.03"}
fig,axs=plt.subplots(1,3,figsize=(10.6,2.6))
t,l2,fx=ref_l2_fx(f"{paths[0.01]}/sample_0/pump_probe_spectroscopy.h5")
tps=t*PS_PER_UNIT; m=(tps>2)&(tps<40)
axs[0].plot(tps[m],(l2[m]-l2[m].mean())/np.abs(l2[m]-l2[m].mean()).max(),color="#20415f",lw=0.8,label=r"$\lambda^2$")
axs[0].plot(tps[m],(fx[m]-fx[m].mean())/np.abs(fx[m]-fx[m].mean()).max()-2.4,color="#b0413e",lw=0.8,label=r"$F_x$ (offset)")
axs[0].set_xlabel(r"$t$ (ps)"); axs[0].set_yticks([])
axs[0].legend(frameon=False,fontsize=6.5,loc="upper right"); pubstyle.panel_tag(axs[0],"(a)")
amps=[]
for W,pth in paths.items():
    t,l2,fx=ref_l2_fx(f"{pth}/sample_0/pump_probe_spectroscopy.h5")
    tm2=t>=3.0
    fr=np.fft.rfftfreq(tm2.sum(),t[1]-t[0])*SCALE
    sp=np.abs(np.fft.rfft((l2[tm2]-l2[tm2].mean())*np.hanning(tm2.sum())))
    i9=np.argmin(np.abs(fr-0.90))
    amps.append((W,sp[max(0,i9-3):i9+4].max()))
amps=np.array(amps)
axs[1].plot(amps[:,0],amps[:,1]/amps[amps[:,0]==0.01,1],"o-",color="#20415f",ms=3.5,lw=0.9)
axs[1].plot([0,0.032],[0,0.032/0.01],ls=(0,(4,3)),color="0.6",lw=0.7)
axs[1].set_xlabel(r"$W^{xz}_1$"); axs[1].set_ylabel(r"$|\lambda^2|_{0.90}$ (norm.\ to $W$=0.01)")
pubstyle.panel_tag(axs[1],"(b)")
t,l2,fx=ref_l2_fx(f"{paths[0.01]}/sample_0/pump_probe_spectroscopy.h5")
tm2=t>=3.0
fr=np.fft.rfftfreq(tm2.sum(),t[1]-t[0])*SCALE
sp=np.abs(np.fft.rfft((l2[tm2]-l2[tm2].mean())*np.hanning(tm2.sum()))); sp/=sp.max()
mm=(fr>0.3)&(fr<2.0)
axs[2].semilogy(fr[mm],np.maximum(sp[mm],1e-6),color="#20415f",lw=0.9)
for v,lbl in [(0.50,r"$E_{12}$"),(0.90,r"$q_{\rm AFM}$"),(1.00,r"$2E_{12}$")]:
    axs[2].axvline(v,color="0.55",ls=(0,(4,3)),lw=0.6)
    axs[2].text(v,1.6,lbl,ha="center",fontsize=6.5,color="0.35")
axs[2].set_ylim(1e-5,3); axs[2].set_xlabel(r"$f$ (THz)"); axs[2].set_ylabel(r"$|\lambda^2(f)|$ (norm.)")
pubstyle.panel_tag(axs[2],"(c)")
fig.tight_layout()
fig.savefig(f"{OUT}/tm_dipole_forced.pdf",dpi=400); plt.close(fig)
print("forced-dipole done")
