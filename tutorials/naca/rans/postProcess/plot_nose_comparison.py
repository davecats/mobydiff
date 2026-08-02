#!/usr/bin/env python3
"""Before/after figure for the level-12 nose band: leading-edge Cf of
the level-11 baseline and the nose-refined production case against
OpenFOAM. Reads the two npz files postprocess.sh writes to
assets/referenceStats/ and the OpenFOAM wall-shear sampling; run from
the case directory. Series hues are Okabe-Ito steps, validated for
colour-vision deficiency (all-pairs dE 11.0 deutan), and each curve
also carries its own dash pattern so identity never rests on colour.
"""
import numpy as np, os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
new=np.load("assets/referenceStats/cpcf_c11_nose_final.npz")
old=np.load("assets/referenceStats/cpcf_c11_aoa5_final.npz")
samp="assets/openfoam/postProcessing/surface_sampling/1479"
OF={}
for sd in ("suction","pressure"):
    q=np.loadtxt(os.path.join(samp,f"wallShearStress_{sd}_side.raw"),comments="#")
    o=np.argsort(q[:,0]); q=q[o]; OF[sd]=(q[:,0],2.0*np.linalg.norm(q[:,3:6],axis=1)*np.sign(-q[:,3]))
INK,MUTED="#1a1a1a","#5c5c5c"
fig,axs=plt.subplots(1,2,figsize=(12.2,4.8))
for ax,up,lab in ((axs[0],1.0,"suction side"),(axs[1],0.0,"pressure side")):
    sd="suction" if up==1.0 else "pressure"
    ax.plot(OF[sd][0],OF[sd][1],"--",lw=2.0,color=INK,label="OpenFOAM (body-fitted)",zorder=5)
    for d,c,ls,lb in ((old,"#D55E00","-.","level 11 everywhere (y+ 4.6 at LE)"),
                      (new,"#0072B2","-","+ level-12 nose band (y+ 2.3)")):
        m=d["up"]==up; o=np.argsort(d["xoc"][m])
        ax.plot(d["xoc"][m][o],d["cf"][m][o],ls,lw=1.9,color=c,label=lb,zorder=4)
    ax.axvspan(0,0.056,color="#0072B2",alpha=0.07,zorder=0)
    ax.set_xlim(0,0.14); ax.set_ylim(-0.03 if up==0 else 0,0.045)
    ax.axhline(0,lw=0.8,color=MUTED,alpha=0.5)
    ax.set_title(lab,fontsize=10,color=INK); ax.set_xlabel("x/c",color=INK)
    ax.grid(alpha=0.25,lw=0.6); ax.tick_params(colors=MUTED,labelsize=9)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    for s in ("left","bottom"): ax.spines[s].set_color(MUTED)
axs[0].set_ylabel(r"$C_f$",color=INK); axs[0].legend(fontsize=8.5,labelcolor=INK)
axs[0].annotate("shaded = refined band",xy=(0.062,0.0305),fontsize=7.5,color=MUTED)
fig.suptitle(r"NACA 0012, Re = 4e5, $\alpha$ = 5$^\circ$: leading-edge $C_f$ before/after the level-12 nose band",
             fontsize=11,color=INK)
fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig("assets/figures/cf_nose_refinement.png",dpi=150)
print("wrote assets/figures/cf_nose_refinement.png")
