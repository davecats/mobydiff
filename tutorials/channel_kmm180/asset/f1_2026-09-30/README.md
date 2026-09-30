# KMM180 developed statistics, pre-F1 vs F1 (2026-09-30)

The stretched-line (natural y, dyw+ 0.5) re-validation of numerics-review
step 8 (F1, the flux-form viscous stencil). Two identical two-leg runs on the
two A6000s of istmcetus, one per binary: `~/step7b_ref_binaries` (784ee0b,
pre-F1) and `main` at 0737b88 (F1). Both restart from the archived
`channel_kmm180_restart_dywplus0p5.h5` (mobydiff.scalar, t = 5, the tail of a
cold-start transient), run `legA.ini` to t = 35 with statistics OFF (the
transient is discarded), then `legB.ini` from that snapshot to t = 55 with
statistics ON (20 t.u., sampled every 500 steps; 64 000 steps at dt 3.125e-4).
Tutorial settings otherwise (red-black niter 3, sor 1.5, forcing_x = 1).
`run_leg.sh <dir> <binary> <gpu>` is the driver; `analyse.py pre f1` the
comparison; `kmm_f1_stats.png` from `tools/plot_channel_stats.py`.

| quantity | pre-F1 | F1 | diff |
|---|---|---|---|
| U+ centreline | 18.515 | 18.510 | -0.02 % |
| bulk U+ | 15.834 | 15.830 | -0.03 % |
| u'+ peak (y+ 14.8) | 2.652 | 2.623 | -1.1 % |
| v'+ peak (y+ 54) | 0.829 | 0.834 | +0.6 % |
| w'+ peak (y+ 40/38) | 1.074 | 1.083 | +0.9 % |
| -<u'v'>+ peak (y+ 32) | 0.7192 | 0.7193 | +0.02 % |

Max |F1 - pre| over the folded profiles: U+ 0.041 (0.2 % of scale) at y+ 21,
u'+ 0.036, v'+ 0.010, w'+ 0.019, -<u'v'>+ 0.008. The wall-to-wall asymmetry of
ONE run is of the same size, so the two binaries differ by the sampling scatter
of a 20 t.u. window on a chaotic field and nothing systematic; u_tau = 1.0000
in both (the forced channel's wall stress). Against KMM/Moser Re_tau 180
(U+ centreline 18.3, u'+ peak 2.66 at y+ 15, v'+ 0.84, w'+ 1.08, -<u'v'>+ 0.73)
both runs sit within ~1 %, U+ centreline +1.1 %.

LANDMINE (fixed live, in tools/README of remote-hosts memory): two separate
`mpirun -n 1` on istmcetus both bind to core 0 and one rank sits on the wrong
NUMA node -- GPU utilisation 67-71 % instead of 98 %; `taskset -pc` each
rank onto its GPU's node (GPU0 cpus 0-63, GPU1 64-127).
