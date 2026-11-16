#!/usr/bin/env python3
"""RANS -> DNS ratio of the Clauser parameter at the same x and the same
pressure gradient, from the zero-pressure-gradient pair:

    python3 beta_ratio.py RANS_ZPG.h5 OUT.txt

beta = delta*/tau_w dP_e/dx. A top boundary that imposes a given dP_e/dx(x)
on both (the inviscid response is the same to first order) gives
beta_DNS/beta_RANS = (delta*_DNS/delta*_RANS) (tau_w,RANS/tau_w,DNS).
Both delta* use the pressure-gradient convention and the same edge (the
free stream just above the layer, wall to delta_99, pg_bl.py). The ratio
assumes the RANS errors of the zero-pressure-gradient layer carry over to
the mild-gradient one.
"""
import os
import sys

import numpy as np
from scipy.signal import savgol_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rans_bl as rb  # noqa: E402
import pg_bl  # noqa: E402


def ratio(rans_path):
    d = rb.dns()
    dt = pg_bl.dns_truncated(d)
    q = pg_bl.pg_integrals(rb.load(rans_path))
    x = np.arange(100.0, 550.0 + 1e-9, 5.0)
    ds_d = np.interp(x, dt["x"], savgol_filter(dt["dstar"], 51, 2))
    tw_d = np.interp(x, d["x"], savgol_filter(0.5 * d["cf"] * d["U"][-1] ** 2, 51, 2))   # export c_f is on U_top
    ds_r = np.interp(x, q["x"], q["dstar"])
    tw_r = np.interp(x, q["x"], q["tauw"])
    return dict(x=x, Re_theta_dns=np.interp(x, d["x"], d["Re_theta"]),
                dstar_ratio=ds_d / ds_r, tauw_ratio=tw_r / tw_d, r=(ds_d / ds_r) * (tw_r / tw_d))


def main():
    src, out = sys.argv[1:3]
    t = ratio(src)
    hdr = ("beta_DNS / beta_RANS at the same x and dP_e/dx, from the ZPG pair\n"
           f"RANS: {os.path.basename(src)};  DNS: ../smooth tbl_dns_smooth.mat (same edge, delta* to delta_99 on both)\n"
           "x  Re_theta(DNS)  dstar_DNS/dstar_RANS  tauw_RANS/tauw_DNS  r = beta_DNS/beta_RANS")
    np.savetxt(out, np.c_[t["x"], t["Re_theta_dns"], t["dstar_ratio"], t["tauw_ratio"], t["r"]],
               fmt=["%.1f", "%.1f", "%.4f", "%.4f", "%.4f"], header=hdr)
    s = (t["x"] >= 200)
    print(f"wrote {out}: r on x 200-550 = {t['r'][s].min():.3f} .. {t['r'][s].max():.3f}, mean {t['r'][s].mean():.3f}")


if __name__ == "__main__":
    main()
