#!/usr/bin/env python3
"""Initial condition for the conjugate turbulent-channel validation.

TWO jobs, both of which exist only to shorten a transient nobody measures:

1. THE VELOCITY. Map the developed KMM180 DNS field
   (tutorials/channel_kmm180, 288 x 136 x 288 at Re_tau = 180 in the SAME
   4pi x 2pi box) into the fluid gap of the extended domain, i.e. shift it up
   by y_lo and interpolate onto this case's grid. The source is
   wall-resolved and stretched in y; the target is uniform at dy+ = 1.5, so
   this is a mild filtering in x, z and a re-interpolation in y -- the large
   scales are preserved exactly and only the smallest have to regenerate,
   which costs a few eddy turnovers instead of the t ~ 100 a cold start
   would. The solid is zeroed (the IBM pins it anyway).
   The staggered interpolation is tools/make_channel_restart.py's, the same
   one that built validation/channel_interface/les_ibm's IC.

2. THE SCALARS. Each conjugate scalar gets the steady SERIES-RESISTANCE
   profile it is heading for: linear through each solid slab, and in the
   fluid the Reynolds analogy of the mapped field's OWN mean velocity. With
   antisymmetric wall values +-1 and no source the steady total flux is

       J = 2 / (2 R_s + R_f),   R_s = d/(kappa_s D_f),  R_f = 2 Pr U_c / u_tau^2 * D_f/D_f

   with R_f taken from the analogy (theta and u obey the same constant-flux
   transport up to Pr). The interface value is then theta_i = 1 - J R_s, and
   it differs per scalar only through kappa_s -- so the three kappa_s = 1
   scalars start from the SAME profile, which is what their identical steady
   mean will confirm.

   This is NOT circular with the validation: the profile comes from the run's
   own velocity field through a textbook analogy, not from the reference data
   the gate compares against, and it is deliberately wrong near the wall by
   the Pr-dependent sublayer offset that the run has to find for itself.

    ./make_cht_ic.py cs_1.h5 IC_cht.h5
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
sys.path.insert(0, os.path.join(ROOT, "validation", "scalar"))
import make_channel_restart as mkr                           # noqa: E402
from scalar_tools import BlockGeometry                       # noqa: E402

SRC = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")
Y_LO, Y_HI = 1.0, 3.0            # the grid-aligned interfaces
GAP = Y_HI - Y_LO
RE, PR = 180.0, 0.71
VEL = ("un", "vn", "wn", "pn")
VAR = {"un": 0, "vn": 1, "wn": 2, "pn": 3}
PERIODIC = (True, False, True)

# (name, kappa_s) -- must match cht180.ini's [scalar.N] blocks
SCALARS = [("k1", 1.0), ("k2", 1.0), ("k3", 1.0),
           ("a1", 0.1), ("a2", 10.0), ("a3", 100.0)]


def global_velocity(src, nodes, lengths, var):
    """The y-shifted KMM180 staggered field on `nodes`, solid zeroed."""
    sp = list(mkr.staggered_positions(src["nodes"], var))
    sp[1] = sp[1] + Y_LO
    dp = mkr.staggered_positions(nodes, var)
    g = mkr.interp_field(src["fields"][VEL[var]], sp, dp, PERIODIC, lengths)
    yp = dp[1]
    g[:, (yp < Y_LO) | (yp > Y_HI), :] = 0.0        # array axes are (z, y, x)
    return g


def scalar_profile(yc, umean, kappa_s, d, uc):
    """The steady series-resistance profile for one kappa_s (see the header).

    theta = +1 at y = 0 and -1 at y = ly, linear through each slab, Reynolds
    analogy in between, antisymmetric about the centreline.
    """
    df = 1.0/(RE*PR)                                # fluid diffusivity
    rs = d/(kappa_s*df)                             # one slab's resistance
    rf = 2.0*PR*uc/1.0                              # analogy: J = dtheta_half/(Pr Uc)
    j = 2.0/(2.0*rs + 2.0*rf)
    ti = 1.0 - j*rs                                 # interface temperature
    yc = np.asarray(yc)
    th = np.empty_like(yc)
    lo, hi = yc < Y_LO, yc > Y_HI
    mid = ~(lo | hi)
    th[lo] = 1.0 - (1.0 - ti)*(yc[lo] - (Y_LO - d))/d
    th[hi] = -1.0 + (1.0 - ti)*((Y_HI + d) - yc[hi])/d
    # Reynolds analogy on the fluid: theta falls from +ti to -ti with the
    # shape of the mean velocity, antisymmetrised about the centreline.
    ym = 0.5*(Y_LO + Y_HI)
    shape = np.interp(np.abs(yc[mid] - ym), np.abs(umean[0] - ym), umean[1])
    shape = np.clip(shape/max(uc, 1e-30), 0.0, 1.0)
    th[mid] = np.sign(ym - yc[mid])*ti*shape
    return th, j, ti


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("template", help="a solver-minted 1-step snapshot of THIS case")
    ap.add_argument("out")
    a = ap.parse_args()

    shutil.copyfile(a.template, a.out)
    with h5py.File(a.out, "r+") as f:
        geo = BlockGeometry(f)
        nodes = (f["x"][...], f["y"][...], f["z"][...])
        lengths = tuple(float(f.attrs[k]) for k in ("lx", "ly", "lz"))
        nb = geo.nb[0]
        blocks = f["blocks"][...]
        if int(blocks[:, 3].max()) != 0:
            raise SystemExit("single-level cases only")

        src = mkr.load_source(SRC)
        glob = {nm: global_velocity(src, nodes, lengths, VAR[nm]) for nm in VEL}
        glob["pn"][...] = 0.0                       # the first projection sets it

        for nm in VEL:
            rows = np.zeros(f[nm].shape)
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                rows[bid] = glob[nm][oz:oz + nb, oy:oy + nb, ox:ox + nb]
            f[nm][...] = rows

        # the mapped field's own mean streamwise velocity, per y row
        ynode = nodes[1]
        yc = 0.5*(ynode[:-1] + ynode[1:])
        umean = glob["un"].mean(axis=(0, 2))
        uc = float(umean.max())
        d = Y_LO                                    # slab thickness

        print(f"mapped KMM180 -> {glob['un'].shape}  max|u| = "
              f"{max(np.abs(glob[n]).max() for n in ('un','vn','wn')):.4f}"
              f"   U_centre = {uc:.4f}")
        for nm, ks in SCALARS:
            if nm not in f:
                print(f"   (skipping {nm}: not in the template)")
                continue
            th, j, ti = scalar_profile(yc, (yc, umean), ks, d, uc)
            rows = np.zeros(f[nm].shape)
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                rows[bid] = np.broadcast_to(th[oy:oy + nb][None, :, None],
                                            (nb, nb, nb))
            f[nm][...] = rows
            print(f"   {nm}: kappa_s = {ks:<7g} J = {j:.6f}  theta_interface = {ti:.6f}")

        f.attrs["step"] = np.int32(0)
        f.attrs["t_current"] = 0.0
    print(f"{a.out}: written")


if __name__ == "__main__":
    main()
