#!/usr/bin/env python3
"""Settle the (K, lambda_sf) -> (kappa_s/kappa_f, rho c_s/rho c_f) mapping of
Neuhauser's conjugate pipe dataset, from the dataset alone.

The result (2026-09-15):

    kappa_s/kappa_f = lambda_sf          ->  [scalar.N] solid_k
    rho c_s/rho c_f = 1/(K^2 lambda_sf)  ->  [scalar.N] solid_rhocp
    alpha_s/alpha_f = K^2 lambda_sf^2

i.e. lambda_sf is the CONDUCTIVITY ratio and K is the FLUID-to-SOLID effusivity
ratio K^2 = (kappa rho c)_f/(kappa rho c)_s -- the inverse of the usual Tiselj
activity ratio, which is why the dataset's "IF corresponds to K -> infinity"
reads backwards at first sight and is right.

Four independent checks, --check {coeff,flux,script,controls,all}:

  coeff     the file's own Nek coefficients.  diffusionCoeff/transportCoeff are
            vdiff/vtrans = CONDUCTIVITY and rho c_p, NOT diffusivity; the fluid
            pair is exactly (nu/Pr, 1), which is what fixes the identification.
  flux      flux continuity in the data: the jump in the phi- and time-averaged
            dT/dr across r = 0.5 must be kappa_s/kappa_f.
  controls  theta'_rms at the interface must be ~0 for MBC (isothermal
            fluctuations -> ibm_wall = dirichlet) and largest for IF (zero
            fluctuating flux -> adiabatic), with the K sweep monotone between.

Needs zarr >= 3 (zarr 2.17 cannot open these files): ~/ibmc/bin/python.
Reads straight out of the distributed tars -- extracting only the arrays it
needs takes under a second and ~1 GB of scratch, not the full 7.2 GB.
"""
import argparse, os, subprocess, sys
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "neuhauser_data/10.35097-26za20q32xsz43yk/data/dataset")
NU = 1.0 / 5300.0
IPR = 1          # Pr axis: [0.025, 0.71]


def extract(case, which, names, dest):
    """Pull just `names` out of the case's tar. zarr v3 stores one directory
    per array, so a member name is the array name."""
    stem = "joined_datasets.interp.zarr" if which == "interp" else "joined_datasets.zarr"
    tar = os.path.join(DATA, case, stem + ".tar")
    root = os.path.join(dest, stem)
    todo = [n for n in names if not os.path.isdir(os.path.join(root, n))]
    if todo or not os.path.isfile(os.path.join(root, "zarr.json")):
        os.makedirs(dest, exist_ok=True)
        subprocess.run(["tar", "-xf", tar, "-C", dest, stem + "/zarr.json"]
                       + [f"{stem}/{n}" for n in todo], check=True)
    return root


def load(root, name):
    import zarr
    return zarr.open_array(os.path.join(root, name), mode="r")[...]


def coords(root):
    return {n: load(root, n) for n in ("type", "K", "lambda_sf", "bccode", "isc")}


def check_coeff(scratch):
    root = extract("cht_short", "joined", [
        "type", "K", "lambda_sf", "bccode", "isc", "Pr",
        "diffusionCoeff", "diffusionCoeffSolid",
        "transportCoeff", "transportCoeffSolid"], scratch)
    c = coords(root)
    Pr = load(root, "Pr")
    dF, dS = load(root, "diffusionCoeff"), load(root, "diffusionCoeffSolid")
    tF, tS = load(root, "transportCoeff"), load(root, "transportCoeffSolid")

    print("--- check 1: the file's own Nek coefficients ---")
    print("fluid (diffusionCoeff, transportCoeff) vs (nu/Pr, 1):")
    ok = True
    for ip, pr in enumerate(Pr):
        good = np.isclose(dF[0, ip], NU / pr) and tF[0, ip] == 1
        ok &= good
        print(f"  Pr={pr:<6} {dF[0,ip]:.8e} vs {NU/pr:.8e}, rho c_p = {tF[0,ip]}"
              f"   {'OK' if good else 'MISMATCH'}")
    print("  => diffusionCoeff is the CONDUCTIVITY (vdiff), not the diffusivity.\n")

    hdr = (f"{'isc':>3} {'type':5} {'K':>5} {'lam':>5} | {'k_s/k_f':>9}"
           f" {'rc_s/rc_f':>10} {'a_s/a_f':>9} | k_s/k_f == lam   rc == 1/(K^2 lam)")
    print(hdr)
    for i in np.flatnonzero(c["bccode"] == 0):
        if str(c["type"][i]) != "CHT":
            continue
        for ip in range(len(Pr)):
            kr, cr = dS[i, ip] / dF[i, ip], tS[i, ip] / tF[i, ip]
            a, b = np.isclose(kr, c["lambda_sf"][i]), \
                   np.isclose(cr, 1.0 / (c["K"][i] ** 2 * c["lambda_sf"][i]))
            ok &= a and b
            print(f"{c['isc'][i]:>3} {str(c['type'][i]):5} {c['K'][i]:>5}"
                  f" {c['lambda_sf'][i]:>5} | {kr:>9.4f} {cr:>10.4f} {kr/cr:>9.4f} |"
                  f" {str(a):>13}   {str(b):>16}   (Pr={Pr[ip]})")
    return ok


def radial_gradient(root, ir, i, phi):
    """phi- and time-averaged dT/dr at radial index ir for scalar i."""
    import zarr
    gx = zarr.open_array(os.path.join(root, "d(t)dx"), mode="r")[:, ir, :, i, IPR]
    gy = zarr.open_array(os.path.join(root, "d(t)dy"), mode="r")[:, ir, :, i, IPR]
    return np.nanmean(gx * np.cos(phi)[None, :] + gy * np.sin(phi)[None, :])


def check_flux(scratch):
    root = extract("cht_short", "interp", [
        "type", "K", "lambda_sf", "bccode", "isc", "r", "phi",
        "d(t)dx", "d(t)dy"], scratch)
    c, r, phi = coords(root), load(root, "r"), load(root, "phi")
    i_in = int(np.max(np.flatnonzero(r < 0.5)))
    i_out = int(np.min(np.flatnonzero(r > 0.5)))
    print("\n--- check 2: flux continuity across the interface ---")
    print(f"r = {r[i_in]:.8f} (fluid) and {r[i_out]:.8f} (solid)")
    print(f"{'isc':>3} {'type':5} {'K':>5} {'lam':>5} | {'dT/dr|f':>11} {'dT/dr|s':>11}"
          f" | {'ratio':>8} {'lam_sf':>7} {'K*sqrt(lam)':>11}  verdict")
    ok = True
    for i in np.flatnonzero(c["bccode"] == 0):
        if str(c["type"][i]) != "CHT":
            continue
        gf = radial_gradient(root, i_in, i, phi)
        gs = radial_gradient(root, i_out, i, phi)
        ratio, lam, K = gf / gs, c["lambda_sf"][i], c["K"][i]
        # 1e-3 covers the polar interpolation residual; the rival reading
        # differs by a factor 2 or 4, so the test is not delicate.
        good = abs(ratio / lam - 1) < 1e-3
        ok &= good
        print(f"{c['isc'][i]:>3} {str(c['type'][i]):5} {K:>5} {lam:>5} |"
              f" {gf:>11.4f} {gs:>11.4f} | {ratio:>8.4f} {lam:>7} {K*np.sqrt(lam):>11.4f}"
              f"  {'lam_sf' if good else 'NEITHER'}")
    return ok


def check_controls(scratch):
    root = extract("cht_short", "interp", [
        "type", "K", "lambda_sf", "bccode", "isc", "r", "phi", "t", "t*t"], scratch)
    import zarr
    c, r = coords(root), load(root, "r")
    i_in = int(np.max(np.flatnonzero(r < 0.5)))
    T = zarr.open_array(os.path.join(root, "t"), mode="r")
    TT = zarr.open_array(os.path.join(root, "t*t"), mode="r")
    print("\n--- check 3: the controls bracket the K sweep ---")
    print(f"{'isc':>3} {'type':6} {'K':>5} {'lam':>5} | {'theta_rms':>11} {'/isc0':>7}")
    rms, base = {}, None
    for i in np.flatnonzero(c["bccode"] == 0):
        t1, t2 = T[:, i_in, :, i, IPR], TT[:, i_in, :, i, IPR]
        v = np.sqrt(max(np.nanmean(t2 - t1 ** 2), 0.0))
        base = v if base is None else base
        rms[int(c["isc"][i])] = v / base
        print(f"{c['isc'][i]:>3} {str(c['type'][i]):6} {c['K'][i]:>5}"
              f" {c['lambda_sf'][i]:>5} | {v:>11.4f} {v/base:>7.3f}")
    # MBC -> isothermal fluctuations -> dirichlet; IF -> zero fluctuating flux
    # -> adiabatic; and the K sweep must sit monotonically between them.
    ok = rms[5] < 0.05 and rms[6] > 1.5 and rms[4] < rms[0] < rms[3]
    print(f"  MBC ~ 0 (dirichlet), IF largest (adiabatic), K sweep monotone: "
          f"{'OK' if ok else 'FAILED'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", default="all",
                    choices=["coeff", "flux", "controls", "all"])
    ap.add_argument("--scratch", default="/tmp/neuhauser_pin_down",
                    help="where the extracted arrays land (kept between runs)")
    args = ap.parse_args()

    checks = ["coeff", "flux", "controls"] if args.check == "all" else [args.check]
    ok = all({"coeff": check_coeff, "flux": check_flux,
              "controls": check_controls}[k](args.scratch) for k in checks)
    print("\nMAPPING: solid_k = lambda_sf,  solid_rhocp = 1/(K^2 * lambda_sf),"
          "  alpha_s/alpha_f = K^2 * lambda_sf^2")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
