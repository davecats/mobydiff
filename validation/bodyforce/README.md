# Body-force gates: the correction kernel visits the force's SUPPORT only

The volumetric body force (2026-10-08; `bodyforce_type sup_blocks` /
`sup_var_lo:hi`) is applied over its support only -- the trip (the CaNS/
SIMSON trip of the boundary-layer tutorial) touches `v` on its block list,
every other source the whole domain. The profile of HoreKa job 5184720 had
the whole-domain pass at 8.5 % of the boundary-layer step; the gate is that
the restricted kernel changes nothing.

```bash
REF=/abs/ref/build_cpu_nofma/moby_solve NEW=/abs/build_cpu_nofma/moby_solve \
    [RANKS=1] ./run_gates.sh
```

Three legs, both binaries nofma, each side preparing its own case file:

1. `trip_small.ini` (a 512 x 176 x 48 boundary layer, trip ON on 7 of 32
   blocks, 20 steps): REF vs NEW, `un vn wn pn` at tolerance 0;
2. its twin with `trip_amp = 0` vs NEW: `vn` MUST differ (the gate is not
   vacuous -- it once ran the GENERIC case instead of the boundary layer,
   which is why a foreign `[case.<x>]` section is a config error now);
3. `force_sine.ini` (a generic `[force]` sine profile on three components,
   128 x 64 x 16, 10 steps): REF vs NEW at tolerance 0 -- the whole-domain
   support of every non-trip source.

`trip_small.ini` keeps plain Jacobi `niter = 6` on purpose: it is a
bit-exactness gate of the force kernel, not a physics case.

RE-VALIDATED 2026-10-21, commit 532ef54 (banner, nofma builds), REF
`~/pass2_ref_binaries` (d821341), workstation CPU, 1 rank, 203 s: ALL PASS --
trip ON max_abs 0 on `un vn wn pn`, the trip-off twin moves `vn`, the sine
profile max_abs 0. Settings changed for the rerun: none.
