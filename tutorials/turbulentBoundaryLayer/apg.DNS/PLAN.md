# apg.DNS — plan (2026-11-16)

**CLOSED 2026-11-20: the run is done; the results are in [`README.md`](README.md).**
This file stays as the design record and the run log (the "Update" sections).

DNS of the `../apg.RANS` layer: the ZPG DNS of `../smooth` with its top outlet
under the pressure table found in RANS. The β target is +0.5 from x = 200 to the
end of the read-out range. Nothing below has run yet; the numbers come from
`../apg.RANS/production_170182.h5`, `../smooth.RANS/production_62365.h5`, the
ZPG DNS statistics (`../smooth/production_stats*.h5`) and its last field
(`../smooth/restart_field.h5`).

## Status (2026-11-16 16:40)

- **Decisions taken (user):** D1 Lz = 52, D2 β_DNS = 0.5, D3 go.
- **Start field:** `start_field.h5` (10.7 GB).
  - Made by `assets/tile_restart.py` in job 5242756.
  - The tool was tested end to end on a small case before use: every tiled
    copy steps forward identically to the untiled run on the global restart
    path. The only difference from a block-table restart is the outlet
    planes, which are not carried and restart from zero-gradient.
- **Smoke test (job 5242764, `smoke/`):** 2500 steps on 8 A100, 0.221 s/step.
  - The top-face pressure equals table × 1.087 to 1e-7.
  - L2_div is 2e-4 just after the switch, against 6e-6 in the ZPG run, and
    falling. Watch it decay in Run A.
- **Run A:** job 5242802 (accelerated, 2 nodes, 24 h, self-chaining to 20
  Nov 11:45), in `/hkfs/work/workspace/scratch/xt8786-apgdns/run`.
- **RANS check** of the scaled table: `rans_check/` (CPU, t 5100 → 7100).
  If it misses, the corrected table goes in between chunks (replace
  `ptop.txt`).

### Update (2026-11-16 20:00)

- **Run A, chunks 1–4** (t 15 500 → 16 500) all ran clean, about 46 min
  each. L2_div was back at 6–9e-6 (the ZPG level) by step 822 000.
- **RANS check** (coarse grid, t 7100, steady):
  - Table × 1.087 gives β_RANS +0.592 on x 200–550, not the linear +0.543.
    β grows faster than the table value because δ*/τ_w grows with it.
  - With the ZPG ratio r, the predicted β_DNS is +0.545 (+0.52 … +0.56).
    That would be 9 % over the target.
- **First DNS window** (t 16 000–16 500, 500–1000 t.u. after the switch,
  `dns_pg.py`):
  - β is +0.48 … +0.49 on x 200–350 and falls to +0.38 at x 550; the mean
    over x 200–550 is +0.461.
  - Re_θ flattens at x > 450, so the outer layer downstream has not adjusted
    yet; one flow-through takes about 760 t.u.
- **Decision:** no table change on this evidence. The two predictions miss
  in opposite directions, and the DNS is still developing. The window at
  t 17 000–17 500 decides, and if a correction is needed, `dns_pg.py
  --correct` sets it.

### Update (2026-11-16 22:40): transient over, no correction

- **Chunks 1–8** (t 15 500 → 17 500) all ran clean. L2_div is 5.4e-6.
- **β** on x 200–550 (`dns_pg.py`):

  | window | mean | min | max |
  |---|---|---|---|
  | t 16 500–17 000 | +0.496 | +0.459 | +0.554 |
  | t 17 000–17 500 | +0.489 | +0.459 | +0.520 |
  | t 16 500–17 500 (1000 t.u.) | +0.492 | +0.460 | +0.536 |

  - The two 500 t.u. windows agree to their noise, so the transient is
    over, about 1000 t.u. after the switch.
  - The mean is within ±0.03 of 0.5, so section 3.4 does not apply: no
    correction, and Run A continues. The statistics window starts at
    t = 16 500.
  - Re_θ at x = 550 is about 1040; c_f(550) = 0.0039.
- **The RANS → DNS ratio under the APG is r ≈ 0.83**, not 0.92: the scaled
  table gives β_RANS 0.592 and β_DNS 0.492. So the ZPG ratio overshoots
  under a pressure gradient. The DNS correction was the right arbiter, and
  the 1.087 scaling happened to land.
- **Transition region** (section 2.4; APG t 16 500–17 500 against ZPG
  `../smooth/production_stats.h5` minus `_t0`):

  | x | 5 | 10 | 15 | 20 | 40 | 60 | 80 | 100 | 160 | 200 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | c_f APG / ZPG | 1.001 | 0.966 | 0.949 | 0.959 | 0.965 | 0.962 | 0.970 | 0.964 | 0.932 | 0.942 |
  | U_e APG / ZPG | 1.000 | 1.000 | 1.000 | 0.999 | 0.997 | 0.993 | 0.987 | 0.978 | | |
  | Re_θ ZPG / APG | 182 / 182 | | | 216 / 214 | 238 / 238 | 264 / 266 | 290 / 297 | 317 / 328 | 394 / 432 | 439 / 502 |

  - The transition has not moved. Re_θ agrees to 1 % up to x = 60, and the
    u′_rms peak (0.184 at x = 20, about 0.14 behind it) to 1–2 % up to
    x = 100.
  - c_f is 3–5 % low from the trip onward, against the 1–2 % the RANS
    predicted at x = 20–60. The deficit starts AT the trip (x ≥ 8), where
    U_e and dp_w/dx equal the ZPG run, and not upstream of it (x ≤ 5:
    identical to 1e-3). So it comes from the trip at Lz = 52, not from the
    APG.
  - Candidate cause: the forcing is (amp/nmodes)(√2 cos φ₀ + Σ 2 cos …).
    nmodes 32 keeps the spectral density of the spanwise modes, and the
    total rms to 1.5 %, but makes the spanwise-uniform mode √2 weaker.
  - Window-to-window noise of c_f at x = 10–20 is up to 1.7 %.
  - By x = 80 the ratio (0.970) is the RANS APG value (0.967). The read-out
    range from x = 200 is downstream of the effect, so no action. A ZPG run
    at Lz = 52 would separate trip from Lz if it is ever needed.

## 0. Decisions needed from the user

| # | question | recommendation |
|---|---|---|
| D1 | spanwise box: keep Lz = 26 (the ZPG/CaNS box) or 52 | **52** (section 1.2): 2 nodes instead of 1, same wall time |
| D2 | β target in the DNS: the RANS table as is (β_DNS ≈ 0.45 by the stored ratio r) or β_DNS = 0.5 | **0.5**: scale the table by 1/r̄ = 1.087, check it in RANS first (section 2.3), then one in-DNS correction if needed (section 3.4) |
| D3 | run on legacy HoreKa A100 under hk-project-exasim before the 20 Nov noon stop | yes (section 3.1) |

## 1. Mesh

### 1.1 Wall and outer resolution: keep the x and y grid

The APG lowers u_τ by 20–25 % downstream of x = 150, so in wall units the
grid is finer than in the ZPG case. Upstream of the ramp nothing changes; the
largest u_τ is the trip overshoot at x ≈ 70, where it is the ZPG value.

The APG DNS values below are estimated as the ZPG DNS value times the
APG/ZPG ratio of the two RANS fields. Grid: Δx = 0.203, Δz = 0.191, Δy_wall
= 0.0105, cans y stretching 2.822.

| x | δ99 ZPG DNS | δ99 APG est. | c_f APG est. | Δx⁺ | Δz⁺ | Δy⁺ wall (cell centre) | Δy⁺ at δ99 | Δy/δ99 at δ99 | Lz/δ99 | Ly/δ99 |
|---|---|---|---|---|---|---|---|---|---|---|
| 70 (DNS c_f peak) | 5.2 | 5.3 | 0.0062 | 5.0 | 4.7 | 0.13 | 2.1 | 0.016 | 4.9 | 19 |
| 200 | 8.3 | 9.2 | 0.0044 | 4.0 | 3.7 | 0.10 | 2.7 | 0.015 | 2.8 | 11 |
| 400 | 12.1 | 14.8 | 0.0038 | 3.4 | 3.2 | 0.09 | 3.6 | 0.014 | 1.8 | 6.8 |
| 550 | 14.5 | 18.4 | 0.0037 | 3.3 | 3.1 | 0.08 | 4.1 | 0.014 | **1.41** | 5.4 |
| 640 | 16.4 | 21.2 | 0.0036 | 3.2 | 3.0 | 0.08 | 4.5 | 0.014 | 1.22 | 4.7 |

- **Wall resolution:** Δx⁺ ≤ 5.1, Δz⁺ ≤ 4.8 and Δy⁺_w ≤ 0.13 everywhere. These
  are the ZPG maxima, all at the trip.
- **Outer resolution:** Δy at the edge is 1.4 % of δ99.
- **Box height:** Ly/δ99 ≥ 4.7 at the outlet. That is ample for the top outlet,
  where V_top is +0.02 … +0.03, i.e. outflow.
- **Box length:** the outlet zone (about 2 δ99, about 45 δ*₀ here) stays inside
  the 100 δ*₀ already excluded, so the read-out range is x ≤ 550 as in the ZPG
  case.
- **Restart:** no re-interpolation is needed. The ZPG DNS field is on the same
  x/y grid.

### 1.2 Spanwise width: Lz = 26 is narrow for the APG outer layer

Lz/δ99 at x = 550 drops from 1.80 (ZPG) to 1.41 (APG).

The spanwise correlation R_uu(Δz) of u in the ZPG DNS was measured from its
last field (`restart_field.h5`), averaged over an x band:

- **x = 530–570:** R_uu(Lz/2) is between −0.24 and +0.03 at y/δ99 = 0.05–0.8.
- **x = 380–400:** R_uu(Lz/2) reaches −0.40 at y/δ99 = 0.4 and −0.66 at 0.8.

This is one snapshot, so it is noisy. Even so, the outer large-scale u
structures already span the box width in the ZPG case. An APG strengthens
exactly those structures: the outer peak of u′u′ is the APG signature. In
this box the APG layer would be narrower in δ99 units than any station of the
ZPG case.

**Recommendation: Lz = 52, nz = 272.**
- Δz stays the same, so the cell count doubles to 334 M.
- On 8 A100 that is 42 M cells per GPU, the same load as the ZPG run on 4.
- nz = 272 is even, as the red-black projection needs.

**The start field** is the ZPG field tiled twice in z. The copy is exact and
divergence-free.
- A 1e-3 random perturbation on one copy breaks the exact Lz/2 symmetry; the
  first projection removes its divergence.
- The copies decorrelate as the trip-born flow and the APG adjustment sweep
  through, about one flow-through. That time falls inside the transient that
  has to be discarded anyway (section 3.3).
- Monitor R_uu(Lz/2) during the transient to confirm it.

**The trip with Lz = 52.** The trip's spanwise signal is a Fourier series in
2π/Lz (`bodyforce.f90`), so doubling Lz halves its wavenumbers. To keep the
ZPG forcing spectrum:
- set `trip_nmodes = 32`, which gives harmonics 1–16 of 2π/52, i.e. the same
  band down to a 3.25 wavelength, with twice the mode density;
- set `trip_amp = 0.18854 · √2 = 0.26664`. The effective amplitude is
  amp/√nmodes, so the rms stays 0.047 and the spectral density per
  wavenumber is unchanged.

Section 3.3 has a check that the transition does not move.

If Lz = 26 is kept (D1), the ZPG field is used directly and no tool is needed.
The trip stays as it is.

## 2. Boundary conditions: RANS → DNS

### 2.1 What carries over unchanged

These are the same in the ZPG DNS and in the RANS:
- **Inlet:** Blasius at Re_δ*,0 = 450.
- **Wall:** no-slip.
- **x = 650:** the predicted outlet face.
- **Top (y = 100):** an outlet.
- **Pressure gradient:** set only through the top outlet's pressure, with the
  same key on the same face:

```ini
[boundary]
y_max_p_profile = table
y_max_p_value   = 1.0        ; or 1.087, see D2 / section 2.3
y_max_p_table   = ptop.txt   ; = ../apg.RANS/ptop.txt (p_top 0 -> -0.179 at the inlet, 0 at x = 650)
```

- **Grid:** the table is interpolated linearly in x, so the finer DNS x grid
  needs nothing extra.
- **Datum:** p_top(650) = 0 is the x outlet's datum, as in the RANS.
- **Device:** the feature is on `main` (8238e75). On the GPU it was exercised
  by every RANS run.
- **Smoke test (section 3.2):** check, on the DNS case, that the top-row
  pressure follows the table and that V_top responds.

### 2.2 What does not carry over

These RANS keys are dropped:
- `[turbulence] model = rans` and the whole `[rans]` section (`tu`,
  `nut_ratio`, `ambient_sustain`, `kpin_box`, both `ktrip_box`);
- `trip_enabled = false`.

The DNS takes the ZPG DNS production block verbatim:
- the trip `trip_amp 0.18854, x0 10, lx 4, ly 1, ts 4, nmodes 16` (or the
  Lz = 52 pair from section 1.2);
- the numerics: `convection = divergence`, red-black with `niter 6` and
  `sor 1.5`, `dtmax 0.02`.

The time step is capped by dtmax, with a Courant number of 0.41 at the cap,
so `cflmax` is inactive. The free stream of the DNS is laminar. The RANS
ambient Tu = 0.1 % was only a model need.

### 2.3 The β the table gives in the DNS

The table was tuned so that the RANS has β = 0.5. The DNS layer differs from
the RANS in the ZPG pair: its δ* is 2–5 % smaller and its τ_w 3–6 % higher.
The same dP_e/dx therefore gives β_DNS = r · β_RANS with r = 0.90–0.94
(`../apg.RANS/assets/beta_ratio_zpg.txt`, mean 0.92), so about 0.45.

r is measured on the ZPG pair and may drift under the APG. For a DNS target
of 0.5 (D2):

1. **Scale the table** by 1/0.92 with `y_max_p_value = 1.087`. This scales g
   exactly and keeps the datum.
2. **Check it in RANS.** Run a coarse RANS (800 × 192, from
   `../apg.RANS/coarse_121694.h5`) with that value until steady, about
   1500–2000 t.u. The thickening feedback makes β superlinear in g, so this
   check is needed before β_RANS · r is trusted.
   - **Accept** if β_RANS · r(x) stays in 0.5 ± 0.02 on x 200–550.
   - **Otherwise** re-run `top_optimise.py` from the it14 coefficients with
     target 0.5/0.92 = 0.543, which should take 3–5 iterates.
   - **Where to run it:** on the CPU build (614 k cells), since the local GPU
     carries the sailplane SD7003 chain until about 17 Nov 03:00.
3. **Correct once in the DNS** (section 3.4) if the measured β_DNS misses.

### 2.4 The transition region

In the RANS the k-source boxes impose the transition, so the pressure field
cannot move it. In the DNS the transition is physical. The pressure gradient
reaches upstream of the ramp start (x = 80), into the DNS's late transitional
stage (c_f overshoot peaking at x ≈ 70).

RANS APG against ZPG:

| x | 20 | 40 | 60 | 70 | 80 | 100 | 120 | 160 | 200 |
|---|---|---|---|---|---|---|---|---|---|
| dp_w/dx APG (ZPG) ·10⁴ | +4.4 (+3.6) | +0.5 (−1.0) | +2.0 (−0.3) | +3.0 (+0.2) | +2.3 (−0.8) | +3.6 (−0.3) | +4.4 (−0.1) | +4.9 (0.0) | +4.4 (0.0) |
| β APG | 0.18 | 0.02 | 0.05 | 0.07 | 0.09 | 0.17 | 0.25 | 0.41 | 0.48 |
| c_f APG / ZPG | 0.993 | 0.988 | 0.981 | 0.976 | 0.967 | 0.946 | 0.926 | 0.889 | 0.868 |

At x ≤ 20, both cases have the same adjustment pressure gradient behind the
trip and the Blasius inlet.

- **What changes.** The forcing at x = 10 sees the same flow to within 1 %.
  The overshoot region at x = 40–80 sees a weak APG: β ≤ 0.09, and the RANS
  c_f is 1–3 % lower there.
- **What the RANS cannot say.** An APG on a non-equilibrium, just-transitioned
  layer lengthens and strengthens the overshoot relaxation. Re_θ at x = 200
  (500 in the RANS) and the β ramp 100–200 may come out differently in the
  DNS.
- **What is fixed.** The trip keys, the ramp position (80–200) and the target
  range (from x = 200) stay as they are. Moving the ramp would mean a new
  RANS optimisation, for a weak effect.
- **What is measured** (section 3.3): DNS APG against DNS ZPG for c_f(x),
  Re_θ(x) and the u′_rms peak over x ≤ 100. They should agree to within
  sampling up to x ≈ 40 and separate by 1–3 % in c_f by x = 80. A larger
  shift means the APG has moved the transition. Then decide at the
  correction step:
  - either include x = 120–200 in the correction fit;
  - or start the read-out range further downstream.

## 3. Running it on HoreKa

### 3.1 Budget and calendar (measured 2026-11-16 15:40 on legacy HoreKa)

**Budget** (`kit_project_usage`, hk-project-exasim):

| resource | granted | used | left |
|---|---|---|---|
| GPU | 146 814 GPU-h | 138 915 | **7 899 GPU-h** |
| CPU | 14.3 M CPU-h | 11.8 M | 2.6 M CPU-h |

The account has 14 157 GPU-h of the user's jobs in the last 30 days.

**Calendar:**
- **Job stop:** legacy jobs end on **20 Nov 2026 at noon**. exasim is a
  legacy-only association; HoreKa 2 has only hk-project-p0024333, and its
  GPUs did not serve jobs as of 11 Nov (`~/HOREKA.md`).
- **Power outage:** 23 Nov to about 4–7 Dec, with no data access. Copy the
  results off by 22 Nov.
- **Workspace:** `roughBoundaryLayer` expires in 5 d 17 h, about 22 Nov.
  Allocate a new one (`ws_allocate apgdns 60`) rather than run inside it.
- **Queue:** `accelerated` had 5 idle nodes, MaxTime 2 days. Two of the
  user's own jobs are running there (not ours, do not touch).

**Cost.** The ZPG DNS ran at 0.18–0.193 s/step on 4 A100. dt = 0.02, so 1
t.u. is 50 steps.

| run | t.u. | Lz = 26, 1 node | Lz = 52, 2 nodes (est. 0.21 s/step) |
|---|---|---|---|
| A: transient + first window (3.3) | 2 500 | 6.7 h, 27 GPU-h | 7.3 h, 58 GPU-h |
| B: after the correction (3.4) | 3 500 | 9.4 h, 38 GPU-h | 10.2 h, 82 GPU-h |
| total | 6 000 | ~16 h, ~65 GPU-h | ~18 h, ~140 GPU-h |

The budget is not the limit; the calendar is. From a submission on 17 Nov
evening to the 20 Nov noon stop, about 60 h of wall time remain.

**Already on HoreKa:** the ZPG start field is in
`.../optimiseBlockRefinement/rough_tbl_run/smooth_775000.h5`. It is the same
step and size as the local `../smooth/restart_field.h5`, so no 5.3 GB upload
is needed.

### 3.2 Before the production job (workstation, 17 Nov)

1. **Inis**, from `../smooth/production_stats.ini` plus sections 1–2:
   - `apg_dns.ini` with the grid, trip, `[boundary]` table and `[restart]`;
   - `stats_sample_interval 25`;
   - `stats_write_interval 12500` (250 t.u.), with a new `stats_file`
     `apg_stats.h5`, so that nothing accumulates onto the ZPG statistics;
   - `ptop.txt` copied from `../apg.RANS`.
2. **Start field**, Lz = 52 only. Write `assets/postpro/tile_restart.py`:
   - assemble the four blocks to global (z, y, x);
   - tile ×2 in z and add the 1e-3 perturbation to one copy;
   - write legacy global-3D (no `blocks` dataset, as `make_global_restart.py`
     does), readable on any rank count;
   - set the attrs nz = 272 and lz = 52 and the z line;
   - drop the outlet planes (CLAUDE.md: an IC tool writes the plane or
     deletes it).

   Check it locally on a small cut, then run it on HoreKa as a `cpuonly` job
   (10.7 GB out).
3. **Case file:** `mpirun -n 8` (or 4) `moby_prepare apg_dns.ini` on the CPU
   build, at the solve rank count; it is about 35 kB. Ship it with `ssh horeka
   "cat > ..."`, because plain scp fails.
4. **Code on HoreKa:** a separate worktree at `4b0dcc7`, so the shared
   `moby-2to1-code` clone is not touched while the user's other jobs build
   from it:

   ```bash
   git -C moby-2to1-code fetch && git -C moby-2to1-code worktree add --detach <W>/moby-apgdns 4b0dcc7
   ```

   Build inside the job, as `submit_outlet_tbl.sh` does.
5. **Smoke test** on `dev_accelerated`, 1 h, about 2 000 steps from the start
   field. Check:
   - the top-row p equals the table;
   - V_top moves toward the RANS V_top;
   - no warnings, and L2(div) at the ZPG level;
   - the s/step and the GPU memory.

### 3.3 Run A: transient and first window

**The job:**
- `accelerated`, account hk-project-exasim, 2 nodes × 4 ranks (or 1 × 4).
- Self-chained 24 h legs with `--dependency=afterok`, whose time limits end
  before 20 Nov 12:00.
- `[restart]` from the start field at t = 15 500. Statistics are on from the
  first step; the leg script copies the cumulative stats file at every write
  to `apg_stats_<step>.h5`.
- Windows are differences of these copies, so the settling is read
  afterwards (the DNS analogue of the RANS chunk marching), and no phase
  boundary has to be guessed in advance.

**Read-out** (workstation, stats copied back, 206 MB each, 412 MB with Lz =
52). Per 250 t.u. window:
- β(x), c_f, H and Re_θ with the Griffin edge on the DNS mean pressure. The
  statistics carry P, accumulator 9. `pg_bl.pg_integrals` needs only a
  loader for the stats layout.
- R_uu(Lz/2) from the latest field, for Lz = 52.

**The transient is over when** two successive windows agree in β(x) to the
window noise on x 200–550. The RANS needed about 1750 t.u. from the ZPG
field to 0.2 % in c_f; expect 1500–2000 t.u. in the DNS.

Then compare the transition region against the ZPG DNS (section 2.4).

### 3.4 Correction and Run B

**If** the first post-transient window (1000 t.u.) misses β = 0.5 by more
than ±0.03 on x 200–550, correct the table once:
- Use the RANS response model of `top_optimise.py`: δβ = κ (δ*/τ_w) R δa on
  the DNS state, with R the Laplace response and κ from the RANS.
- Fit the new coefficients, write the new table, and restart from the last
  DNS field. About 1000 t.u. of re-settling follows (smaller step, faster).
- This needs a small tool (`dns_correct_table.py`) that reuses the
  optimiser's functions on the DNS window.
- A second correction does not fit before 20 Nov.

**Otherwise** Run A simply continues. The statistics window is everything
after the transient, at least 2000 t.u. The ZPG case's 1250 t.u. gave c_f to
±1.4 % per station.

### 3.5 After the stop (20–22 Nov)

- **Copy back:** the stats copies, the last field (5.3 / 10.7 GB) and the
  logs, to the workstation or to the HoreKa 2 workspace from its login node
  (the legacy mounts stay readable there until December).
- **Write-up:** `apg.DNS/README.md` with β(x) DNS against RANS, c_f, H,
  Re_θ, U⁺, the stress profiles at stations in Re_θ against the ZPG DNS, the
  measured r_APG against the stored r_ZPG, and the transition-region
  comparison.
- **Records:** a worklog entry and a row in the tutorial index.

## 4. Order of work

| when | step |
|---|---|
| 16–17 Nov | D1–D3; RANS check of the scaled table (CPU); inis; `tile_restart.py` (if Lz = 52); case file; worktree on HoreKa; smoke test |
| 17 Nov evening | submit Run A (chained legs) |
| 18 Nov | read the windows; the transition comparison; the correction decision → Run B |
| 18–20 Nov noon | statistics |
| 20–22 Nov | copy back before the outage; write-up |

## 5. Risks

- **The queue.** `accelerated` starts were within hours this month (the
  sailplane calib job started a day before its estimate), but they are not
  guaranteed. If Run A starts after about 19 Nov, there is no correction step;
  run the table chosen in D2 to the stop.
- **No GPU after 20 Nov.** The next GPU time is HoreKa 2 (p0024333, H100,
  not yet serving) or local. A 334 M-cell DNS does not fit the RTX 3060.
- **The tiled start field** carries the exact Lz/2 symmetry until the trip
  and the perturbation break it. R_uu(Lz/2) shows when; windows before that
  are discarded.
- **The top outlet under a DNS free stream** is unchanged from the ZPG run
  (same face, same predictor); only its datum varies along x.
