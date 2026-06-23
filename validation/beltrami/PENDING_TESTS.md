# §6e both-worlds (MOBY_AGGLOM_MPFEED) — results so far + PENDING tests

Branch `claude/blocks`. Session 2026-06-23. The §6e "both-worlds" combo is
implemented behind **`MOBY_AGGLOM_MPFEED`** (`lin` = stock (3/4,1/4) linear
prolong cheap pre-test; `mp`/`1` = mean-preserving Piece-1 prolong). It refills
the slaved high halo with a tangentially-varying prolong ONCE, before
`composite_qs_setup` freezes the predictor `qs` — predictor side only, never the
in-projection exchange loop. Independent of `MOBY_AGGLOM`, so it composes with
per-cell / flat / additive pressure.

## RESULTS ALREADY OBTAINED (keep — do not rerun)

Bit-exact baselines (CPU, flag unset) — reproduce the references exactly:
- off Test B refined = 1.0340e-3 (single block 0.0)
- flat Test B refined = 4.1607e-4 (single block 0.0)

Test B (corrector-alone) — INSENSITIVE to MP feed (the predictor-base shift
contaminates ||900003-900002||, so B is NOT a useful metric for this combo):
- per-cell+MP 1.0445e-3, flat+LIN 4.1694e-4, flat+MP 4.1783e-4, additive+MP 3.7938e-4

DECISIVE metric — agglom_trajectory.py gpu niter=10 --interval 200, peak ratio:
| config | peak | note |
|---|---|---|
| off (per-cell, production default) | 1.03788 | = reference 1.038 |
| flat-only agglom | 1.04328 | production agglom base |
| **per-cell + MP** | **1.02859** | WINNER — below off AND flat |
| flat + MP | 1.04116 | only marginally below flat |
| additive edge + MP | 1.05429 | WORSE (edge pressure agglom regresses) |
| flat + LIN (cheap pre-test) | 1.04167 | ~ flat -> linear prolong barely moves it; the MEAN-PRESERVING P_mp is what matters |

interface_decay gate (CPU, 200 steps) — ALL PASS (u,v,w,p contract):
off, per-cell+MP, flat+MP, additive+MP, flat+LIN.

CONCLUSION (preliminary): the both-worlds win is **per-cell pressure +
MP-fed predictor velocity** (1.0286), NOT agglomerated pressure. The MP-fed
predictor alone drops the edge trajectory below flat-only; adding pressure
agglomeration on top HURTS (the shared pressure fights the per-face MP velocity,
exactly the §6g conflict). The cheap linear-prolong feed is insufficient.

## 6-CASE SWEEP (32^3 refined_fast, both Re, L2-vs-exact metric) — 2026-06-23

| case | Re | peak L2 vs exact | peak ratio | status |
|---|---|---|---|---|
| off (per-cell)  | 100 | 0.500 | 1.187 | stable |
| flat agglom     | 100 | 0.356 | 1.305 | stable |
| per-cell+MP     | 100 | 0.332 | 1.087 | stable (BEST) |
| off (per-cell)  | 400 | 0.215 | 1.233 | stable |
| flat agglom     | 400 | --    | --    | **BLEW UP t~2.4 (max|u| 1e116, dt->1e-119)** |
| per-cell+MP     | 400 | 0.210 | 1.282 | stable |

DECISIVE: **flat pressure agglomeration DIVERGES at Re=400** (stable at Re=100,
the validated production fix — but NOT robust). per-cell + MP-fed predictor
(MOBY_AGGLOM_MPFEED=mp, NO agglomeration) is stable at BOTH Re and most accurate
(best L2 at Re=100; tied-best at Re=400 though higher Linf). => the both-worlds
win is the MP-fed predictor velocity with per-cell pressure; agglomeration is not
just unnecessary but HARMFUL at higher Re. (peak-ratio and L2 disagree on
off-vs-flat at Re=100; trust L2. At Re=400 pcmp ~ off on L2, pcmp Linf higher.)
Field-error plots: scratchpad/plots/{r100,r400}_{off,flat,pcmp}.png.
TODO: confirm flat-Re400 blow-up is not niter-starvation (try niter=30); check
per-cell+MP long-time (t=40) stability at both Re.

## PENDING TESTS (stopped 2026-06-23 mid-run to avoid resource thrash)

1. **Whole-Beltrami Test C (t=8)**, the accuracy/stability confirmation, GPU:
   - off (baseline, expect ~3.7e-2 — was mid-run, NOT captured)
   - flat (MOBY_AGGLOM=1)
   - per-cell+MP (MOBY_AGGLOM_MPFEED=mp) — expect a DROP toward the single-block
     1.36e-4 if the MP-fed velocity is genuinely 2nd order (the key accuracy claim)
   - flat+MP (MOBY_AGGLOM=1 MOBY_AGGLOM_MPFEED=mp)
   Command: `python3 tools/beltrami_regression.py gpu --tests C --tag <t> [env]`
2. **GPU off ABC bit-exact** confirmation (flag unset must match the prior
   baseline on GPU as well as CPU). Was mid-run, NOT captured.
3. **GPU interface_decay** for the variants (only the CPU decay gate was run).
4. **Longer-time stability** of per-cell+MP (t=20 trajectory is bounded; a
   t=40 whole-Beltrami would confirm no slow blow-up — the original concern the
   agglomeration was meant to fix; per-cell+MP passes the 200-step decay gate but
   a long run is the stronger check).
5. **Fast-case validation** (see refined_fast.ini below): confirm the 32^3
   case reproduces the per-cell+MP < off < flat ranking, and find the Re/tfinal
   that compresses the trajectory without inverting the (late-time) ranking.

## Fast validation case (engineered + VALIDATED this session)

`validation/beltrami/refined_fast.ini` — 32^3, nb=4, central 2x2x2-block patch
(interfaces+edges+corners all dirs). ~8x fewer cells + ~2x CFL dt than refined.ini.
Tool: `tools/agglom_trajectory.py --src refined_fast.ini --niter 10 --interval 100`.

VALIDATED (GPU, niter=10, t=20, ~1 min/run vs ~20 min for 64^3): reproduces the
ranking EXACTLY with a ~20x larger spread —
  per-cell+MP 1.086  <  off 1.187  <  flat 1.305    (cf 64^3: 1.029 < 1.038 < 1.043)
All bounded, no NaN. Use re=100, t_final=20 (the late-rise carries the signal).

Re LEVER = DEAD END (do not retry): re=400/t=10 (~29 s) INVERTS the ranking
(per-cell+MP 1.282 > off 1.233) — high re kills the viscous decay the late-time
discrimination needs, and the early transient (where per-cell+MP's floor is
higher) dominates. Only fewer cells is a safe speed lever, not Re/time.

Cheapest mechanism check for the MP prolong alone (does NOT rank per-cell vs
flat): `beltrami_regression --tests A --mp` (div-of-exact 1.44 -> ~0.02 in one
substep) + interface_decay (~1 min). The trajectory is still needed to rank the
pressure variants.
