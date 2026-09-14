# Pre-registered readings — the between-iteration divergence halo, with peers

Written **before** the job ran, at `3c2903a`. Task 1 of
`docs/next_session_divergence_halos.md`.

## What is being changed

Between projection iterations the solver currently refreshes the velocity halos
with `exchange_halos(c, blk, [U,V,W], interp=.false.)` — the whole 26-direction
same-level shell, three components. At **one rank** it instead calls
`sync_divergence_halos`, which writes only `q(nb+1)` per dimension of the
component NORMAL to that face. That routine stopped at one rank because with
peers the same planes have to come over MPI and the entry list was not
partitioned for it.

Now it is. The entry enumeration gains a **third round**: pure `+axis` face
same-level COPY entries are emitted FIRST, then the remaining same-level copies,
then the cross-level entries. A divergence round is therefore a prefix of the
copy prefix, on both ends of every message, and carries **one variable per
entry** (the face normal) instead of three. The same prefix drives the
same-rank kernel, so `dsSlot` and its three little kernels are deleted: one
mechanism, one local kernel, peers or not.

## Why the reduced set is complete — by construction, not by gate

The claim to check before writing code was that the reduced set delivers every
halo value the next divergence reads. It does, and the entry list says so:

1. **What is read.** `jacobi_compute_phi` forms the divergence over `i,j,k =
   1..nb`, so the only halo cells it touches are `q(nb(d)+1, ·, ·, VAR_d)` at
   *interior* tangential indices `1..nb`. `jacobi_apply` reads `phi`, not the
   velocity halo; `interface_correct` is documented interior-tangential-only;
   `apply_bc` works on physical ghosts the block owns.
2. **What writes it.** From `entry_boxes`: an entry whose direction has
   `off(d) = +1` has `dstLo(d) = nb(d)+1, ext(d) = 1`; one with `off(d) = -1`
   writes `{0}`; one with `off(d) = 0` writes `lo..hi` which reaches `nb(d)+1`
   only through the tangential extension, i.e. only at *halo* indices of the
   other dims. So within the read range, the plane is written by **pure `+axis`
   face entries and nothing else** — edges and corners land at tangential index
   0 or nb+1, outside it.
3. **Cross-level ones write nothing there anyway.** For a 2:1 `+axis` face
   entry `interface_normal_dim` returns `d`, and `unpack_entries` /
   `copy_local_cross_level` skip `var == rNrm`: the low-side block owns that
   face and reconstructs it. So excluding them is not an approximation.

Hence: pure `+axis` same-level face COPY entries = exactly the writers of the
divergence halo. Everything else the full copy round writes is either read no
earlier than the next substage — after the last iteration's FULL exchange
restores it — or is written with the identical value.

The empirical corroboration is that **`1 rank == 4 ranks` is a standing gate
that passes today**, with the reduced set in use on one side and the full shell
on the other.

## The bucket being targeted (job 5145554, `after` column — the current code)

| case | `proj vel_exchange` | of which between-iteration | step |
|---|---|---|---|
| `rect_jacobi` 4×1 | 12.321 ms (7.8 %) | 15 of 18 calls | 158.75 ms |
| `rect_jacobi` 8×2 | 7.409 (8.8 %) | 15 of 18 | 83.88 |
| `refined_yp82` 4×1 | 6.243 (8.3 %) | 15 of 18 | 75.00 |
| `refined_yp82` 8×2 | 4.353 (10.2 %) | 15 of 18 | 42.65 |

Per block (`nb = 64 44 48`) a full same-level copy round moves 16632 points × 3
components = **49896 values**; the `+axis` faces at one component each are
2112 + 3072 + 2816 = **8000**. A **6.24x** cut in bytes, with the launch count
unchanged at three kernels per round (pack, local copy, unpack).

## What the job must return

| result | conclusion |
|---|---|
| `proj vel_exchange` **−40 to −60 %**, step **−3.5 to −5.5 %** on both cases at both rank counts | as designed; take it |
| the new `exchange MB per round` line shows `div nv=1` at **1/6.0 to 1/6.5** of `copy-only nv=3` | the partition selected the intended entries. If it is much larger the predicate is too loose; if much smaller, entries are missing and a gate should have failed |
| `proj vel_exchange` falls but the step does not | the saving was in a bucket that overlapped something else; report the bucket and do not quote a step win |
| step falls by materially MORE than ~6 % | distrust it and find what else changed — the launch count is unchanged, so there is no second mechanism available |
| bucket barely moves | the div prefix is nearly the whole copy prefix (wrong predicate), or the round is falling back to the full exchange. Check the MB line first |
| any nonzero `max_abs` | a real difference. The 1-rank path changes too (one kernel, tangentially extended, instead of three interior-only ones), so the 1==4 gate and Pass G both bear on it |

**Gate flags: PRODUCTION, not nofma.** No expression moves: the values written
are copied, never recomputed, and every surviving kernel's arithmetic is
textually what it was. By the rule in `run_mapgate.sh` the production comparison
is then the tighter test. If it fails at 1e-16 the rule was wrong here and
`submit_nofma_gate.sh` is the fallback — worth recording either way.

## What this job cannot say

- **Nothing about 16 ranks** (`accelerated` queue). The per-round message is
  smaller, so `mpi_wait` should benefit more where the network matters more,
  but that is an extrapolation and must be labelled one.
- **Nothing about the launch floor.** Three kernels per round before, three
  after; the saving is bytes and wait, and if the measured saving looks like a
  launch saving it is being misread.
- **Nothing about red-black.** It uses per-colour copy-only exchanges, not this
  path, and is untouched.
