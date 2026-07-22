#!/usr/bin/env python3
"""Term attribution of the v1 interface instability from diag_terms.bin
(the diag/v1-terms worktree dump: per substage, tag 1 = before the
momentum predictor, 2 = after predictor + halo sync, 3 = after the
projection correction; q_u and q_w with one halo layer per block).

For each dumped step and substage:
  pred = q2 - q1, corr = q3 - q2          (raw increments)
  p'   = q1 - q1(base step)               (the perturbation)
  e    = p'/|p'|                          (mode direction, fluid cells)
  G_pred = e . (pred - pred_base)         (perturbation-projected growth)
  G_corr = e . (corr - corr_base)
  D      = e . nu Lap(p') dt_sub          (exact discrete diffusion of the
                                           perturbation, dumped halos)
  G_conv = G_pred - D                     (the predictor remainder:
                                           convection + IBM, IBM = 0 in
                                           fluid cells)
Background subtraction (pred_base/corr_base from the first dumped step)
removes the quasi-steady base-flow increments.
"""
import numpy as np

NB = 10          # 0:nb+1 per dim
NI = 512         # interior cells per block (8^3)
DT = 5.0e-5
NU = 1.0/4.0e5
RK_GAMMA = np.array([8.0/15.0, 5.0/12.0, 3.0/4.0])   # dt_sub per stage


def read_all(path):
    raw = open(path, "rb")
    recs = {}
    order = []
    while True:
        hdr = np.fromfile(raw, dtype=np.int32, count=4)
        if hdr.size < 4:
            break
        step, stage, tag, nd = map(int, hdr)
        blocks = []
        for _ in range(nd):
            meta = np.fromfile(raw, dtype=np.int32, count=5)
            qu = np.fromfile(raw, dtype=np.float64, count=NB**3).reshape(
                NB, NB, NB, order="F")
            qw = np.fromfile(raw, dtype=np.float64, count=NB**3).reshape(
                NB, NB, NB, order="F")
            blocks.append((tuple(meta), qu, qw))
        recs[(step, stage, tag)] = blocks
        order.append((step, stage, tag))
    return recs, order


def interior(blocks):
    """Concatenated interior u,w over all blocks."""
    return np.concatenate([np.concatenate(
        [b[1][1:9, 1:9, 1:9].ravel(), b[2][1:9, 1:9, 1:9].ravel()])
        for b in blocks])


def lap(blocks, levels_h):
    """nu Lap(field) on interior cells using the dumped halo layer,
    per block at its own level spacing (x and z; y uniform coarse and
    the mode is y-uniform after span-averaged forcing — include it)."""
    out = []
    for meta, qu, qw in blocks:
        lev = meta[1]
        hx, hz = levels_h[lev]
        hy = 0.1875/8.0
        for q in (qu, qw):
            c = q[1:9, 1:9, 1:9]
            d = ((q[2:10, 1:9, 1:9] - 2*c + q[0:8, 1:9, 1:9])/hx**2
                 + (q[1:9, 2:10, 1:9] - 2*c + q[1:9, 0:8, 1:9])/hy**2
                 + (q[1:9, 1:9, 2:10] - 2*c + q[1:9, 1:9, 0:8])/hz**2)
            out.append(NU*d.ravel())
    return np.concatenate(out)


def main():
    levels_h = {l: (128.0/(384*2**l), 96.0/(288*2**l)) for l in range(12)}
    recs, order = read_all("diag_terms.bin")
    steps = sorted({s for s, _, _ in recs})
    base_step = steps[0]
    print(f"{len(steps)} dumped steps: {base_step} .. {steps[-1]}")

    base_q1 = {st: interior(recs[(base_step, st, 1)]) for st in (1, 2, 3)}
    base_pred = {st: interior(recs[(base_step, st, 2)]) - base_q1[st]
                 for st in (1, 2, 3)}
    base_corr = {st: interior(recs[(base_step, st, 3)])
                 - interior(recs[(base_step, st, 2)]) for st in (1, 2, 3)}
    # fluid mask from the base state (penalized cells ~ 1e-26)
    fluid = np.abs(base_q1[1]) > 1e-15

    # perturbation-vs-base BLOCK subtraction needs matching blocks — the
    # dump order is fixed (same diagIds every time), so vectors align.
    print(f"{'step':>7} {'stage':>5} {'|p1|':>10} {'G_pred':>11} "
          f"{'G_corr':>11} {'D_diff':>11} {'G_conv=P-D':>11}")
    for step in steps[1:]:
        for st in (1, 2, 3):
            q1 = interior(recs[(step, st, 1)])
            q2 = interior(recs[(step, st, 2)])
            q3 = interior(recs[(step, st, 3)])
            p1 = (q1 - base_q1[st])[fluid]
            nrm = np.linalg.norm(p1)
            if nrm < 1e-14:
                continue
            e = p1/nrm
            gp = float(e @ ((q2 - q1)[fluid] - base_pred[st][fluid]))
            gc = float(e @ ((q3 - q2)[fluid] - base_corr[st][fluid]))
            # diffusion of the perturbation (halo-inclusive)
            pblocks = [((m1), qu1 - qu0, qw1 - qw0)
                       for (m1, qu1, qw1), (m0, qu0, qw0)
                       in zip(recs[(step, st, 1)], recs[(base_step, st, 1)])]
            dl = lap(pblocks, levels_h)[fluid]
            dd = float(e @ dl)*DT*RK_GAMMA[st-1]
            if step % 200 == 0 or step > steps[-1] - 40:
                print(f"{step:>7} {st:>5} {nrm:>10.3e} {gp:>11.3e} "
                      f"{gc:>11.3e} {dd:>11.3e} {gp - dd:>11.3e}")


if __name__ == "__main__":
    main()
