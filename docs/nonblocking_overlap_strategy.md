# Nonblocking Halo-Exchange Overlap Strategy

These notes collect the proposed communication/computation overlap strategy for a
possible later optimization. The idea should be revisited after profiling shows
that halo exchange is a relevant cost.

## Motivation

The current MPI path already has a nonblocking halo API:

```fortran
call start_halo_exchange(c, f, vars)
call finish_halo_exchange(c, f)
```

At the moment the code mostly uses the blocking wrapper:

```fortran
call exchange_halos(c, f, vars)
```

This is simpler and should remain the default until profiling shows that MPI
communication is a bottleneck. If communication becomes important, the main
opportunity is the red-black pressure solver, because it exchanges velocity
halos after every color sweep.

## Preferred Refactor

Use a core/shell split, but keep the numerical stencil formula in one place.

- Core points: local points whose stencil does not touch MPI halos.
- Shell points: the one-cell layer close to MPI boundaries.

The core can be computed while halo messages are in flight. The shell is
computed after `finish_halo_exchange`, when the received halos are available.

The code should avoid scattering processor-boundary tests inside the stencil.
Instead, precompute the loop regions once:

```fortran
type stencil_work_type
    integer(C_INT) :: coreLo(3), coreHi(3)
    integer(C_INT) :: nShell = 0_C_INT
    integer(C_INT), allocatable :: shellI(:), shellJ(:), shellK(:)
end type stencil_work_type
```

Then the pressure solver can expose two loop shapes:

```fortran
call redblack_sweep_core(...)
call redblack_sweep_shell(...)
```

Both should use the same point-update logic, so the red-black formula is not
duplicated.

## Pressure Solver Sketch

For each red/black color:

```fortran
call start_halo_exchange(c, f, [VAR_U, VAR_V, VAR_W])

call redblack_sweep_core(...)

call finish_halo_exchange(c, f)

call redblack_sweep_shell(...)
call apply_bc(f, dns, g, bc)
```

For a sequence of red-black sweeps, the next exchange can be started as soon as
the shell and boundary-condition updates for the current color are complete.

## Momentum Solver

The same core/shell machinery could also be used for the momentum predictor:

```fortran
call start_halo_exchange(c, f, [VAR_U, VAR_V, VAR_W])
call momentum_core(...)
call finish_halo_exchange(c, f)
call momentum_shell(...)
```

However, this should not be the first target unless profiling shows momentum
halo exchange is important. The red-black pressure loop is much more likely to
benefit because halo exchange happens repeatedly inside the SOR iterations.

## Deferred Optimization: Color-Filtered Exchange

In red-black SOR, after a color sweep only one checkerboard color has changed.
In principle, only the updated color values that can be read by the opposite
color need to be transferred.

This should be treated as a second-stage optimization:

```fortran
call start_halo_exchange(c, f, [VAR_U, VAR_V, VAR_W], color=color)
```

It would require color-aware packing lists and extra tag/state handling in
`comm.f90`. That is worthwhile only if profiling shows the message volume is a
real problem after the simpler core/shell overlap is implemented.

## Suggested Decision Point

Before implementing this refactor, profile a representative multi-rank case and
check:

- time spent in `MPI_Waitall` or equivalent synchronization,
- time spent in host/device updates around halo boxes,
- red-black sweep time versus communication time,
- whether communication time grows with rank count enough to justify overlap.

If MPI wait time is small compared with GPU kernel time, keep the current
blocking `exchange_halos` calls for readability.
