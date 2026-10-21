! CONJUGATE heat transfer at the immersed interface (increments C1-C3,
! split out of scalar.f90 on 2026-10-21): the one face coefficient the
! scheme consists of, the banded solid (F2), the cut-cell volume fraction
! (C3), the tangential term and its indicator (C2), the Peclet rate, and
! the conjugate initialisation. The transport kernel (scalar.f90) and the
! statistics (scalar_stats.f90) share these definitions.
module scalar_conjugate
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P, VAR_S0, SOLID_FACE_THRESHOLD
    use :: blocks, only: block_set_type, FACE_COARSE, FACE_FINE
    use :: boundary, only: boundary_type
    use :: comm, only: comm_type, comm_allreduce_sum, comm_allreduce_max
    use :: ibmm, only: ibm_type
    use :: io, only: read_dwall_blocks
    use :: scalar_types, only: SC_IBM_CONJUGATE, scalar_conjugate_enabled, scalar_type
    use :: scalar_wallfn, only: eddy_diffusivity, plain_face_diffusivity
    implicit none
    private
    ! Grazing-arm guard (strategy doc Section 8, LaTeX note Section 6.6). The
    ! level-set weight w is exact for a PLANE; for a curved interface the two
    ! cell centres measure to different nearest points and w carries a
    ! relative error O(curvature*h/a), with a = cos(angle between the arm and
    ! the interface normal). Arms nearly TANGENT to the surface therefore
    ! carry an unreliable w -- and little of the interface flux -- so below
    ! this direction cosine the face falls back to the plain harmonic mean.
    real(C_DOUBLE), parameter :: CONJ_MIN_COSINE = 5.0d-2
    ! Smallest magnitude a SOLID cell's signed distance is allowed to take,
    ! so that `phi < 0` is an EXACT material test even for a cell centre
    ! sitting exactly on the surface (see init_scalar_conjugate).
    real(C_DOUBLE), parameter :: CONJ_MIN_DISTANCE = 1.0d-300
    ! Smallest |grad phi| at which the interface normal is still trusted
    ! (increment C2). phi is a distance function, so |grad phi| = 1 wherever
    ! the nearest surface point is unique; it collapses on a MEDIAL AXIS,
    ! where the central difference straddles two different nearest points and
    ! the normal it returns is meaningless. Below this the tangential
    ! correction and its indicator switch off at that face, which returns the
    ! C1 baseline there rather than a spurious flux.
    real(C_DOUBLE), parameter :: CONJ_MIN_GRADPHI = 1.0d-2

    public :: band_face_diffusivity, band_face_total, band_material, check_conjugate_refinement, conjugate_face_diffusivity
    public :: conjugate_face_local_k, conjugate_face_total, conjugate_tangential, init_scalar_conjugate, init_scalar_solid_fields
    public :: plane_box_fraction, scalar_conjugate_indicator, scalar_conjugate_peclet_rate, scalar_conjugate_to_device, scalar_volume_fraction

contains

    ! Face diffusivity of the conjugate interface scheme -- the ONE new
    ! coefficient the whole increment consists of (strategy doc Section 1,
    ! LaTeX note Eq. (6.baseline)). phiL/phiR are the SIGNED distances at the
    ! two cell centres (positive in the fluid), invd the inverse
    ! centre-to-centre distance, dm = 1/(Re Pr) the FLUID molecular
    ! diffusivity and ks = kappa_s the solid-to-fluid conductivity ratio.
    !
    !   uncut  ->  kappa*dm: one material across the face, i.e. today's
    !              kernel line with kappa = 1 in the fluid;
    !   cut    ->  the DISTANCE-WEIGHTED HARMONIC MEAN dm/(w/k_L + (1-w)/k_R)
    !              on the level-set fraction w = phiL/(phiL - phiR),
    !              i.e. two resistances in series.
    !
    ! THE OBLIQUITY LEMMA is what makes w usable at any interface
    ! orientation: phi stores the PERPENDICULAR distance to the surface, not
    ! the distance along the arm, so both arms are shortened by the same
    ! direction cosine a = n.e_d and it CANCELS in the ratio. w is then
    ! exactly the fraction delta_L/h_d that the series resistance needs, with
    ! no normal ever computed and no per-arm cut record ever stored -- which
    ! is why this increment adds no dataset. Exact for a plane at any angle.
    !
    ! The cut test is the MATERIAL DISAGREEMENT of the two centres rather
    ! than the strategy doc's phiL*phiR < 0. They are the same test (the sign
    ! of phi IS the marker), but the marker form cannot be defeated by a cell
    ! centre lying exactly on the surface -- init_scalar_conjugate keeps a
    ! solid cell's phi strictly negative for precisely that reason.
    !
    ! A contact resistance rc (per unit area, same non-dimensionalisation)
    ! adds to the series: R = h_d(w/k_L + (1-w)/k_R)/dm + rc, hence the
    ! rc*dm*invd term. rc = 0 gives exactly the harmonic mean.
    real(C_DOUBLE) function conjugate_face_diffusivity(dm, phiL, phiR, ks, rc, invd) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, phiL, phiR, ks, rc, invd

        real(C_DOUBLE) :: w, kl, kr
        logical :: sl, sr

        sl = phiL < 0.0d0
        sr = phiR < 0.0d0
        kl = merge(ks, 1.0d0, sl)
        if (sl .eqv. sr) then
            d = dm*kl
            return
        end if
        kr = merge(ks, 1.0d0, sr)
        w = conjugate_face_weight(phiL, phiR, invd)
        d = dm/(w/kl + (1.0d0 - w)/kr + rc*dm*invd)
    end function conjugate_face_diffusivity

    ! The level-set fraction of a cut arm, with the grazing guard: a =
    ! |phiL - phiR|/h_d is the direction cosine (a free by-product of the
    ! obliquity lemma) and w is unreliable for arms nearly TANGENT to the
    ! surface, which carry little of the interface flux anyway.
    real(C_DOUBLE) function conjugate_face_weight(phiL, phiR, invd) result(w)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: phiL, phiR, invd

        if (abs(phiL - phiR)*invd < CONJ_MIN_COSINE) then
            w = 0.5d0
        else
            w = phiL/(phiL - phiR)
        end if
    end function conjugate_face_weight

    ! Diffusivity of the material the FACE MIDPOINT lies in -- k_loc of the
    ! exact cut-face flux (LaTeX note Eq. (main)), needed only by the C2
    ! tangential correction. The midpoint is at arm fraction 1/2, so it is on
    ! the low side iff w > 1/2, and with w = phiL/(phiL - phiR) that test is
    ! simply the sign of phiL + phiR -- no division, and it degrades
    ! gracefully where the grazing guard has already replaced w by 1/2.
    real(C_DOUBLE) function conjugate_face_local_k(dm, phiL, phiR, ks) result(k)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, phiL, phiR, ks

        k = dm*merge(ks, 1.0d0, phiL + phiR < 0.0d0)
    end function conjugate_face_local_k

    ! Which material a cell centre is in: 0 fluid, 1 shell, 2 outer band.
    ! dsh <= 0 disables the split, so every solid cell reads 1 and the
    ! arithmetic below collapses onto conjugate_face_diffusivity's.
    integer function band_material(ph, dsh) result(m)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: ph, dsh

        if (ph >= 0.0d0) then
            m = 0
        else if (dsh <= 0.0d0 .or. ph > -dsh) then
            m = 1
        else
            m = 2
        end if
    end function band_material

    real(C_DOUBLE) function band_kappa(m, ks, ko) result(k)
        !$omp declare target
        integer, intent(in) :: m
        real(C_DOUBLE), intent(in) :: ks, ko

        if (m == 0) then
            k = 1.0d0
        else if (m == 1) then
            k = ks
        else
            k = ko
        end if
    end function band_kappa

    ! conjugate_face_diffusivity generalized to the three bands. Same face
    ! coefficient, same level-set weight, evaluated on whichever iso-surface
    ! the face crosses: phi = 0 between fluid and shell, phi = -dsh between
    ! shell and outer. The contact resistance belongs to the BODY SURFACE, so
    ! it enters only at the former.
    !
    ! A fluid|outer face means the band is thinner than the arm, which
    ! init_scalar_conjugate rejects; if one ever reached here it falls into
    ! the phi = 0 branch and reads the outer kappa, i.e. an insulated wall.
    real(C_DOUBLE) function band_face_diffusivity(dm, phiL, phiR, ks, ko, dsh, rc, invd) &
            result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, phiL, phiR, ks, ko, dsh, rc, invd

        integer :: ml, mr
        real(C_DOUBLE) :: w, kl, kr, lev, rcf

        ml = band_material(phiL, dsh)
        mr = band_material(phiR, dsh)
        kl = band_kappa(ml, ks, ko)
        if (ml == mr) then
            d = dm*kl
            return
        end if
        kr = band_kappa(mr, ks, ko)
        ! An insulating partner makes the series resistance infinite. Said
        ! as a branch, because the harmonic mean would divide by zero.
        if (kl <= 0.0d0 .or. kr <= 0.0d0) then
            d = 0.0d0
            return
        end if
        lev = merge(0.0d0, -dsh, min(ml, mr) == 0)
        rcf = merge(rc, 0.0d0, min(ml, mr) == 0)
        w = conjugate_face_weight(phiL - lev, phiR - lev, invd)
        d = dm/(w/kl + (1.0d0 - w)/kr + rcf*dm*invd)
    end function band_face_diffusivity

    ! FLUID VOLUME FRACTION of a box cut by a plane -- the 3D sibling of the
    ! face-area fraction (increment C3, strategy doc Section 6).
    !
    ! phic is the signed distance at the box CENTRE (positive in the fluid)
    ! and a_i = |n_i| h_i the projections of the box edges on the interface
    ! normal. Measuring from the corner the plane reaches first,
    ! s = phic + (a1+a2+a3)/2, the fluid volume is the volume of the box
    ! {0 <= y_i <= a_i} below the plane sum(y_i) = s, which is the standard
    ! inclusion-exclusion over the corner simplices:
    !
    !   V(s) = [ s+^3 - sum_i (s-a_i)+^3 + sum_i<j (s-a_i-a_j)+^3
    !            - (s-a1-a2-a3)+^3 ] / 6 ,      x+ = max(x, 0)
    !
    ! divided by a1 a2 a3. It clips itself: a box the plane misses has
    ! s <= 0 or s >= sum(a), giving exactly 0 or 1, so no interface test is
    ! needed around it.
    !
    ! DEGENERATE DIRECTIONS ARE NOT A SPECIAL CASE, THEY ARE THE COMMON CASE.
    ! A grid-aligned wall has two of the a_i equal to zero (up to the float
    ! noise of the distance field), and the 3D form then divides a numerator
    ! that has cancelled to nothing by a denominator that is nothing. The
    ! limit a_3 -> 0 IS the 2D form and a_2 -> 0 the 1D form
    ! (f = clip(1/2 + phic/h), the strategy doc's simple estimate -- which is
    ! therefore not an approximation for a grid-aligned interface but the
    ! closed form itself), so the routine drops directions whose a_i is
    ! negligible against the largest and evaluates the reduced form. The
    ! threshold is a pure conditioning choice: the cancellation in the full
    ! form costs about eps*(sum a)/a_min relative, and dropping a direction
    ! costs about a_i/sum(a) -- at 1e-8 both are ~1e-8, far below the O(h)
    ! error of the normal itself.
    real(C_DOUBLE) function plane_box_fraction(phic, a1, a2, a3) result(f)
        real(C_DOUBLE), intent(in) :: phic, a1, a2, a3

        real(C_DOUBLE), parameter :: DEGENERATE = 1.0d-8
        real(C_DOUBLE) :: a(3), s, tot, amax
        integer :: m

        amax = max(a1, max(a2, a3))
        m = 0
        if (a1 > DEGENERATE*amax) then; m = m + 1; a(m) = a1; end if
        if (a2 > DEGENERATE*amax) then; m = m + 1; a(m) = a2; end if
        if (a3 > DEGENERATE*amax) then; m = m + 1; a(m) = a3; end if
        if (m == 0) then
            ! The normal has no component along any box edge: the plane is
            ! parallel to the box and the whole cell is on one side.
            f = merge(1.0d0, 0.0d0, phic > 0.0d0)
            return
        end if

        tot = sum(a(1:m))
        s = phic + 0.5d0*tot
        ! Clip FIRST, and exactly. Outside the band the inclusion-exclusion
        ! sum is an identity that holds analytically and cancels numerically:
        ! its terms grow like s^3 while the answer stays 6 a1 a2 a3, so a cell
        ! a few h from the interface would come back as 1 - 1e-13 instead of
        ! 1. Nothing downstream would break, but "which cells are cut" is a
        ! classification, and a classification must not be decided by
        ! round-off (the 2026-08-05 body-heat lesson, one level down).
        if (s <= 0.0d0) then
            f = 0.0d0
            return
        else if (s >= tot) then
            f = 1.0d0
            return
        end if
        select case (m)
        case (1)
            f = (pos(s) - pos(s - a(1)))/a(1)
        case (2)
            f = (pos(s)**2 - pos(s - a(1))**2 - pos(s - a(2))**2 &
                 + pos(s - tot)**2)/(2.0d0*a(1)*a(2))
        case default
            f = (pos(s)**3 &
                 - pos(s - a(1))**3 - pos(s - a(2))**3 - pos(s - a(3))**3 &
                 + pos(s - a(1) - a(2))**3 + pos(s - a(1) - a(3))**3 &
                 + pos(s - a(2) - a(3))**3 &
                 - pos(s - tot)**3)/(6.0d0*a(1)*a(2)*a(3))
        end select
        f = min(max(f, 0.0d0), 1.0d0)
    end function plane_box_fraction

    real(C_DOUBLE) function pos(x)
        real(C_DOUBLE), intent(in) :: x

        pos = max(x, 0.0d0)
    end function pos

    ! The fluid volume fraction of ONE cell, from the signed distance at its
    ! centre and the face-normal it implies. |grad phi| = 1 for a distance
    ! function, so the central differences of phi ARE the interface normal
    ! (the identity the whole scheme rests on, strategy doc Section 8); where
    ! they collapse -- a medial axis, where the two nearest surface points
    ! differ -- there is no usable plane, and the cell falls back to the
    ! pointwise C1 marker, which is what the rest of the scheme does there.
    real(C_DOUBLE) function scalar_volume_fraction(phic, gx, gy, gz, hx, hy, hz) result(f)
        real(C_DOUBLE), intent(in) :: phic, gx, gy, gz   ! phi and its gradient
        real(C_DOUBLE), intent(in) :: hx, hy, hz         ! the cell's edge lengths

        real(C_DOUBLE) :: gn

        gn = sqrt(gx*gx + gy*gy + gz*gz)
        if (gn < CONJ_MIN_GRADPHI) then
            f = merge(1.0d0, 0.0d0, phic > 0.0d0)
            return
        end if
        f = plane_box_fraction(phic, abs(gx)*hx/gn, abs(gy)*hy/gn, abs(gz)*hz/gn)
    end function scalar_volume_fraction

    ! The tangential scalar s_t = e_d.grad T - (n.e_d)(n.grad T) at a cut
    ! face (increment C2; LaTeX note Section 7.3, construction 1). It is the
    ! part of the coordinate-direction derivative that the interface normal
    ! does NOT carry, and it is the one piece of the exact cut-face flux the
    ! C1 baseline drops.
    !
    ! gt* is the face-centred temperature gradient and gp* the face-centred
    ! grad phi, each given as (face-normal direction d, the two tangential
    ! directions) -- so the routine never needs to know WHICH axis d is.
    !
    ! WHY grad phi IS the normal: phi is a signed DISTANCE, so |grad phi| = 1
    ! and no separate normal field is needed (the identity sst%wnorm already
    ! exploits). It is normalised here anyway, because a discrete central
    ! difference only reproduces |grad phi| = 1 to truncation error.
    !
    ! SIGN-INVARIANT under n -> -n (the projection is quadratic in n), so it
    ! does not matter that grad phi points into the fluid while the note's
    ! own n points across the arm.
    !
    ! THE STRADDLING CORRECTION, and why it is here. The note (Section 7.3)
    ! argues that ordinary central differences suffice because the jump in
    ! grad T is purely NORMAL, so the projection removes it. That is right
    ! for the two TANGENTIAL differences, which weight the two sides (1/2,
    ! 1/2), and WRONG for the arm difference, which weights them (w, 1-w).
    ! For a piecewise-linear field with a kink mu in the normal derivative,
    !
    !     grad T^ = G_L + mu [ n/2 + (1/2 - w) n_d e_d ]
    !     =>  s_t^raw = s_t + mu n_d (1/2 - w)(1 - n_d^2),
    !
    ! an O(mu) bias -- NOT small, and growing with the conductivity contrast
    ! that sets mu. C2 measured it: the raw estimate makes the interface flux
    ! WORSE than dropping the term altogether at every tangential ratio below
    ! ~0.1 (validation/conjugate/README.md, gate 1).
    !
    ! The bias is removable in closed form with no wider stencil. Estimating
    ! mu from the face's own normal flux, mu = -q (1/k_R - 1/k_L) with
    ! q = (grad_d T - s_t) k_face/a, makes the n_d cancel and leaves a linear
    ! equation for s_t:
    !
    !     s_t = (s_t^raw - c grad_d T)/(1 - c),
    !     c = kappa_face (1/kappa_R - 1/kappa_L)(1/2 - w)(1 - n_d^2).
    !
    ! It is UNCONDITIONALLY well posed: whenever c > 0 the geometry forces
    ! kappa_face <= 1/w (or 1/(1-w)) and hence c <= 1/2, so 1 - c >= 1/2 for
    ! every material pair and every cut position. And it is inert exactly
    ! where it should be -- equal materials, w = 1/2, a grid-aligned face,
    ! and a purely tangential field (s_t^raw = grad_d T gives s_t = grad_d T
    ! back, whatever c is).
    !
    ! DEVIATION from the note, which offers same-side least squares as the
    ! fallback: that needs the far cell's own far neighbour, i.e. a two-deep
    ! halo -- the same comm-layer change the TVD increment is blocked on.
    ! This costs four lines and no data.
    real(C_DOUBLE) function conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
            kfrac, phiL, phiR, ks, invd) result(st)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: gtd, gt1, gt2      ! grad T at the face
        real(C_DOUBLE), intent(in) :: gpd, gp1, gp2      ! grad phi at the face
        real(C_DOUBLE), intent(in) :: kfrac              ! k_face/dm, dimensionless
        real(C_DOUBLE), intent(in) :: phiL, phiR, ks, invd

        real(C_DOUBLE) :: gn, dn, nd, w, cc

        gn = sqrt(gpd*gpd + gp1*gp1 + gp2*gp2)
        if (gn < CONJ_MIN_GRADPHI) then
            st = 0.0d0                                   ! medial axis: no usable normal
            return
        end if
        nd = gpd/gn
        dn = (gpd*gtd + gp1*gt1 + gp2*gt2)/gn            ! n.grad T
        st = gtd - nd*dn                                 ! the raw projection
        w = conjugate_face_weight(phiL, phiR, invd)
        cc = kfrac*(1.0d0/merge(ks, 1.0d0, phiR < 0.0d0) &
                  - 1.0d0/merge(ks, 1.0d0, phiL < 0.0d0)) &
            *(0.5d0 - w)*(1.0d0 - nd*nd)
        st = (st - cc*gtd)/(1.0d0 - cc)
    end function conjugate_tangential

    ! THE CONJUGATE FAMILY'S FACE DIFFUSIVITY: the C1 coefficient plus, where
    ! the caller says so (addNut: a FLUID cell's UNCUT face under a
    ! turbulence model), the eddy part on nutFace, the 1/2 (L + R) face eddy
    ! viscosity the caller formed (this family has no wall function, so it
    ! never needs the pair) -- nu_t enters neither the solid nor a cut face
    ! (strategy doc Section 12: at DNS resolution nu_t -> 0 at the wall
    ! anyway, and ibm_aware already zeroes it in solid cells, but the
    ! interface coefficient is molecular by construction and must not be
    ! polluted). One definition for the transport kernel and the statistics
    ! kernels, for the same reason as plain_face_diffusivity. The banded
    ! solid (F2) has its own twin, band_face_total, and the CALLER
    ! dispatches on the scalar once around its faces: a per-face dispatch
    ! inside one helper cost the conjugate kernel 5 % (measured 2026-10-20),
    ! the per-scalar branch nothing.
    real(C_DOUBLE) function conjugate_face_total(dm, phL, phR, ks, rc, invd, &
            addNut, nutFace, pr, prt, prtModel, re) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, phL, phR, ks, rc, invd, nutFace, pr, prt, re
        logical, intent(in) :: addNut
        integer(C_INT), intent(in) :: prtModel

        d = conjugate_face_diffusivity(dm, phL, phR, ks, rc, invd)
        if (addNut) d = d + eddy_diffusivity(nutFace, pr, prt, prtModel, re)
    end function conjugate_face_total

    real(C_DOUBLE) function band_face_total(dm, phL, phR, ks, ko, dsh, rc, invd, &
            addNut, nutFace, pr, prt, prtModel, re) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, phL, phR, ks, ko, dsh, rc, invd, nutFace, pr, prt, re
        logical, intent(in) :: addNut
        integer(C_INT), intent(in) :: prtModel

        d = band_face_diffusivity(dm, phL, phR, ks, ko, dsh, rc, invd)
        if (addNut) d = d + eddy_diffusivity(nutFace, pr, prt, prtModel, re)
    end function band_face_total

    ! Build the signed distance phi = +dwall (fluid) / -dwall (solid) at
    ! every cell centre INCLUDING THE GHOSTS, which is the entire geometric
    ! input of the conjugate scheme.
    !
    ! Nothing here is new computation: the magnitude is the RANS wall
    ! distance, use-associated from its two existing producers
    ! (fill_body_distance_analytic for the analytic indicator, the case
    ! file's dwall_blocks tiles for the file path), and the sign is the
    ! existing cell-centred IBM marker. That is the whole point of the
    ! baseline -- no new dataset, no moby_prepare change, no case-file format
    ! change (strategy doc Section 8 / Section 11).
    !
    ! THE MARKER: set_ibm_coeff writes SOLID/Re at every cell whose CENTRE is
    ! inside the body and a graded coefficient only at fluid-centred cells,
    ! so `|coef(VAR_P)| > SOLID_FACE_THRESHOLD` is exactly "this cell centre
    ! is in the solid" -- on both geometry paths, and ghost-inclusive on
    ! both, which is what lets a cut face on a block boundary see both signs.
    !
    ! HOST CODE writing a device-mapped array: the caller must push phi (and
    ! blk%q, if init_scalar_solid_fields ran) to the device afterwards. It
    ! must also have refreshed the HOST copy of ibm%coef first -- the
    ! analytic coefficients are computed on the device.
    subroutine init_scalar_conjugate(sc, dns, blk, bc, ibm, c)
        type(scalar_type), intent(inout) :: sc
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        type(ibm_type), intent(in) :: ibm
        type(comm_type), intent(in) :: c

        real(C_DOUBLE), allocatable :: dwall(:,:,:,:)
        integer :: i, j, k, b, nx, ny, nz
        logical :: found
        real(C_DOUBLE) :: counts(2)

        if (.not. scalar_conjugate_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        allocate(dwall(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))

        call read_dwall_blocks(dwall, found, dns, blk, c%has_terminal)
        if (.not. found) then
            if (c%has_terminal) print *, "error: ibm_wall = conjugate needs the wall", &
                " distance, but the case file has no dwall_blocks; re-run", &
                " moby_prepare (it writes them by default): ", trim(dns%case_file)
            error stop "conjugate needs dwall_blocks"
        end if

        ! A solid cell is given a STRICTLY negative phi, so that `phi < 0` is
        ! an exact material test even for a cell centre lying exactly on the
        ! surface (where dwall is 0). The kernel's cut test is that material
        ! disagreement, which is why the degenerate case cannot leak into it.
        counts = 0.0d0
        do b = 1, int(blk%nBlocks)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    if (abs(ibm%coef(i,j,k,VAR_P,b)) > SOLID_FACE_THRESHOLD) then
                        sc%phi(i,j,k,b) = -max(dwall(i,j,k,b), CONJ_MIN_DISTANCE)
                        if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny &
                            .and. k >= 1 .and. k <= nz) counts(1) = counts(1) + 1.0d0
                    else
                        sc%phi(i,j,k,b) = dwall(i,j,k,b)
                    end if
                end do
            end do
        end do
        end do
        deallocate(dwall)

        call fill_volume_fraction(sc, blk)

        counts(2) = real(count_cut_faces(sc, blk), C_DOUBLE)
        call comm_allreduce_sum(c, counts)
        if (c%has_terminal) then
            print '(a,i0,a,i0,a)', &
                " conjugate interface: ", nint(counts(1)), " solid cells, ", &
                nint(counts(2)), " cut faces"
            ! Say it out loud: which flux the run is applying at those faces.
            print '(a,l1,a,i0,a)', "    tangential correction: ", &
                sc%nTangential > 0_C_INT, " (", sc%nTangential, " scalars)"
        end if

        call check_conjugate_refinement(sc, blk, c)
        call check_scalar_bands(sc, blk, c)
    end subroutine init_scalar_conjugate

    ! The second material band (F2): report how the solid splits, and check
    ! the ONE precondition the band arithmetic has.
    !
    ! phi is 1-Lipschitz, so an arm of length h changes it by at most h and a
    ! face can only ever join ADJACENT bands -- unless the band is thinner
    ! than the arm, in which case a face joins fluid to the outer material
    ! directly and the shell has holes. The check is stated as the thing
    ! itself (does any face skip a band?), not as a proxy on the spacing, so
    ! it holds on stretched and refined grids without a margin argument.
    subroutine check_scalar_bands(sc, blk, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c

        integer :: i, j, k, b, is, nx, ny, nz, m0
        real(C_DOUBLE) :: dsh, counts(3)

        if (sc%nBanded == 0_C_INT) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do is = 1, int(sc%n)
            dsh = sc%solidDepth(is)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE .or. dsh <= 0.0d0) cycle
            counts = 0.0d0
            do b = 1, int(blk%nBlocks)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        m0 = band_material(sc%phi(i,j,k,b), dsh)
                        if (m0 == 1) counts(1) = counts(1) + 1.0d0
                        if (m0 == 2) counts(2) = counts(2) + 1.0d0
                        if (abs(m0 - band_material(sc%phi(i-1,j,k,b), dsh)) == 2 .or. &
                            abs(m0 - band_material(sc%phi(i+1,j,k,b), dsh)) == 2 .or. &
                            abs(m0 - band_material(sc%phi(i,j-1,k,b), dsh)) == 2 .or. &
                            abs(m0 - band_material(sc%phi(i,j+1,k,b), dsh)) == 2 .or. &
                            abs(m0 - band_material(sc%phi(i,j,k-1,b), dsh)) == 2 .or. &
                            abs(m0 - band_material(sc%phi(i,j,k+1,b), dsh)) == 2) &
                            counts(3) = counts(3) + 1.0d0
                    end do
                end do
            end do
            end do
            call comm_allreduce_sum(c, counts)
            if (c%has_terminal) print '(a,i0,a,f0.4,a,i0,a,i0,a)', &
                "    scalar ", is, ": solid band at depth ", dsh, " -- ", &
                nint(counts(1)), " shell cells, ", nint(counts(2)), " outer cells"
            if (nint(counts(3)) > 0) then
                if (c%has_terminal) print '(a,i0,a,i0,a)', &
                    " error: [scalar.", is, &
                    "] solid_thickness is thinner than the grid there: ", &
                    nint(counts(3)), " cells have a face joining the fluid" // &
                    " straight to the outer band"
                error stop "[scalar.N] solid_thickness below the cell size"
            end if
        end do
    end subroutine check_scalar_bands

    ! The fluid volume fraction at every cell (increment C3), from phi and
    ! its central differences. Pure geometry, so it is built ONCE here, right
    ! after phi, and is thereafter read by the transport kernel (the
    ! capacity), by the explicit time-step limiter (the same capacity) and by
    ! the statistics -- one number, three consumers.
    !
    ! GHOSTS: the outermost layer has no central difference of its own, so it
    ! keeps the pointwise marker value (0 or 1). Nothing reads it -- the
    ! capacity is a CELL property and only interior cells are advanced -- but
    ! it makes the field it rides out to the snapshot well defined everywhere.
    ! The inner ghost layer gets the real fraction, which is what a checker
    ! rebuilding a rank-independent statement needs.
    subroutine fill_volume_fraction(sc, blk)
        type(scalar_type), intent(inout) :: sc
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: gx, gy, gz, hx, hy, hz

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do b = 1, int(blk%nBlocks)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    if (i == 0 .or. i == nx+1 .or. j == 0 .or. j == ny+1 &
                        .or. k == 0 .or. k == nz+1) then
                        sc%vfrac(i,j,k,b) = merge(0.0d0, 1.0d0, sc%phi(i,j,k,b) < 0.0d0)
                        cycle
                    end if
                    gx = (sc%phi(i+1,j,k,b) - sc%phi(i-1,j,k,b))*sc%cdx(i,b)
                    gy = (sc%phi(i,j+1,k,b) - sc%phi(i,j-1,k,b))*sc%cdy(j,b)
                    gz = (sc%phi(i,j,k+1,b) - sc%phi(i,j,k-1,b))*sc%cdz(k,b)
                    hx = blk%x(i+1,VAR_U,b) - blk%x(i,VAR_U,b)
                    hy = blk%y(j+1,VAR_V,b) - blk%y(j,VAR_V,b)
                    hz = blk%z(k+1,VAR_W,b) - blk%z(k,VAR_W,b)
                    sc%vfrac(i,j,k,b) = scalar_volume_fraction(sc%phi(i,j,k,b), &
                        gx, gy, gz, hx, hy, hz)
                end do
            end do
        end do
        end do
    end subroutine fill_volume_fraction

    ! Interior cut faces of one rank (the two cell centres in different
    ! materials), for the init report. Counted on the low side of each cell
    ! so no face is counted twice.
    integer function count_cut_faces(sc, blk) result(n)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        logical :: s0

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        n = 0
        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    s0 = sc%phi(i,j,k,b) < 0.0d0
                    if ((sc%phi(i-1,j,k,b) < 0.0d0) .neqv. s0) n = n + 1
                    if ((sc%phi(i,j-1,k,b) < 0.0d0) .neqv. s0) n = n + 1
                    if ((sc%phi(i,j,k-1,b) < 0.0d0) .neqv. s0) n = n + 1
                end do
            end do
        end do
        end do
    end function count_cut_faces

    ! LANDMINE (strategy doc Section 12): refine_body's one-block
    ! 26-neighbour buffer is what keeps cut cells away from 2:1 block
    ! interfaces. The cut-face coefficient is a SAME-LEVEL two-point arm; on
    ! a coarse/fine block face the halo value is a restriction or a prolong
    ! of the other level, so neither phi nor the flux means what the scheme
    ! assumes -- and phi's own halo there is evaluated at THIS block's level's
    ! ghost coordinate, which is not a neighbouring cell centre at all. Away
    ! from the interface that is harmless (only the SIGN is read, and a point
    ! well inside one material has the right one), which is precisely what
    ! this check enforces. The buffer also keeps curvature*h small, which is
    ! what makes the level-set weight w trustworthy. CHECK it -- do not
    ! assume it.
    subroutine check_conjugate_refinement(sc, blk, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c

        integer :: i, j, k, b, d, side, n(3), lo(3), hi(3), p(3), q(3)
        integer(C_INT) :: kind
        real(C_DOUBLE) :: bad(1)

        n = int(blk%nb(1:3))
        bad = 0.0d0

        ! Every 2:1 face of every block: the cell pair straddling it is
        ! (0, 1) on the low side and (n, n+1) on the high side along d, at
        ! every interior index of the two tangential directions.
        do b = 1, int(blk%nBlocks)
            do d = 1, 3
                do side = 0, 1
                    kind = merge(blk%physHigh(d,b), blk%physLow(d,b), side == 1)
                    if (kind /= FACE_COARSE .and. kind /= FACE_FINE) cycle
                    lo = 1; hi = n
                    lo(d) = merge(n(d), 0, side == 1); hi(d) = lo(d)
                    do k = lo(3), hi(3)
                        do j = lo(2), hi(2)
                            do i = lo(1), hi(1)
                                p = [i, j, k]; q = p; q(d) = p(d) + 1
                                if (face_crosses_material(sc, sc%phi(p(1),p(2),p(3),b), &
                                    sc%phi(q(1),q(2),q(3),b))) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call comm_allreduce_sum(c, bad)
        if (bad(1) <= 0.0d0) return
        if (c%has_terminal) print '(a,i0,a)', &
            " error: ibm_wall = conjugate: ", nint(bad(1)), " material faces sit on a" // &
            " 2:1 block face. The conjugate face coefficient is a same-level arm and" // &
            " cannot read across a refinement interface -- for the body surface, and" // &
            " equally for a solid_thickness band boundary. Use [blocks] refine_body =" // &
            " true (its one-block 26-neighbour buffer keeps the surface strictly inside" // &
            " the finest level), or refine the region explicitly."
        error stop "conjugate interface crosses a 2:1 block face"
    end subroutine check_conjugate_refinement

    ! Do these two cell centres hold different MATERIALS -- for any conjugate
    ! scalar? The body surface (phi = 0) is one material boundary and every
    ! scalar shares it; a solid_thickness band adds a second, at phi = -d, and
    ! d is per scalar. Both are same-level arms, so both are equally unable to
    ! sit on a 2:1 block face, and the 2:1 precondition must test both.
    !
    ! Without a band this reduces to the sign test it replaces: band_material
    ! returns 0 in the fluid and 1 everywhere in the solid when dsh = 0.
    logical function face_crosses_material(sc, phiL, phiR) result(bad)
        type(scalar_type), intent(in) :: sc
        real(C_DOUBLE), intent(in) :: phiL, phiR

        integer :: is

        bad = (phiL < 0.0d0) .neqv. (phiR < 0.0d0)
        if (bad .or. sc%nBanded == 0_C_INT) return
        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) cycle
            if (sc%solidDepth(is) <= 0.0d0) cycle
            if (band_material(phiL, sc%solidDepth(is)) &
                /= band_material(phiR, sc%solidDepth(is))) then
                bad = .true.
                return
            end if
        end do
    end function face_crosses_material

    ! The solid's own initial temperature ([scalar.N] solid_init, default =
    ! `initial`). COLD START ONLY: on a restart the solid field is part of
    ! the saved state and comes from the file like every other cell.
    ! Ghost-inclusive, so the halo layer is consistent before the first
    ! exchange. HOST code writing a device-mapped blk%q -- the caller pushes.
    subroutine init_scalar_solid_fields(sc, blk)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk

        integer :: is, i, j, k, b, nx, ny, nz

        if (.not. scalar_conjugate_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) cycle
            do b = 1, int(blk%nBlocks)
            do k = 0, nz+1
                do j = 0, ny+1
                    do i = 0, nx+1
                        if (sc%phi(i,j,k,b) < 0.0d0) &
                            blk%q(i,j,k,VAR_S0+is,b) = sc%solidInit(is)
                    end do
                end do
            end do
            end do
        end do
    end subroutine init_scalar_solid_fields

    ! The conjugate body's contribution to the explicit diffusive (Peclet)
    ! rate, max'ed into dns%peclet_rate by the caller once the interface
    ! exists.
    !
    ! TWO THINGS GATE 1 TAUGHT, both of which cost a NaN before they were
    ! understood, and neither of which a per-material bound can see:
    !
    ! 1. THE LIMIT IS NOT max OVER MATERIALS OF alpha = kappa/C. A cut face
    !    carries k_face up to max(kappa_L, kappa_R) -- that bound IS the
    !    strategy doc Section 7 argument that the interface is not stiff --
    !    but it feeds the cell on the OTHER side, whose capacity belongs to
    !    the other material. So a fluid cell against a kappa_s = 1000 solid
    !    sees 1000x the fluid rate even when alpha_s = alpha_f exactly.
    !    Building the rate from the ACTUAL face coefficients is both correct
    !    and far less conservative than max(kappa)/min(C): at w = 1/2 and
    !    kappa_s = 1000 the true penalty is 2x, not 1000x.
    !
    ! 2. THE CONVENTION IS THE GERSHGORIN ONE: HALF THE DIAGONAL OVER THE
    !    CAPACITY. The operator has zero row sum, so its spectral radius is
    !    bounded by 2*A_ii/C_i, and this RK3's real-axis limit of 2.51 gives
    !    dt <= 1.256 C_i/A_ii. The rate reported is A_ii/(2 C_i): for a
    !    regular cell that is exactly sum_d alpha/h_d^2, the number
    !    `pecletmax` bounds everywhere (step.f90 RK3_DIFFUSION_LIMIT), so the
    !    limit on pecletmax is the same 0.628 here as in the rest of the
    !    solver and a cut cell needs no rule of its own.
    !    HISTORY, because the gates carry numbers taken under it: until
    !    2026-10-02 `pecletmax` bounded the single-direction alpha/h^2, a
    !    regular cell reported diag/6 and a cut cell a special share --
    !    diag/3 in C1 (measured: (kappa_s, C_s) = (0.01, 0.01) at w = 0.95
    !    blew up at pecletmax 0.3; (1000, 1000) at w = 0.80 at 0.4), diag/2
    !    from C3 on (an OBLIQUE interface at kappa_s = 10^3 went to NaN at 0.4
    !    with diag/3). diag/2 is the present rule, so the cut-cell step of
    !    every C3 result is unchanged; what changed is the BULK, which went
    !    from diag/6 to diag/2. The old bulk convention survived only where
    !    one direction carried the diagonal or the Courant limit bound first:
    !    on an isotropic grid the extreme mode is excited from round-off
    !    (validation/diffusion_limit/).
    !
    ! The capacity here is the C3 FLUID-FRACTION-WEIGHTED one, the same the
    ! transport kernel divides by -- and it moves the limit at exactly the
    ! cells with the least margin, in both directions: a fluid cut cell
    ! against a heavy solid gains capacity (a larger step), a solid cut cell
    ! loses it (a smaller one). Reading it from sc%vfrac rather than
    ! re-deriving it is what keeps the limiter and the kernel talking about
    ! the same cell.
    real(C_DOUBLE) function scalar_conjugate_peclet_rate(sc, blk, dns, c) result(rate)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        integer :: i, j, k, b, is, nx, ny, nz
        real(C_DOUBLE) :: dm, ks, rc, cc, phc, diag, r(1)
        real(C_DOUBLE) :: dsh, ko, csb
        logical :: tangOn, banded

        rate = 0.0d0
        if (.not. scalar_conjugate_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        r = 0.0d0

        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) cycle
            dm = 1.0d0/(dns%re*sc%pr(is))
            ks = sc%solidK(is)
            rc = sc%contactR(is)
            tangOn = sc%tangCorr(is) /= 0_C_INT
            ! The banded solid (F2) reaches the limiter through the same two
            ! numbers the kernel uses: the face coefficients, and the cell's
            ! own capacity. An INSULATED outer band therefore contributes a
            ! rate of exactly zero -- which is right, it is inert.
            dsh = sc%solidDepth(is)
            ko = sc%outerK(is)
            banded = dsh > 0.0d0
            do b = 1, int(blk%nBlocks)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        phc = sc%phi(i,j,k,b)
                        csb = sc%solidC(is)
                        if (banded) then
                            if (band_material(phc, dsh) == 2) csb = sc%outerC(is)
                        end if
                        cc = sc%vfrac(i,j,k,b) &
                            + (1.0d0 - sc%vfrac(i,j,k,b))*csb
                        if (banded) then
                            diag = (band_face_diffusivity(dm, sc%phi(i-1,j,k,b), phc, &
                                        ks, ko, dsh, rc, sc%invDx(i,b))*sc%invDx(i,b) &
                                  + band_face_diffusivity(dm, phc, sc%phi(i+1,j,k,b), &
                                        ks, ko, dsh, rc, sc%invDx(i+1,b))*sc%invDx(i+1,b)) &
                                    *blk%d1x(i,VAR_P,b) &
                                 + (band_face_diffusivity(dm, sc%phi(i,j-1,k,b), phc, &
                                        ks, ko, dsh, rc, sc%invDy(j,b))*sc%invDy(j,b) &
                                  + band_face_diffusivity(dm, phc, sc%phi(i,j+1,k,b), &
                                        ks, ko, dsh, rc, sc%invDy(j+1,b))*sc%invDy(j+1,b)) &
                                    *blk%d1y(j,VAR_P,b) &
                                 + (band_face_diffusivity(dm, sc%phi(i,j,k-1,b), phc, &
                                        ks, ko, dsh, rc, sc%invDz(k,b))*sc%invDz(k,b) &
                                  + band_face_diffusivity(dm, phc, sc%phi(i,j,k+1,b), &
                                        ks, ko, dsh, rc, sc%invDz(k+1,b))*sc%invDz(k+1,b)) &
                                    *blk%d1z(k,VAR_P,b)
                        else
                            diag = (conjugate_face_diffusivity(dm, sc%phi(i-1,j,k,b), phc, &
                                        ks, rc, sc%invDx(i,b))*sc%invDx(i,b) &
                                  + conjugate_face_diffusivity(dm, phc, sc%phi(i+1,j,k,b), &
                                        ks, rc, sc%invDx(i+1,b))*sc%invDx(i+1,b)) &
                                    *blk%d1x(i,VAR_P,b) &
                                 + (conjugate_face_diffusivity(dm, sc%phi(i,j-1,k,b), phc, &
                                        ks, rc, sc%invDy(j,b))*sc%invDy(j,b) &
                                  + conjugate_face_diffusivity(dm, phc, sc%phi(i,j+1,k,b), &
                                        ks, rc, sc%invDy(j+1,b))*sc%invDy(j+1,b)) &
                                    *blk%d1y(j,VAR_P,b) &
                                 + (conjugate_face_diffusivity(dm, sc%phi(i,j,k-1,b), phc, &
                                        ks, rc, sc%invDz(k,b))*sc%invDz(k,b) &
                                  + conjugate_face_diffusivity(dm, phc, sc%phi(i,j,k+1,b), &
                                        ks, rc, sc%invDz(k+1,b))*sc%invDz(k+1,b)) &
                                    *blk%d1z(k,VAR_P,b)
                        end if
                        ! The C2 tangential correction is an EXPLICIT spatial
                        ! operator at the same cut faces, so it enters the
                        ! same rate. It is not diffusive -- it couples this
                        ! cell to its TANGENTIAL neighbours with a sign the
                        ! diagonal does not balance -- so what is added is
                        ! its Gershgorin row sum, which bounds both the real
                        ! and the imaginary parts of its spectrum. Its size
                        ! grows with the conductivity contrast (LaTeX note
                        ! Section 8.3): at a face whose midpoint is in a
                        ! kappa_s = 1000 solid, |k_loc - k_face| is ~1000 dm
                        ! against a principal term of ~2 dm. That penalty is
                        ! real, and gate 4 measures it.
                        if (tangOn) diag = diag &
                            + (face_corr_rate(sc, is, i-1, j, k, b, 1, dm, ks, rc) &
                             + face_corr_rate(sc, is, i,   j, k, b, 1, dm, ks, rc)) &
                                *blk%d1x(i,VAR_P,b) &
                            + (face_corr_rate(sc, is, i, j-1, k, b, 2, dm, ks, rc) &
                             + face_corr_rate(sc, is, i, j,   k, b, 2, dm, ks, rc)) &
                                *blk%d1y(j,VAR_P,b) &
                            + (face_corr_rate(sc, is, i, j, k-1, b, 3, dm, ks, rc) &
                             + face_corr_rate(sc, is, i, j, k,   b, 3, dm, ks, rc)) &
                                *blk%d1z(k,VAR_P,b)
                        r(1) = max(r(1), diag/(2.0d0*cc))
                    end do
                end do
            end do
            end do
        end do

        call comm_allreduce_max(c, r)
        rate = r(1)
    end function scalar_conjugate_peclet_rate

    ! Gershgorin row sum of the C2 tangential correction at ONE face, the
    ! face between cell (i,j,k) and its neighbour one step up in direction d
    ! (so the caller passes the LOW cell of the face). Zero at an uncut face,
    ! where the correction is identically zero.
    !
    ! The correction contributes s_t (k_loc - k_face) to this face's flux and
    !   s_t = (1 - n_d^2) grad_d T - n_d n_1 grad_1 T - n_d n_2 grad_2 T,
    ! whose stencil coefficients are +-invd (twice, the arm) and +-cd/2 (four
    ! times, the two cells' tangential central differences) -- hence the 2s
    ! below. Multiplied by the caller's 1/dx it is the row sum this face adds.
    ! HOST ONLY: it runs once, at init, beside the face diffusivities.
    real(C_DOUBLE) function face_corr_rate(sc, is, i, j, k, b, d, dm, ks, rc) result(rt)
        type(scalar_type), intent(in) :: sc
        integer, intent(in) :: is, i, j, k, b, d
        real(C_DOUBLE), intent(in) :: dm, ks, rc

        integer :: t1, t2, o1, o2, o3, p1, q1, r1, p2, q2, r2
        real(C_DOUBLE) :: phl, phr, invd, cd1, cd2, gpd, gp1, gp2, gn, nd, n1, n2

        rt = 0.0d0
        o1 = merge(1, 0, d == 1)
        o2 = merge(1, 0, d == 2)
        o3 = merge(1, 0, d == 3)
        phl = sc%phi(i,j,k,b)
        phr = sc%phi(i+o1,j+o2,k+o3,b)
        if ((phl < 0.0d0) .eqv. (phr < 0.0d0)) return

        t1 = merge(2, 1, d == 1)     ! increasing order, as in the kernel
        t2 = merge(2, 3, d == 3)
        p1 = merge(1, 0, t1 == 1); q1 = merge(1, 0, t1 == 2); r1 = merge(1, 0, t1 == 3)
        p2 = merge(1, 0, t2 == 1); q2 = merge(1, 0, t2 == 2); r2 = merge(1, 0, t2 == 3)
        select case (d)
        case (1);      invd = sc%invDx(i+1,b)
        case (2);      invd = sc%invDy(j+1,b)
        case default;  invd = sc%invDz(k+1,b)
        end select
        select case (t1)
        case (1);      cd1 = sc%cdx(i,b)
        case (2);      cd1 = sc%cdy(j,b)
        case default;  cd1 = sc%cdz(k,b)
        end select
        select case (t2)
        case (1);      cd2 = sc%cdx(i,b)
        case (2);      cd2 = sc%cdy(j,b)
        case default;  cd2 = sc%cdz(k,b)
        end select

        gpd = (phr - phl)*invd
        gp1 = 0.5d0*((sc%phi(i+p1,j+q1,k+r1,b) - sc%phi(i-p1,j-q1,k-r1,b)) &
                   + (sc%phi(i+o1+p1,j+o2+q1,k+o3+r1,b) &
                    - sc%phi(i+o1-p1,j+o2-q1,k+o3-r1,b)))*cd1
        gp2 = 0.5d0*((sc%phi(i+p2,j+q2,k+r2,b) - sc%phi(i-p2,j-q2,k-r2,b)) &
                   + (sc%phi(i+o1+p2,j+o2+q2,k+o3+r2,b) &
                    - sc%phi(i+o1-p2,j+o2-q2,k+o3-r2,b)))*cd2
        gn = sqrt(gpd*gpd + gp1*gp1 + gp2*gp2)
        if (gn < CONJ_MIN_GRADPHI) return          ! the correction is off there too
        nd = gpd/gn; n1 = gp1/gn; n2 = gp2/gn

        rt = abs(conjugate_face_local_k(dm, phl, phr, ks) &
               - conjugate_face_diffusivity(dm, phl, phr, ks, rc, invd)) &
            *2.0d0*((1.0d0 - nd*nd)*invd + abs(nd*n1)*cd1 + abs(nd*n2)*cd2)
    end function face_corr_rate

    ! e_face = h_d s_t/(T_R - T_L) at cut faces, reduced to max and rms
    ! (strategy doc Section 3, LaTeX note Eq. (indicator)). This is the
    ! MEASUREMENT the increment exists for: the C1 baseline drops the
    ! tangential term of the exact cut-face flux, and
    !
    !     (q_n^num - q_n)/q_n = h_d s_t/(T_R - T_L - h_d s_t)
    !
    ! so e_face IS the relative error the baseline makes in the
    ! interface-normal flux at that face, to leading order. It costs the same
    ! s_t the correction needs, so a run can report the error it is making
    ! whether or not it is correcting it -- which is the whole point.
    !
    ! e_face = s_t/gtd with gtd = (T_R - T_L)/h_d: the metric cancels, so
    ! nothing here depends on the mesh spacing.
    !
    ! Each face is visited ONCE, as the LOW face of its high cell, which
    ! covers every face except a block's outermost high faces (the low face
    ! of the neighbouring block's first cell -- present unless that face is a
    ! physical boundary, where cut faces would be pathological anyway).
    !
    ! TWO PASSES, because the denominator can legitimately vanish: a cut face
    ! carrying no normal flux has e_face = infinity while meaning nothing.
    ! Pass 1 finds the largest |grad_d T| over cut faces; pass 2 accumulates
    ! only where the denominator is above 1e-6 of it, and reports how many
    ! faces that excluded so the reader can see whether the statistic is
    ! representative.
    subroutine scalar_conjugate_indicator(sc, blk, c, step)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c
        integer(C_INT), intent(in) :: step

        real(C_DOUBLE), parameter :: DENOM_FRACTION = 1.0d-6
        real(C_DOUBLE) :: gmax(1), emax(1), acc(3), rms
        integer :: is

        if (.not. scalar_conjugate_enabled(sc)) return

        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) cycle
            call indicator_pass(sc, blk, is, -1.0d0, gmax(1), emax(1), acc)
            call comm_allreduce_max(c, gmax)
            call indicator_pass(sc, blk, is, DENOM_FRACTION*gmax(1), gmax(1), &
                emax(1), acc)
            call comm_allreduce_max(c, emax)
            call comm_allreduce_sum(c, acc)
            rms = 0.0d0
            if (acc(2) > 0.0d0) rms = sqrt(acc(1)/acc(2))
            if (c%has_terminal) print '(a,i8,a,a,a,i9,a,es11.4,a,es11.4,a,i0,a)', &
                " conjugate indicator step", step, "  '", trim(sc%name(is)), &
                "': cut faces", nint(acc(2) + acc(3)), &
                "  max|e_face| =", emax(1), "  rms =", rms, &
                "  (", nint(acc(3)), " below the denominator floor)"
        end do
    end subroutine scalar_conjugate_indicator

    ! One reduction sweep of the indicator over the LOW face of every
    ! interior cell in all three directions. thresh < 0 asks for the
    ! denominator scale alone (pass 1); otherwise the face is accumulated
    ! when |grad_d T| > thresh and counted as skipped when it is not.
    ! acc = (sum e^2, faces used, faces skipped).
    subroutine indicator_pass(sc, blk, is, thresh, gmax, emax, acc)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: is
        real(C_DOUBLE), intent(in) :: thresh
        real(C_DOUBLE), intent(out) :: gmax, emax, acc(3)

        integer :: i, j, k, b, d, t1, t2, nx, ny, nz, nBlocks, var
        integer :: il, jl, kl, o1, o2, o3, p1, q1, r1, p2, q2, r2
        real(C_DOUBLE) :: gtd, gt1, gt2, gpd, gp1, gp2, invd, cd1, cd2, st, e
        real(C_DOUBLE) :: gmx, emx, esum, nuse, nskip, ks, rc

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        var = VAR_S0 + is
        ! The face weight and the straddling correction need the solid's
        ! properties; the diffusivity itself is wanted only as the ratio
        ! k_face/dm, so it is formed with dm = 1.
        ks = sc%solidK(is)
        rc = sc%contactR(is)
        gmx = 0.0d0; emx = 0.0d0; esum = 0.0d0; nuse = 0.0d0; nskip = 0.0d0

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& reduction(max: gmx, emx) reduction(+: esum, nuse, nskip) &
        !$omp& map(to: nx, ny, nz, var, thresh, ks, rc, blk%q, sc%phi, &
        !$omp& sc%invDx, sc%invDy, sc%invDz, sc%cdx, sc%cdy, sc%cdz) &
        !$omp& private(i,j,k,b,d,t1,t2,il,jl,kl,o1,o2,o3,p1,q1,r1,p2,q2,r2, &
        !$omp& gtd,gt1,gt2,gpd,gp1,gp2,invd,cd1,cd2,st,e)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do d = 1, 3
                        ! The low face in direction d, and the unit offsets of
                        ! d and of the two directions tangential to it.
                        o1 = merge(1, 0, d == 1)
                        o2 = merge(1, 0, d == 2)
                        o3 = merge(1, 0, d == 3)
                        il = i - o1; jl = j - o2; kl = k - o3
                        if ((sc%phi(il,jl,kl,b) < 0.0d0) .eqv. &
                            (sc%phi(i,j,k,b) < 0.0d0)) cycle
                        ! The two tangential directions in INCREASING order,
                        ! the transport kernel's convention (2026-10-20: a y
                        ! face used (z, x) here, which changes the round-off
                        ! of conjugate_tangential's dot products).
                        t1 = merge(2, 1, d == 1)
                        t2 = merge(2, 3, d == 3)
                        p1 = merge(1, 0, t1 == 1)
                        q1 = merge(1, 0, t1 == 2)
                        r1 = merge(1, 0, t1 == 3)
                        p2 = merge(1, 0, t2 == 1)
                        q2 = merge(1, 0, t2 == 2)
                        r2 = merge(1, 0, t2 == 3)
                        select case (d)
                        case (1);      invd = sc%invDx(i,b)
                        case (2);      invd = sc%invDy(j,b)
                        case default;  invd = sc%invDz(k,b)
                        end select
                        select case (t1)
                        case (1);      cd1 = sc%cdx(i,b)
                        case (2);      cd1 = sc%cdy(j,b)
                        case default;  cd1 = sc%cdz(k,b)
                        end select
                        select case (t2)
                        case (1);      cd2 = sc%cdx(i,b)
                        case (2);      cd2 = sc%cdy(j,b)
                        case default;  cd2 = sc%cdz(k,b)
                        end select
                        ! The same face-centred gradients the transport
                        ! kernel forms for the correction, in the same
                        ! low -> high order.
                        gtd = (blk%q(i,j,k,var,b) - blk%q(il,jl,kl,var,b))*invd
                        gt1 = 0.5d0*((blk%q(il+p1,jl+q1,kl+r1,var,b) &
                                    - blk%q(il-p1,jl-q1,kl-r1,var,b)) &
                                   + (blk%q(i+p1,j+q1,k+r1,var,b) &
                                    - blk%q(i-p1,j-q1,k-r1,var,b)))*cd1
                        gt2 = 0.5d0*((blk%q(il+p2,jl+q2,kl+r2,var,b) &
                                    - blk%q(il-p2,jl-q2,kl-r2,var,b)) &
                                   + (blk%q(i+p2,j+q2,k+r2,var,b) &
                                    - blk%q(i-p2,j-q2,k-r2,var,b)))*cd2
                        gpd = (sc%phi(i,j,k,b) - sc%phi(il,jl,kl,b))*invd
                        gp1 = 0.5d0*((sc%phi(il+p1,jl+q1,kl+r1,b) &
                                    - sc%phi(il-p1,jl-q1,kl-r1,b)) &
                                   + (sc%phi(i+p1,j+q1,k+r1,b) &
                                    - sc%phi(i-p1,j-q1,k-r1,b)))*cd1
                        gp2 = 0.5d0*((sc%phi(il+p2,jl+q2,kl+r2,b) &
                                    - sc%phi(il-p2,jl-q2,kl-r2,b)) &
                                   + (sc%phi(i+p2,j+q2,k+r2,b) &
                                    - sc%phi(i-p2,j-q2,k-r2,b)))*cd2
                        st = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                            conjugate_face_diffusivity(1.0d0, sc%phi(il,jl,kl,b), &
                                sc%phi(i,j,k,b), ks, rc, invd), &
                            sc%phi(il,jl,kl,b), sc%phi(i,j,k,b), ks, invd)
                        gmx = max(gmx, abs(gtd))
                        if (thresh < 0.0d0) cycle
                        if (abs(gtd) > thresh) then
                            e = st/gtd
                            emx = max(emx, abs(e))
                            esum = esum + e*e
                            nuse = nuse + 1.0d0
                        else
                            nskip = nskip + 1.0d0
                        end if
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        gmax = gmx
        emax = emx
        acc(1) = esum
        acc(2) = nuse
        acc(3) = nskip
    end subroutine indicator_pass

    ! Push the host-filled signed distance to the device (a no-op on the CPU
    ! build and without a conjugate scalar).
    subroutine scalar_conjugate_to_device(sc)
        type(scalar_type), intent(in) :: sc

        if (.not. scalar_conjugate_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target update to(sc%phi, sc%vfrac)
#endif
    end subroutine scalar_conjugate_to_device

end module scalar_conjugate
