module channel_profile
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W
    implicit none

    private

    real(C_DOUBLE), parameter :: PI = 3.1415926535897932384626433832795d0
    integer, parameter :: LARGE_DISTURBANCE_NMODES = 5
    real(C_DOUBLE), parameter :: LUCHINI_KAPPA = 0.392d0
    real(C_DOUBLE), parameter :: LUCHINI_B = 4.48d0
    real(C_DOUBLE), parameter :: LUCHINI_A0 = -7.374d0
    real(C_DOUBLE), parameter :: LUCHINI_A1 = -0.4930d0
    real(C_DOUBLE), parameter :: LUCHINI_A2 = 0.02450d0
    real(C_DOUBLE), parameter :: LUCHINI_B1 = 0.05736d0
    real(C_DOUBLE), parameter :: LUCHINI_B2 = 0.01101d0
    real(C_DOUBLE), parameter :: LUCHINI_C = 0.03385d0
    real(C_DOUBLE), parameter :: LUCHINI_D = 3.109d0

    public :: initialise_channel_fields

contains

    subroutine initialise_channel_fields(f, dns, g, n_walls, mean_sine_amp, large_amp, noise_amp)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: mean_sine_amp, large_amp, noise_amp

        integer :: i, j, k, var, nx, ny, nz
        real(C_DOUBLE) :: stream_x, span_x, wall_y, envelope, mean_u
        real(C_DOUBLE) :: large, noise

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

        f%q = 0.0d0
        f%qs = 0.0d0
        f%oldrhs = 0.0d0

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do var = VAR_U, VAR_W
                        stream_x = g%x(i,var)
                        span_x = g%z(k,var)
                        wall_y = max(0.0d0, min(g%y(j,var), dns%leng(2)))

                        call mean_profile(dns, n_walls, wall_y, mean_sine_amp, mean_u)
                        envelope = disturbance_envelope(n_walls, var, wall_y, dns%leng(2))

                        large = large_disturbance(large_amp, var, stream_x, wall_y, span_x, &
                            dns%leng(1), dns%leng(2), dns%leng(3), n_walls)
                        noise = noise_amp*envelope * &
                            deterministic_noise(i + int(dns%localSize(1,0)) - 1, &
                                                j + int(dns%localSize(2,0)) - 1, &
                                                k + int(dns%localSize(3,0)) - 1, var)

                        f%q(i,j,k,var) = stream_profile(var, mean_u) + large + noise
                    end do
                end do
            end do
        end do
    end subroutine initialise_channel_fields

    subroutine mean_profile(dns, n_walls, wall_y, mean_sine_amp, mean_u)
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: wall_y, mean_sine_amp
        real(C_DOUBLE), intent(out) :: mean_u

        real(C_DOUBLE) :: wall_length

        wall_length = max(dns%leng(2), 1.0d-12)
        mean_u = turbulent_channel_profile(dns, n_walls, wall_y)
        if (n_walls == 2) then
            mean_u = mean_u + mean_sine_amp*sin(2.0d0*PI*wall_y/wall_length)
        else
            mean_u = mean_u + mean_sine_amp*sin(0.5d0*PI*wall_y/wall_length)
        end if
    end subroutine mean_profile

    real(C_DOUBLE) function turbulent_channel_profile(dns, n_walls, wall_y) result(value)
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: wall_y

        real(C_DOUBLE) :: outer_length, wall_distance, zplus, zouter

        if (n_walls == 2) then
            outer_length = max(0.5d0*dns%leng(2), 1.0d-12)
            wall_distance = max(0.0d0, min(wall_y, dns%leng(2) - wall_y))
        else
            outer_length = max(dns%leng(2), 1.0d-12)
            wall_distance = max(0.0d0, min(wall_y, dns%leng(2)))
        end if

        zouter = min(1.0d0, wall_distance/outer_length)
        zplus = dns%re*zouter

        value = luchini_wall_law(zplus) + luchini_plane_duct_wake(zouter)
    end function turbulent_channel_profile

    real(C_DOUBLE) function luchini_wall_law(zplus) result(value)
        real(C_DOUBLE), intent(in) :: zplus

        real(C_DOUBLE) :: z, pade, denom

        z = max(0.0d0, zplus)
        denom = 1.0d0 + LUCHINI_B1*z + LUCHINI_B2*z*z
        pade = ((LUCHINI_A0 + LUCHINI_A1*z + LUCHINI_A2*z*z)/denom)*exp(-LUCHINI_C*z)
        value = pade + log(z + LUCHINI_D)/LUCHINI_KAPPA + LUCHINI_B
        if (z <= 0.0d0) value = 0.0d0
    end function luchini_wall_law

    real(C_DOUBLE) function luchini_plane_duct_wake(zouter) result(value)
        real(C_DOUBLE), intent(in) :: zouter

        real(C_DOUBLE) :: z

        z = max(0.0d0, min(1.0d0, zouter))
        value = z - 0.57d0*z**7
    end function luchini_plane_duct_wake

    real(C_DOUBLE) function stream_profile(var, laminar_u) result(value)
        integer, intent(in) :: var
        real(C_DOUBLE), intent(in) :: laminar_u

        value = merge(laminar_u, 0.0d0, var == VAR_U)
    end function stream_profile

    real(C_DOUBLE) function disturbance_envelope(n_walls, var, wall_y, wall_length) result(value)
        integer, intent(in) :: n_walls, var
        real(C_DOUBLE), intent(in) :: wall_y, wall_length

        real(C_DOUBLE) :: length

        length = max(wall_length, 1.0d-12)
        if (n_walls == 1 .and. var /= VAR_V) then
            value = sin(0.5d0*PI*wall_y/length)
        else
            value = sin(PI*wall_y/length)
        end if
    end function disturbance_envelope

    real(C_DOUBLE) function large_disturbance(amplitude, var, stream_x, wall_y, span_x, &
            stream_length, wall_length, span_length, n_walls) result(value)
        real(C_DOUBLE), intent(in) :: amplitude, stream_x, wall_y, span_x
        real(C_DOUBLE), intent(in) :: stream_length, wall_length, span_length
        integer, intent(in) :: var, n_walls

        integer :: mx, my, mz
        real(C_DOUBLE) :: stream_phase, wall_phase, span_phase, phase_shift, weight, norm
        real(C_DOUBLE) :: wall_mode

        phase_shift = 0.37d0*real(var - VAR_U, C_DOUBLE)
        value = 0.0d0
        norm = 0.0d0
        do mz = 1, LARGE_DISTURBANCE_NMODES
            span_phase = 2.0d0*PI*real(mz, C_DOUBLE)*span_x/max(span_length, 1.0d-12)
            do my = 1, LARGE_DISTURBANCE_NMODES
                wall_mode = real(my, C_DOUBLE)
                if (n_walls == 1 .and. var /= VAR_V) wall_mode = wall_mode - 0.5d0
                wall_phase = PI*wall_mode*wall_y/max(wall_length, 1.0d-12)
                do mx = 1, LARGE_DISTURBANCE_NMODES
                    stream_phase = 2.0d0*PI*real(mx, C_DOUBLE)*stream_x/max(stream_length, 1.0d-12)
                    weight = 1.0d0/real(mx + my + mz, C_DOUBLE)
                    value = value + weight * sin(stream_phase + phase_shift) * &
                        sin(wall_phase) * cos(span_phase - phase_shift)
                    norm = norm + weight
                end do
            end do
        end do

        if (norm > 0.0d0) value = amplitude*value/norm
    end function large_disturbance

    real(C_DOUBLE) function deterministic_noise(i, j, k, var) result(noise)
        integer, intent(in) :: i, j, k, var

        integer(int64) :: n

        n = int(i, int64)*73856093_int64 + int(j, int64)*19349663_int64 + &
            int(k, int64)*83492791_int64 + int(var, int64)*2654435761_int64
        n = ieor(n, ishft(n, -13))
        n = mod(abs(n*1274126177_int64), 2147483647_int64)
        noise = 2.0d0*(real(n, C_DOUBLE)/2147483647.0d0) - 1.0d0
    end function deterministic_noise

end module channel_profile
