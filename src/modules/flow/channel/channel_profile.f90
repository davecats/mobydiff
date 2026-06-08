module channel_profile
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W
    implicit none

    private

    real(C_DOUBLE), parameter :: PI = 3.1415926535897932384626433832795d0
    integer, parameter :: LARGE_DISTURBANCE_NMODES = 5

    public :: initialise_channel_fields

contains

    subroutine initialise_channel_fields(f, dns, g, n_walls, large_amp, noise_amp)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: large_amp, noise_amp

        integer :: i, j, k, var, nx, ny, nz
        real(C_DOUBLE) :: stream_x, span_x, wall_y, envelope, laminar_u
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

                        call laminar_profile(dns, n_walls, wall_y, laminar_u, envelope)

                        large = large_disturbance(large_amp, envelope, stream_x, span_x, &
                            dns%leng(1), dns%leng(3))
                        noise = noise_amp*envelope * &
                            deterministic_noise(i + int(dns%localSize(1,0)) - 1, &
                                                j + int(dns%localSize(2,0)) - 1, &
                                                k + int(dns%localSize(3,0)) - 1, var)

                        f%q(i,j,k,var) = stream_profile(var, laminar_u) + &
                            disturbance_component(var, large) + noise
                    end do
                end do
            end do
        end do
    end subroutine initialise_channel_fields

    subroutine laminar_profile(dns, n_walls, wall_y, laminar_u, envelope)
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: wall_y
        real(C_DOUBLE), intent(out) :: laminar_u, envelope

        if (n_walls == 2) then
            laminar_u = 0.5d0*dns%re*dns%forcing(VAR_U)*wall_y*(dns%leng(2) - wall_y)
            envelope = sin(PI*wall_y/dns%leng(2))
        else
            laminar_u = dns%re*dns%forcing(VAR_U)*(dns%leng(2)*wall_y - 0.5d0*wall_y*wall_y)
            envelope = sin(0.5d0*PI*wall_y/dns%leng(2))
        end if
    end subroutine laminar_profile

    real(C_DOUBLE) function stream_profile(var, laminar_u) result(value)
        integer, intent(in) :: var
        real(C_DOUBLE), intent(in) :: laminar_u

        value = merge(laminar_u, 0.0d0, var == VAR_U)
    end function stream_profile

    real(C_DOUBLE) function large_disturbance(amplitude, envelope, stream_x, span_x, &
            stream_length, span_length) result(value)
        real(C_DOUBLE), intent(in) :: amplitude, envelope, stream_x, span_x
        real(C_DOUBLE), intent(in) :: stream_length, span_length

        integer :: mode
        real(C_DOUBLE) :: stream_phase, span_phase, weight, norm

        stream_phase = 2.0d0*PI*stream_x/max(stream_length, 1.0d-12)
        span_phase = 2.0d0*PI*span_x/max(span_length, 1.0d-12)

        value = 0.0d0
        norm = 0.0d0
        do mode = 1, LARGE_DISTURBANCE_NMODES
            weight = 1.0d0/real(mode, C_DOUBLE)
            value = value + weight * &
                sin(real(mode, C_DOUBLE)*stream_phase) * &
                cos(real(mode, C_DOUBLE)*span_phase)
            norm = norm + weight
        end do

        if (norm > 0.0d0) value = amplitude*envelope*value/norm
    end function large_disturbance

    real(C_DOUBLE) function disturbance_component(var, large) result(value)
        integer, intent(in) :: var
        real(C_DOUBLE), intent(in) :: large

        if (var == VAR_U) then
            value = large
        else if (var == VAR_V) then
            value = 0.5d0*large
        else if (var == VAR_W) then
            value = -0.5d0*large
        else
            value = 0.0d0
        end if
    end function disturbance_component

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
