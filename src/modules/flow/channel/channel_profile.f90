module channel_profile
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W
    implicit none

    private

    real(C_DOUBLE), parameter :: PI = 3.1415926535897932384626433832795d0

    public :: initialise_channel_fields
    public :: channel_span_dir, channel_wall_coordinate, channel_plane_coord

contains

    subroutine initialise_channel_fields(f, dns, g, stream_dir, wall_dir, n_walls, large_amp, noise_amp)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        integer(C_INT), intent(in) :: stream_dir, wall_dir
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: large_amp, noise_amp

        integer :: i, j, k, var, nx, ny, nz
        integer(C_INT) :: span_dir
        real(C_DOUBLE) :: stream_x, span_x, wall_y, envelope, laminar_u
        real(C_DOUBLE) :: large, noise

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        span_dir = channel_span_dir(stream_dir, wall_dir)

        f%q = 0.0d0
        f%qs = 0.0d0
        f%oldrhs = 0.0d0

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do var = VAR_U, VAR_W
                        stream_x = staggered_coord(g, stream_dir, i, j, k, var)
                        span_x = staggered_coord(g, span_dir, i, j, k, var)
                        wall_y = channel_wall_coordinate(g, wall_dir, i, j, k, var)
                        wall_y = max(0.0d0, min(wall_y, dns%leng(wall_dir)))

                        call laminar_profile(dns, stream_dir, wall_dir, n_walls, &
                            wall_y, laminar_u, envelope)

                        large = large_disturbance(large_amp, envelope, stream_x, span_x, &
                            dns%leng(stream_dir), dns%leng(span_dir))
                        noise = noise_amp*envelope * &
                            deterministic_noise(i + int(dns%localSize(1,0)) - 1, &
                                                j + int(dns%localSize(2,0)) - 1, &
                                                k + int(dns%localSize(3,0)) - 1, var)

                        f%q(i,j,k,var) = stream_profile(var, stream_dir, laminar_u) + &
                            disturbance_component(var, stream_dir, wall_dir, span_dir, large) + noise
                    end do
                end do
            end do
        end do
    end subroutine initialise_channel_fields

    integer(C_INT) function channel_span_dir(stream_dir, wall_dir) result(span_dir)
        integer(C_INT), intent(in) :: stream_dir, wall_dir
        integer :: dir

        span_dir = 3_C_INT
        do dir = 1, 3
            if (dir /= stream_dir .and. dir /= wall_dir) then
                span_dir = int(dir, C_INT)
                return
            end if
        end do
    end function channel_span_dir

    real(C_DOUBLE) function channel_wall_coordinate(g, wall_dir, i, j, k, var) result(y)
        type(grid_type), intent(in) :: g
        integer(C_INT), intent(in) :: wall_dir
        integer, intent(in) :: i, j, k, var

        y = staggered_coord(g, wall_dir, i, j, k, var)
    end function channel_wall_coordinate

    real(C_DOUBLE) function channel_plane_coord(g, wall_dir, n) result(coord)
        type(grid_type), intent(in) :: g
        integer(C_INT), intent(in) :: wall_dir
        integer, intent(in) :: n

        select case (wall_dir)
        case (1)
            coord = 0.5d0*(g%xNode(n-1) + g%xNode(n))
        case (2)
            coord = 0.5d0*(g%yNode(n-1) + g%yNode(n))
        case default
            coord = 0.5d0*(g%zNode(n-1) + g%zNode(n))
        end select
    end function channel_plane_coord

    real(C_DOUBLE) function staggered_coord(g, dir, i, j, k, var) result(coord)
        type(grid_type), intent(in) :: g
        integer(C_INT), intent(in) :: dir
        integer, intent(in) :: i, j, k, var

        select case (dir)
        case (1)
            coord = g%x(i,var)
        case (2)
            coord = g%y(j,var)
        case default
            coord = g%z(k,var)
        end select
    end function staggered_coord

    subroutine laminar_profile(dns, stream_dir, wall_dir, n_walls, wall_y, laminar_u, envelope)
        type(dns_type), intent(in) :: dns
        integer(C_INT), intent(in) :: stream_dir, wall_dir
        integer, intent(in) :: n_walls
        real(C_DOUBLE), intent(in) :: wall_y
        real(C_DOUBLE), intent(out) :: laminar_u, envelope

        if (n_walls == 2) then
            laminar_u = 0.5d0*dns%re*dns%forcing(stream_dir)*wall_y*(dns%leng(wall_dir) - wall_y)
            envelope = sin(PI*wall_y/dns%leng(wall_dir))
        else
            laminar_u = dns%re*dns%forcing(stream_dir)* &
                (dns%leng(wall_dir)*wall_y - 0.5d0*wall_y*wall_y)
            envelope = sin(0.5d0*PI*wall_y/dns%leng(wall_dir))
        end if
    end subroutine laminar_profile

    real(C_DOUBLE) function stream_profile(var, stream_dir, laminar_u) result(value)
        integer, intent(in) :: var
        integer(C_INT), intent(in) :: stream_dir
        real(C_DOUBLE), intent(in) :: laminar_u

        value = merge(laminar_u, 0.0d0, var == stream_dir)
    end function stream_profile

    real(C_DOUBLE) function large_disturbance(amplitude, envelope, stream_x, span_x, &
            stream_length, span_length) result(value)
        real(C_DOUBLE), intent(in) :: amplitude, envelope, stream_x, span_x
        real(C_DOUBLE), intent(in) :: stream_length, span_length

        value = amplitude*envelope * &
            sin(2.0d0*PI*stream_x/max(stream_length, 1.0d-12)) * &
            cos(2.0d0*PI*span_x/max(span_length, 1.0d-12))
    end function large_disturbance

    real(C_DOUBLE) function disturbance_component(var, stream_dir, wall_dir, span_dir, large) result(value)
        integer, intent(in) :: var
        integer(C_INT), intent(in) :: stream_dir, wall_dir, span_dir
        real(C_DOUBLE), intent(in) :: large

        if (var == stream_dir) then
            value = large
        else if (var == wall_dir) then
            value = 0.5d0*large
        else if (var == span_dir) then
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
