module integrators

implicit none
contains


    subroutine vec_mag(n_dim, rel_vec, rel_mag)
        integer, intent(in) :: n_dim
        real, dimension(n_dim), intent(in) :: rel_vec
        real, intent(out) :: rel_mag
        rel_mag = sqrt(sum(rel_vec**2))
    end subroutine vec_mag


    subroutine euler_forward(n, masses, positions, velocities, dt, poses, vels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3) :: accels
        real, dimension(n, 3), intent(out) :: poses, vels
        integer :: i, j
        real :: dt, rel_mag
    !    f2py intent(in) n, masses, positions, velocities, dt
    !    f2py intent(out) poses, vels
        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0.) then
                        rel = positions(j, :) - positions(i, :)
                        call vec_mag(3, rel, rel_mag)
                        accels(i, :) = accels(i, :) + masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
        end do
        poses = positions + velocities * dt
        vels = velocities + accels * dt
    end subroutine euler_forward


    subroutine adjust_timestep(dt, err, q, t_tot, dt_new)
        implicit none
        real, intent(in) :: dt, err, t_tot
        integer, intent(in) :: q
        real, intent(out) :: dt_new
        real :: facmax=100., facmin=0.01
        real :: fac
        fac = 0.38**(1. / (1. + q))
        dt_new = dt * fac * err**(-1. / (1. + q))
        if (dt_new > facmax * dt) then
            dt_new = facmax * dt
        else if (dt_new < facmin * dt) then
            dt_new = facmin * dt
        else if (dt_new / t_tot < 1e-12) then
            dt_new = t_tot * 1e-12
        end if
    end subroutine adjust_timestep


    subroutine kdk(n, masses, positions, velocities, dt, poses, vels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3) :: accels
        real, dimension(n, 3), intent(out) :: poses, vels
        integer :: i, j
        real, intent(in) :: dt
        real :: rel_mag
        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0) then
                        rel = positions(j, :) - positions(i, :)
                        call vec_mag(3, rel, rel_mag)
                        accels(i, :) = accels(i, :) + masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
            vels(i, :) = velocities(i, :) + accels(i, :) * dt / 2.  ! Kick
            poses(i, :) = positions(i, :) + vels(i, :) * dt  ! Drift
        end do
        accels = 0.  ! Reset accels for second kick

        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0) then
                        rel = positions(j, :) - positions(i, :)
                        call vec_mag(3, rel, rel_mag)
                        accels(i, :) = accels(i, :) + masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
            vels(i, :) = vels(i, :) + accels(i, :) * dt / 2.  ! Kick 2
        end do

    end subroutine kdk


    subroutine rk_euler(n, masses, positions, velocities, dt, vels, accels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3), intent(out) :: vels, accels
        integer :: i, j
        real :: dt, rel_mag
    !    f2py intent(in) masses, positions, velocities, dt
    !    f2py intent(out) poses, vels

        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0.) then
                        rel = positions(j, :) - positions(i, :)
                        call vec_mag(3, rel, rel_mag)
                        accels(i, :) = accels(i, :) + masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
        end do
        vels = velocities + accels * dt

    end subroutine rk_euler


    subroutine rk3_classic(n, masses, positions, velocities, dt, poses, vels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3) :: accels
        real, dimension(n, 3), intent(out) :: poses, vels
        integer :: i, j, k, l
        real, intent(in) :: dt
        real :: rel_mag
        real, dimension(3, 3) :: coefficients
        real, dimension(4) :: weights
        real, dimension(4, n, 3) :: k_poses, k_vels
        real, dimension(n, 3) :: temp_p, temp_v
        coefficients(1, :) = (/0.5, 0., 0./)
        coefficients(2, :) = (/0., 0.5, 0./)
        coefficients(3, :) = (/0., 0., 1./)
        weights(:) = (/1., 2., 2., 1./) / 6.
        temp_p = 0.
        temp_v = 0.

        call rk_euler(n, masses, positions, velocities, dt, k_poses(1, :, :), k_vels(1, :, :))
        do i=1, 3
            temp_p = positions
            temp_v = velocities
            do j=1, n  ! Particles
                do k=1, 3  ! Dimension
                    do l=1, i  ! weight index
                        print *, k_poses(l, j, k)
                        print *, i, j, k, l
                        temp_p(j, k) = temp_p(j, k) + k_poses(l, j, k) * coefficients(i, l)
                    end do
                end do
            end do
            print *, shape(k_poses), "|", shape(reshape(coefficients, (/3, n, 3/)))
!            call rk_euler(n, masses, &
!                    positions + k_poses(:-1, :, :) * coefficients(i, :), &
!                    velocities + k_vels(:-1, :, :) * coefficients(i, :), dt, &
!                    k_poses(i+1, :, :), k_vels(i+1, :, :))
        end do

    end subroutine rk3_classic


    subroutine rkf45(n, masses, positions, velocities, dt, poses, vels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3) :: accels
        real, dimension(n, 3), intent(out) :: poses, vels
        integer :: i, j
        real, intent(in) :: dt
        real :: rel_mag
        real, dimension(5, 5) :: coefficients
        real, dimension(6) :: weights_4, weights_5

        coefficients(1, :) = (/1. / 4., 0., 0., 0., 0./)
        coefficients(2, :) = (/3. / 32., 9. / 32., 0., 0., 0./)
        coefficients(3, :) = (/1932. / 2197., -7200. / 2197., 7296. / 2197., &
                               0., 0./)
        coefficients(4, :) = (/439. / 216., -8., 3680. / 513., -845. / 4104., &
                               0./)
        coefficients(5, :) = (/-8. / 27., 2., -3544. / 2565., 1859. / 4104., &
                            -11. / 40./)
        weights_4(:) = (/25. / 216., 0., 1408. / 2565., 2197. / 4104., &
                         -1. / 5., 0./)
        weights_5(:) = (/16. / 135., 0., 6656. / 12825., 28561. / 56430.,&
                         -9 / 50., 2. / 55./)


    end subroutine rkf45

end module integrators
