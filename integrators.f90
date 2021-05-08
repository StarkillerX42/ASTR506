module integrators
    implicit none
    contains


    function vec_mag(n_dim, rel_vec)
        integer, intent(in) :: n_dim
        real(8), dimension(n_dim), intent(in) :: rel_vec
        real(8) :: vec_mag
        vec_mag = sqrt(sum(rel_vec**2))
    end function vec_mag

    subroutine euler_forward(n, masses, positions, velocities, dt, g, poses, vels)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        real(8), intent(in) :: g
        real(8), dimension(n, 3) :: accels
        real(8), dimension(n, 3), intent(out):: poses, vels
        real(8), intent(in) :: dt
        integer(8) :: i, j
        real(8) :: rel_mag
        accels = 0.
        poses = 0.
        vels = 0.
        !    f2py intent(in) n, masses, positions, velocities, dt
        !    f2py intent(out) poses, vels
        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0.) then
                        rel = positions(j, :) - positions(i, :)
                        rel_mag = vec_mag(3, rel)
                        accels(i, :) = accels(i, :) + g * masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
            poses(i, :) = positions(i, :) + velocities(i, :) * dt
            vels(i, :) = velocities(i, :) + accels(i, :) * dt
        end do
!        print *, "Accels"
!        print *, velocities(2, :) + accels(2, :) * dt
!        poses(:, :) = positions + velocities * dt
!        vels(:, :) = velocities(:, :) + accels(:, :) * dt
    end subroutine euler_forward


    subroutine kdk(n, masses, positions, velocities, dt, g, poses, vels)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        real(8), dimension(n, 3) :: accels
        real(8), intent(in) :: g
        real(8), dimension(n, 3), intent(out) :: poses, vels
        integer(8) :: i, j
        real(8), intent(in) :: dt
        real(8) :: rel_mag
        rel = 0.
        rel_mag = 0.
        accels = 0.
        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0) then
                        rel = positions(j, :) - positions(i, :)
                        rel_mag = vec_mag(3, rel)
                        accels(i, :) = accels(i, :) + g * masses(j) / rel_mag**3 * rel(:)
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
                        rel = poses(j, :) - poses(i, :)
                        rel_mag = vec_mag(3, rel)
                        accels(i, :) = accels(i, :) + g * masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
            vels(i, :) = vels(i, :) + accels(i, :) * dt / 2.  ! Kick 2
        end do

    end subroutine kdk


    subroutine kdk_gal(n, masses, positions, velocities, dt, g, poses, vels)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        real(8), intent(in) :: g
        real(8), dimension(n, 3) :: accels
        real(8), dimension(n, 3), intent(out) :: poses, vels
        integer(8) :: i
        real(8), intent(in) :: dt
        real(8) :: rel_mag
        rel = 0.
        accels = 0.
        rel_mag = 1.
        do i=1, n
            rel = positions(i, :)
            rel_mag = vec_mag(3, rel)
            accels(i, :) = (-log(1 + rel_mag) / rel_mag ** 2 + 1 / rel_mag / (1 + rel_mag)) * rel(:) / rel_mag
!            print*, accels(i, :)
            vels(i, :) = velocities(i, :) + accels(i, :) * dt / 2.  ! Kick
!            print *, vels(i, :)
            poses(i, :) = positions(i, :) + vels(i, :) * dt  ! Drift
        end do
        accels = 0.  ! Reset accels for second kick

        do i=1, n
            rel = poses(i, :)
            rel_mag = vec_mag(3, rel)
            accels(i, :) = (-log(1 + rel_mag) / rel_mag ** 2 + 1 / rel_mag / (1 + rel_mag)) * rel(:) / rel_mag
            vels(i, :) = vels(i, :) + accels(i, :) * dt / 2.  ! Kick 2
        end do

    end subroutine kdk_gal


    subroutine rk_pde(n, masses, positions, velocities, g, vels, accels)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), intent(in) :: g
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        real(8), dimension(n, 3), intent(out) :: vels, accels
        integer(8) :: i, j
        real(8) :: rel_mag
    !    f2py intent(in) masses, positions, velocities, dt
    !    f2py intent(out) poses, vels
        vels = 0.
        accels = 0.

        do i=1, n
            do j=1, n
                if (i /= j) then
                    if (masses(j) /= 0.) then
                        rel = positions(j, :) - positions(i, :)
!                        print *, "rel and accels", rel, accels(i, :)
                        rel_mag = vec_mag(3, rel)
                        accels(i, :) = accels(i, :) + g * masses(j) / rel_mag**3 * rel(:)
!                        print *, "masse and accels", masses(i), accels(i, :)
                    end if
                end if
            end do
        end do
        vels = velocities

    end subroutine rk_pde


    subroutine rk3_classic(n, masses, positions, velocities, dt, g, poses, vels)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        real(8), dimension(n, 3), intent(out) :: poses, vels
        real(8), intent(in) :: g
        integer(8) :: i, j, k, l
        real(8), intent(in) :: dt
        real(8), dimension(3, 3) :: coefficients
        real(8), dimension(4) :: weights
        real(8), dimension(4, n, 3) :: k_poses, k_vels
        real(8), dimension(n, 3) :: temp_p, temp_v
        k_poses = 0.
        k_vels = 0.
        temp_p = 0.
        temp_v = 0.
        coefficients(1, :) = (/0.5, 0., 0./)
        coefficients(2, :) = (/0., 0.5, 0./)
        coefficients(3, :) = (/0., 0., 1./)
        weights(:) = (/1., 2., 2., 1./) / 6.
        temp_p = 0.
        temp_v = 0.

        call rk_pde(n, masses, positions, velocities, g, &
                    k_poses(1, :, :), k_vels(1, :, :))
        do i=1, 3
            temp_p = positions
            temp_v = velocities
            do l=1, i
                temp_p = temp_p + k_poses(i, :, :) * coefficients(i, l) * dt
                temp_v = temp_v + k_vels(i, :, :) * coefficients(i, l) * dt
            end do
            call rk_pde(n, masses, temp_p, temp_v, g, &
                        k_poses(i + 1, :, :), k_vels(i + 1, :, :))
        end do

        poses = positions
        vels = velocities
        do i=1, 4
!            print "(e22.16)", k_poses(i, :, :)
!            print *, k_poses(i, 2, 2)
            poses(:, :) = poses(:, :) + dt * weights(i) * k_poses(i, :, :)
            vels(:, :) = vels(:, :) + dt * weights(i) * k_vels(i, :, :)
        end do

    end subroutine rk3_classic


    subroutine kdkrk3(n, masses, positions, velocities, dt, g, poses_low, &
            vels_low, poses_high, vels_high)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        !        real(8), dimension(n, 3) :: accels
        real(8), dimension(n, 3), intent(out) :: poses_low, vels_low, poses_high, &
                vels_high
        real(8), intent(in) :: g
        integer(8) :: i, j, k, l
        real(8), intent(in) :: dt
        real(8) :: rel_mag
        real(8), dimension(3, 3) :: coefficients
        real(8), dimension(4) :: weights
        real(8), dimension(4, n, 3) :: k_poses, k_vels
        real(8), dimension(n, 3) :: temp_p, temp_v
        k_poses = 0.
        k_vels = 0.
        temp_p = 0.
        temp_v = 0.

        coefficients(1, :) = (/0.5, 0., 0./)
        coefficients(2, :) = (/0., 0.5, 0./)
        coefficients(3, :) = (/0., 0., 1./)
        weights(:) = (/1., 2., 2., 1./) / 6.

!        call rk_pde(n, masses, positions, velocities, &
!                k_poses(1, :, :), k_vels(1, :, :))
!        do i=1, 5
!            temp_p = positions
!            temp_v = velocities
!            do l=1, i
!                temp_p = temp_p + k_poses(i, :, :) * coefficients(i, l) * dt
!                temp_v = temp_v + k_vels(i, :, :) * coefficients(i, l) * dt
!            end do
!            call rk_pde(n, masses, temp_p, temp_v, &
!                    k_poses(i + 1, :, :), k_vels(i + 1, :, :))
!        end do
!
!        poses_low = positions
!        vels_low = velocities
!        poses_high = positions
!        vels_high = velocities
!
!        do i=1, 6
!            poses_low(:, :) = poses_low(:, :) + dt * weights_4(i) * k_poses(i, :, :)
!            vels_low(:, :) = vels_low(:, :) + dt * weights_4(i) * k_vels(i, :, :)
!            poses_high(:, :) = poses_high(:, :) + dt * weights_5(i) * k_poses(i, :, :)
!            vels_high(:, :) = vels_high(:, :) + dt * weights_5(i) * k_vels(i, :, :)
!        end do

    end subroutine kdkrk3


    subroutine rkf45(n, masses, positions, velocities, dt, g, poses_low, &
            vels_low, poses_high, vels_high)
        implicit none
        integer(8), intent(in) :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(3) :: rel
        real(8), dimension(n, 3), intent(in) :: positions, velocities
        real(8), intent(in) :: g
!        real(8), dimension(n, 3) :: accels
        real(8), dimension(n, 3), intent(out) :: poses_low, vels_low, poses_high, &
                vels_high
        integer(8) :: i, j, k, l
        real(8), intent(in) :: dt
        real(8) :: rel_mag
        real(8), dimension(5, 5) :: coefficients
        real(8), dimension(6) :: weights_4, weights_5
        real(8), dimension(6, n, 3) :: k_poses, k_vels
        real(8), dimension(n, 3) :: temp_p, temp_v
        k_poses = 0.
        k_vels = 0.
        temp_p = 0.
        temp_v = 0.

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

        temp_p = 0.
        temp_v = 0.

        call rk_pde(n, masses, positions, velocities, g, &
                    k_poses(1, :, :), k_vels(1, :, :))
        do i=1, 5
            temp_p = positions
            temp_v = velocities
            do l=1, i
                temp_p = temp_p + k_poses(i, :, :) * coefficients(i, l) * dt
                temp_v = temp_v + k_vels(i, :, :) * coefficients(i, l) * dt
            end do
            call rk_pde(n, masses, temp_p, temp_v, g, &
                        k_poses(i + 1, :, :), k_vels(i + 1, :, :))
        end do

        poses_low = positions
        vels_low = velocities
        poses_high = positions
        vels_high = velocities

        do i=1, 6
            poses_low(:, :) = poses_low(:, :) + dt * weights_4(i) * k_poses(i, :, :)
            vels_low(:, :) = vels_low(:, :) + dt * weights_4(i) * k_vels(i, :, :)
            poses_high(:, :) = poses_high(:, :) + dt * weights_5(i) * k_poses(i, :, :)
            vels_high(:, :) = vels_high(:, :) + dt * weights_5(i) * k_vels(i, :, :)
        end do

    end subroutine rkf45

end module integrators
