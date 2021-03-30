module integrators
    implicit none
    contains


    subroutine vec_mag(n_dim, rel_vec, rel_mag)
        integer, intent(in) :: n_dim
        real, dimension(n_dim), intent(in) :: rel_vec
        real, intent(out) :: rel_mag
        rel_mag = sqrt(sum(rel_vec**2))
    end subroutine vec_mag

    subroutine euler_forward(n, masses, positions, velocities, dt, poses, vels, accels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3), intent(out) :: accels
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
!        print *, "Accels"
!        print *, velocities(2, :) + accels(2, :) * dt
        poses(:, :) = positions + velocities * dt
        vels(:, :) = velocities(:, :) + accels(:, :) * dt
    end subroutine euler_forward


    subroutine kdk(n, masses, positions, velocities, dt, poses, vels, accels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3), intent(out) :: accels
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
                        rel = poses(j, :) - poses(i, :)
                        call vec_mag(3, rel, rel_mag)
                        accels(i, :) = accels(i, :) + masses(j) / rel_mag**3 * rel(:)
                    end if
                end if
            end do
            vels(i, :) = vels(i, :) + accels(i, :) * dt / 2.  ! Kick 2
        end do

    end subroutine kdk


    subroutine rk_pde(n, masses, positions, velocities, dt, vels, accels)
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
        vels = velocities

    end subroutine rk_pde


    subroutine rk3_classic(n, masses, positions, velocities, dt, poses, vels, accels)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3) :: null_poses
        real, dimension(n, 3), intent(out) :: poses, vels, accels
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

        call euler_forward(n, masses, positions, velocities, dt, null_poses, &
                k_poses(1, :, :), k_vels(1, :, :))
        do i=1, 3
            temp_p = positions
            temp_v = velocities
            do j=1, n  ! Particles
                do k=1, 3  ! Dimension
                    do l=1, i  ! weight index
                        temp_p(j, k) = temp_p(j, k) + k_poses(l, j, k) * coefficients(i, l)
                        temp_v(j, k) = temp_v(j, k) + k_vels(l, j, k) * coefficients(i, l)
                    end do
                end do
            end do
!            print *, shape(k_poses), "|", shape(reshape(coefficients, (/3, n, 3/)))
            call euler_forward(n, masses, temp_p, temp_v, dt, null_poses, &
                        k_poses(i+1, :, :), k_vels(i+1, :, :))
!            print *, k_poses(i + 1, 2, :)
        end do

        poses = positions
        vels = velocities
        do i=1, 4
            poses(:, :) = poses(:, :) + dt * weights(i) * k_poses(i, :, :)
            vels(:, :) = vels(:, :) + dt * weights(i) * k_vels(i, :, :)
        end do

    end subroutine rk3_classic


    subroutine rkf45(n, masses, positions, velocities, dt, poses_low, &
            vels_low, poses_high, vels_high)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: masses
        real, dimension(3) :: rel
        real, dimension(n, 3), intent(in) :: positions, velocities
        real, dimension(n, 3) :: accels
        real, dimension(n, 3), intent(out) :: poses_low, vels_low, poses_high, &
                vels_high
        integer :: i, j, k, l
        real, intent(in) :: dt
        real :: rel_mag
        real, dimension(5, 5) :: coefficients
        real, dimension(6) :: weights_4, weights_5
        real, dimension(6, n, 3) :: k_poses, k_vels
        real, dimension(n, 3) :: temp_p, temp_v

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


        call rk_pde(n, masses, positions, velocities, dt, k_poses(1, :, :), k_vels(1, :, :))
        do i=1, 5
            temp_p = positions
            temp_v = velocities
            do j=1, n  ! Particles
                do k=1, 3  ! Dimension
                    do l=1, i  ! weight index
                        temp_p(j, k) = temp_p(j, k) + k_poses(l, j, k) * coefficients(i, l)
                        temp_v(j, k) = temp_v(j, k) + k_vels(l, j, k) * coefficients(i, l)
                    end do
                end do
            end do
!            print *, shape(k_poses), "|", shape(reshape(coefficients, (/3, n, 3/)))
            call rk_pde(n, masses, temp_p, temp_v, dt, &
                    k_poses(i+1, :, :), k_vels(i+1, :, :))
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
