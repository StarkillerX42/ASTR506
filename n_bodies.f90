module n_bodies
    use integrators
    implicit none
    contains


    function adjust_timestep(dt, err, q, t_tot)
        implicit none
        real, intent(in) :: dt, err, t_tot
        integer, intent(in) :: q
        real :: adjust_timestep
        real :: fac
!        real :: facmax=6., facmin=0.01
!        real :: fac
        fac = 0.38**(1. / (1. + q))
        adjust_timestep = dt * fac * err**(-1. / (1. + q))
!        if (adjust_timestep > facmax * dt) then
!            adjust_timestep = facmax * dt
!        else if (adjust_timestep < facmin * dt) then
!            adjust_timestep = facmin * dt
        if ((adjust_timestep / t_tot < 1e-12) .or. isnan(adjust_timestep)) then
            adjust_timestep = t_tot * 1e-12
        end if
    end function adjust_timestep

    function center_of_mass(n, masses, positions)
        implicit none
        integer, intent(in) :: n
        real, intent(in), dimension(n) :: masses
        real, intent(in), dimension(n, 3) :: positions
        real, dimension(3) :: center_of_mass
        integer :: i
        center_of_mass = 0.
        do i=1, n
            center_of_mass(:) = center_of_mass(:) + masses(i) * positions(i, :)
        end do
    end function center_of_mass

    function average_vel(n, masses, velocities)
        implicit none
        integer, intent(in) :: n
        real, intent(in), dimension(n) :: masses
        real, intent(in), dimension(n, 3) :: velocities
        real, dimension(3) :: average_vel
        integer :: i
        average_vel = 0.
        do i=1, n
            average_vel(:) = average_vel(:) + masses(i) * velocities(i, :)
        end do
    end function average_vel

    subroutine n_body_model(integrator, n, masses, i_pos, i_vel, t_tot, dt, &
                            g, times, positions, velocities)
        implicit none
        character(len=64), intent(in) :: integrator
        integer :: n
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(n, 3), intent(in) :: i_pos, i_vel
        real(8), dimension(n, 3) :: poses, vels, null_poses, null_vels
        real(8), dimension(n, 3) :: tmp_poses, tmp_vels, accels
        real(8), dimension(3) :: com, avv
        real(8), intent(in) :: t_tot, dt
        real(8), intent(in) :: G
        real(8) :: t_now
        integer :: t, i, j, n_steps
        real(8), intent(inout) :: times(:)
        ! shape n_steps, n, 3
        real(8), dimension(:, :, :), intent(inout) :: positions, velocities
        n_steps = int(t_tot/dt)
!        allocate(times(n_steps))
!        allocate(positions(n_steps, n, 3))
!        allocate(velocities(n_steps, n, 3))


        ! Centers the system inputs
!        com(:) = center_of_mass(n, masses, i_pos)
!        avv(:) = average_vel(n, masses, i_vel)
!        do i=1, n
!            poses(i, :) = i_pos(i, :) - com
!            vels(i, :) = i_vel(i, :) - avv
!        end do

        timeloop: do t=1, n_steps
            times(t) = t_now
            positions(t, :, :) = poses(:, :)
            velocities(t, :, :) = vels(:, :)
            ! Run the integrator of choice for a timestep, using temporary
            ! variables to avoid simultaneous IO for variables
            if (integrator == "euler") then
                call euler_forward(n, masses, poses, vels, dt, tmp_poses, &
                        tmp_vels)
                poses = tmp_poses
                vels = tmp_vels
            else if (integrator == "kdk") then
                call kdk(n, masses, poses, vels, dt, tmp_poses, tmp_vels)
                poses = tmp_poses
                vels = tmp_vels
            else if (integrator == "rk3") then
                call rk3_classic(n, masses, poses, vels, dt, tmp_poses, tmp_vels)
                poses = tmp_poses
                vels = tmp_vels
            else if (integrator == "rk4") then
                call rkf45(n, masses, poses, vels, dt, tmp_poses, tmp_vels, &
                           null_poses, null_vels)
                poses = tmp_poses
                vels = tmp_vels
            else if (integrator == "rk5") then
                call rkf45(n, masses, poses, vels, dt, null_poses, null_vels, &
                        tmp_poses, tmp_vels)
                poses = tmp_poses
                vels = tmp_vels
                t_now = t_now + dt
            end if

            ! Check validity of the integrator's calculation
            do i=1, n
                do j=1, 3
                    if (isnan(poses(i, j))) then
                        print *, "NAN encountered for position at", i, j
                        exit timeloop
                    end if
                    if (isnan(vels(i, j))) then
                        print *, "NAN encountered for velocity at", i, j
                        exit timeloop
                    end if
                end do
            end do

            ! Recenters the system
!            com = center_of_mass(n, masses, poses)
!            avv = average_vel(n, masses, vels)
!            do i=1, n
!                poses(i, :) = poses(i, :) - com(:)
!                vels(i, :) = vels(i, :) - avv(:)
!            end do


        end do timeloop

    end subroutine n_body_model


end module n_bodies