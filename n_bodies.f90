module n_bodies
    use integrators
    implicit none
    contains


    function adjust_timestep(dt, err, q, t_tot)
        implicit none
        real(8), intent(in) :: dt, err, t_tot
        integer(8), intent(in) :: q
        real(8) :: adjust_timestep
        real(8) :: fac
!        real :: facmax=6., facmin=0.01
!        real :: fac
        fac = 0.38**(1. / (1. + q))
        adjust_timestep = dt * fac * err**(-1. / (1. + q))
!        print *, "Adjust timestep", adjust_timestep, dt, fac, err
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
        integer(8), intent(in) :: n
        real(8), intent(in), dimension(n) :: masses
        real(8), intent(in), dimension(n, 3) :: positions
        real(8), dimension(3) :: center_of_mass
        integer(8) :: i
        center_of_mass = 0.
        do i=1, n
            center_of_mass(:) = center_of_mass(:) + masses(i) * positions(i, :)
        end do
    end function center_of_mass

    function average_vel(n, masses, velocities)
        implicit none
        integer(8), intent(in) :: n
        real(8), intent(in), dimension(n) :: masses
        real(8), intent(in), dimension(n, 3) :: velocities
        real(8), dimension(3) :: average_vel
        integer(8) :: i
        average_vel = 0.
        do i=1, n
            average_vel(:) = average_vel(:) + masses(i) * velocities(i, :)
        end do
    end function average_vel


    subroutine n_body_model(integrator, n, n_times, masses, i_pos, i_vel, &
                            t_tot, dt, g, update_com, times, positions, &
                            velocities, &
                            status)
        implicit none
        character(265), intent(in) :: integrator
        integer(8) :: n, n_times, n_times_internal
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(n, 3), intent(in) :: i_pos, i_vel
        real(8), dimension(n, 3) :: poses, vels, null_poses, null_vels
        real(8), dimension(n, 3) :: tmp_poses, tmp_vels, accels
        real(8), dimension(3) :: com, avv
        real(8), intent(in) :: t_tot, dt
        real(8), intent(in) :: G
        logical, intent(in) :: update_com
        real(8) :: t_now
        integer(8) :: t, i, j
        real(8), dimension(n_times + 1), intent(out) :: times
        ! shape n_steps, n, 3
        real(8), dimension(n_times + 1, n, 3), intent(out) :: positions, velocities
        character(256), intent(out) :: status
        status = "Ok"

        if (n_times < int(t_tot/dt)) then
            print *, "SizeError, n_times too small for t_tot/dt"
            status = "SizeError, n_times too small for t_tot/dt"
            stop
        else
            n_times_internal = floor(t_tot / dt)
        end if
!        allocate(times(n_steps))
!        allocate(positions(n_steps, n, 3))
!        allocate(velocities(n_steps, n, 3))
!
!
        ! Centers the system inputs
        if (update_com) then
            com(:) = center_of_mass(n, masses, i_pos)
            avv(:) = average_vel(n, masses, i_vel)
            do i=1, n
                poses(i, :) = i_pos(i, :) - com
                vels(i, :) = i_vel(i, :) - avv
            end do
        end if

        timeloop: do t=1, n_times_internal
!            print *, t
            times(t) = t_now
            do i=1, n
                positions(t, i, :) = poses(i, :)
!                print *, poses(i, :)
                velocities(t, i, :) = vels(i, :)
            end do
!             Run the integrator of choice for a timestep, using temporary
!             variables to avoid simultaneous IO for variables
            if (integrator(1:len("euler")) == "euler") then
                call euler_forward(n, masses, poses, vels, dt, tmp_poses, &
                        tmp_vels)
!                print *, poses(2, :), tmp_poses(2, :)
!                print *, vels(2, :), tmp_vels(2, :)
                do i=1, n
                    poses(i, :) = tmp_poses(i, :)
                    vels(i, :) = tmp_vels(i, :)
                end do
            else if (integrator(1:len("kdk")) == "kdk") then
                call kdk(n, masses, poses, vels, dt, tmp_poses, tmp_vels)
                poses = tmp_poses
                vels = tmp_vels
            else if (integrator(1:len("rk3")) == "rk3") then
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
            end if
            t_now = t_now + dt

!             Check validity of the integrator's calculation
            do i=1, n
                do j=1, 3
                    if (isnan(poses(i, j))) then
                        print *, "NAN encountered for position at", i, j
                        status = "NAN encountered for position"
                        exit timeloop
                    end if
                    if (isnan(vels(i, j))) then
                        print *, "NAN encountered for velocity at", i, j
                        status = "NAN encountered for velocity"
                        exit timeloop
                    end if
                end do
            end do
!
!            print *, vels(2, :)
            ! Recenters the system
            if (update_com) then
                com(:) = center_of_mass(n, masses, i_pos)
                avv(:) = average_vel(n, masses, i_vel)
                do i=1, n
                    poses(i, :) = poses(i, :) - com(:)
                    vels(i, :) = vels(i, :) - avv(:)
                end do
            end if

!            print *, vels(2, :)



        end do timeloop

    end subroutine n_body_model


    subroutine adaptive_n_body_model(integrator, n, n_times, masses, i_pos, &
                                     i_vel, t_tot, tolerance, dt_0, q, G, &
                                     update_com, times, positions, &
                                     velocities, status)
        implicit none
        character(265), intent(in) :: integrator
        integer(8) :: n, n_times, n_times_internal
        real(8), dimension(n), intent(in) :: masses
        real(8), dimension(n, 3), intent(in) :: i_pos, i_vel
        real(8), dimension(n, 3) :: poses, vels, high_poses, high_vels
        real(8), dimension(n, 3) :: low_poses, low_vels, accels
        real(8), dimension(3) :: com, avv
        real(8), intent(in) :: t_tot, tolerance, G, dt_0
        integer(8), intent(in) :: q
        logical, intent(in) :: update_com
        real(8) :: t_now, err_sum, dt
        integer(8) :: t, i, j
        real(8), dimension(n_times + 1), intent(out) :: times
        ! shape n_steps, n, 3
        real(8), dimension(n_times + 1, n, 3), intent(out) :: positions, velocities
        character(256), intent(out) :: status

        status = "Ok"
        dt = dt_0
        if (n_times < int(t_tot/dt)) then
            print *, "SizeError, n_times too small for t_tot/dt"
            status = "SizeError, n_times too small for t_tot/dt"
        else
            n_times_internal = floor(t_tot / dt)
        end if
       ! allocate(times(n_steps))
       ! allocate(positions(n_steps, n, 3))
       ! allocate(velocities(n_steps, n, 3))
!
!
        ! Centers the system inputs
        if (update_com) then
            com(:) = center_of_mass(n, masses, i_pos)
            avv(:) = average_vel(n, masses, i_vel)
            do i=1, n
                poses(i, :) = i_pos(i, :) - com
                vels(i, :) = i_vel(i, :) - avv
            end do
        end if
!        print *, poses(2, :)
!        print *, "Looping for...", n_times_internal
        timeloop: do t=1, n_times_internal
       !     print *, t
            times(t) = t_now
            do i=1, n
                positions(t, i, :) = poses(i, :)
!                print *, poses(i, :), vels(i, :)
                velocities(t, i, :) = vels(i, :)
            end do
       !      Run the integrator of choice for a timestep, using temporary
       !      variables to avoid simultaneous IO for variables
            if (integrator(1:len("kdkrk3")) == "kdkrk3") then
                call kdkrk3(n, masses, poses, vels, dt, low_poses, low_vels, &
                            high_poses, high_vels)
                poses = high_poses
                vels = high_vels
            else if (integrator(1:len("rkf45")) == "rkf45") then
                call rkf45(n, masses, poses, vels, dt, low_poses, low_vels, &
                           high_poses, high_vels)
!                print *, low_poses(2, :), high_poses(2, :)
                do i=1, n
                    poses(i, :) = high_poses(i, :)
                    vels(i, :) = high_vels(i, :)
                end do
            end if
!
       !      Check validity of the integrator's calculation
            do i=1, n
                do j=1, 3
                    if (isnan(poses(i, j))) then
                        print *, "NAN encountered for position at", i, j
                        status = "NAN encountered for position"
                        exit timeloop
                    end if
                    if (isnan(vels(i, j))) then
                        print *, "NAN encountered for velocity at", i, j
                        status = "NAN encountered for velocity"
                        exit timeloop
                    end if
                end do
            end do

       !     print *, vels(2, :)
            ! Recenters the system
            if (update_com) then
                com(:) = center_of_mass(n, masses, i_pos)
                avv(:) = average_vel(n, masses, i_vel)
                do i=1, n
                    poses(i, :) = poses(i, :) - com(:)
                    vels(i, :) = vels(i, :) - avv(:)
                end do
            end if
!
       !     print *, vels(2, :)
            err_sum = 0.
            do i=1, n
                err_sum = err_sum + vec_mag(3, (high_poses - low_poses) &
                        / tolerance)
                err_sum = err_sum + vec_mag(3, (high_vels - low_vels) &
                        / tolerance)
            end do
            err_sum = err_sum / real(n * 3 * 2)**0.5
!            print *, tolerance
!
            dt = adjust_timestep(dt, err_sum, q, t_tot)
!            print *, dt
!
            t_now = t_now + dt
!
!
         end do timeloop

    end subroutine adaptive_n_body_model
end module n_bodies