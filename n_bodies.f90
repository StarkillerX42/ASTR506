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

    subroutine n_bodies(integrator, n, masses, i_pos, i_vel, t_tot, dt, G)
        character, intent(in) :: integrator
        integer, intent(in) :: n
        real, masses(n), intent(in) :: masses
        real, dimension(n, 3), intent(in) :: i_pos
        real,
    end subroutine n_bodies


end module n_bodies