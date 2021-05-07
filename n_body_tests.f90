program n_body_test
    use n_bodies
    implicit none
    integer(8), parameter :: n = 2, n_times = 100000
    real :: q = 1e-3
    real(8), dimension(n) :: masses
    real(8), dimension(n, 3) :: i_positions, i_velocities, poses, vels
    real(8), dimension(n_times, n, 3) :: velocities, positions
    real(8), dimension(n_times) :: times
    real(8) :: dt, t_tot, g, tolerance
    logical :: update_com
    character(256) :: status
    integer(8) :: i, j, k, q_dt
    masses = (/1. - q, q/)
    i_positions(1, :) = (/-q, 0., 0./)
    i_positions(2, :) = (/1. - q, 0., 0./)
    i_velocities(1, :) = (/0., -q, 0./)
    i_velocities(2, :) = (/0., 1. - q, 0./)
    dt = 1e-3
    t_tot = n_times * dt
    tolerance = 1e-10
    q_dt = 4
    g = 1.
    update_com = .true.

    print *, "Euler"
    call n_body_model("euler", n, n_times, masses, i_positions, i_velocities, &
                      t_tot, dt, g, update_com, times, positions, velocities, &
                      status)
    print *, status

    print *, "RK5"
    call n_body_model("rk5", n, n_times, masses, i_positions, i_velocities, &
                  t_tot, dt, g, update_com, times, positions, velocities, &
                  status)
    print *, status
!
    print *, "RKF45"
    status = "Ok"
    call adaptive_n_body_model("rkf45", n, n_times, masses, i_positions, &
            i_velocities, t_tot, tolerance, dt, q_dt, g, update_com, times, &
            positions, velocities, status)
    print *, status
!    do i=1, n_times
!        print *, positions(i, 2, :)
!    end do

end program n_body_test