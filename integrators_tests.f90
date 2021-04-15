program test
    use integrators
    implicit none
    integer, parameter :: n = 2
    real :: q = 1e-3
    real(8), dimension(n) :: masses
    real(8), dimension(n, 3) :: positions, velocities, poses, vels, poses_low, &
            vels_low, poses_high, vels_high, accels, gal_positions, &
            gal_velocities
    real(8) :: dt
    masses = (/1. - q, q/)
    positions(1, :) = (/-q, 0., 0./)
    positions(2, :) = (/1. - q, 0., 0./)
    velocities(1, :) = (/0., -q, 0./)
    velocities(2, :) = (/0., 1. - q, 0./)
    gal_positions(1, :) = (/1., 0., 0./)
    gal_positions(2, :) = (/1., 0., 0./)
    gal_velocities(1, :) = (/0., 0.43948513, 0./)
    gal_velocities(2, :) = (/0., 0.25373686, 0./)
    dt = 1e-3

    call euler_forward(2, masses, positions, velocities, dt, poses, vels)
    print *, "Euler forward:"
    print *, poses
    print *, vels
    print *, ""

    call kdk(2, masses, positions, velocities, dt, poses, vels)
    print *, "Kick Drift Kick:"
    print *, poses
    print *, vels
    print *, ""

    call rk_pde(2, masses, positions, velocities, poses, vels)
    print *, "Runge Kutta PDE"
    print *, poses
    print *, vels
    print *, ""

    call rk3_classic(2, masses, positions, velocities, dt, poses, vels)
    print *, "Runge Kutta Classic 3:"
    print *, poses
    print *, vels

    call rkf45(2, masses, positions, velocities, dt, poses_low, vels_low, &
            poses_high, vels_high)
    print *, "Runge Kutta Fehlberg 4-5:"
    print *, poses_high
    print *, vels_high
    print *, ""

    call kdk_gal(2, masses, gal_positions, gal_velocities, dt, poses, vels)
    print *, "Kick Drift Kick Galaxy:"
    print *, poses
    print *, vels
    print *, ""



end program test
