program test
    use integrators
    implicit none
    integer, parameter :: n = 2
    real :: q = 1e-3
    real, dimension(n) :: masses
    real, dimension(n, 3) :: positions, velocities, poses, vels, poses_low, &
            vels_low, poses_high, vels_high, accels
    real :: dt
    masses = (/1., 0./)
    positions = reshape((/-q, 0., 0., 1. - q, 0., 0./), shape(positions))
    velocities = reshape((/0., -q, 0., 0., 1. - q, 0./), shape(velocities))
    dt = 1e-3

    call euler_forward(2, masses, positions, velocities, dt, poses, vels, accels)
    print *, "Euler forward:"
    print "(e22.16)", poses(1, 1)
    print *, vels

    call kdk(2, masses, positions, velocities, dt, poses, vels, accels)
    print *, "Kick Drift Kick:"
    print *, poses
    print *, vels

    call rk3_classic(2, masses, positions, velocities, dt, poses, vels, accels)
    print *, "Runge Kutta Classic 3:"
    print *, poses
    print *, vels

    call rkf45(2, masses, positions, velocities, dt, poses_low, vels_low, &
            poses_high, vels_high)
    print *, "Runge Kutta Fehlberg 4-5:"
    print *, poses_high
    print *, vels_high

end program test
