program test
    use integrators
    implicit none
    integer, parameter :: n=2
    real :: q = 1e-3
    real, dimension(n) :: masses
    real, dimension(n, 3) :: positions, velocities, poses, vels
    real :: dt
    masses = (/1., 0./)
    positions = reshape((/-q, 0., 0., 1. - q, 0., 0./), shape(positions))
    velocities = reshape((/0., -q, 0., 0., 1. - q, 0./), shape(velocities))
    dt = 1e-3

    call euler_forward(2, masses, positions, velocities, dt, poses, vels)
    print *, "Euler forward:"
    print *, poses
    poses = 0.
    call rk3_classic(2, masses, positions, velocities, dt, poses, vels)
    print *, "RK3:"
    print *, poses
    poses = 0.
    call rkf45(2, masses, positions, velocities, dt, poses, vels)
    print *, "RKF45:"
    print *, poses
end program test
