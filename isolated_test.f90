program test
    use integrators
    implicit none
    integer, parameter :: n = 2
    real :: q = 1e-3
    real, dimension(n) :: masses
    real, dimension(n, 3) :: positions, velocities, poses, vels, poses_low, &
            vels_low, poses_high, vels_high, accels
    real :: dt
    real :: x
    masses = (/1., 0./)
    positions = reshape((/-q, 0., 0., 1. - q, 0., 0./), shape(positions))
    velocities = reshape((/0., -q, 0., 0., 1. - q, 0./), shape(velocities))
    dt = 1e-2
    print *, velocities
    call euler_forward(2, masses, positions, velocities, dt, poses, vels, accels)
    print *, "Euler forward:"
    print *, poses
    print *, vels

end program test