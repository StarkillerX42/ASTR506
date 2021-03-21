program test
    implicit none
    integer, parameter :: n=2
    real, dimension(n) :: masses
    real, dimension(n, 2) :: positions, velocities, poses, vels
    real :: dt
    masses = (/1., 0./)
    positions = transpose(reshape((/0., 0., 1., 0./), shape(positions)))
    velocities = transpose(reshape((/0., 0., 0., 1./), shape(velocities)))
    dt = 1e-3

    call euler_forward(2, masses, positions, velocities, dt, poses, vels)
    print *, poses
    print *, "Hello world!"
end program test

subroutine euler_forward(n, masses, positions, velocities, dt, poses, vels)
    implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: masses
    real, dimension(3) :: rel
    real, dimension(n, 3), intent(in) :: positions, velocities
    real, dimension(n, 3) :: accels
    real, dimension(n, 3), intent(out) :: poses, vels
    integer :: i, j
    real :: dt, rel_mag
!    f2py intent(in) n, masses, positions, velocities, dt
!    f2py intent(out) poses, vels
    do i=1, n
        do j=1, n
            if (i /= j) then
                rel = positions(j, :) - positions(i, :)
                call vec_mag(3, rel, rel_mag)
                accels(i, :) = accels(i, :) + masses(j) / rel_mag**3 * rel(:)
            end if
        end do
    end do
    poses = positions + velocities * dt
    vels = velocities + accels * dt
end subroutine euler_forward

subroutine vec_mag(n_dim, rel_vec, rel_mag)
    integer :: n_dim
    real, dimension(n_dim) :: rel_vec
    real rel_mag
    rel_mag = sqrt(sum(rel_vec**2))
end subroutine vec_mag
