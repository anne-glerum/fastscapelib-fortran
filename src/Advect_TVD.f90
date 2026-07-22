subroutine Advect_TVD()

  use FastScapeContext, only: nx, ny, xl, yl, dt, h2, b2, etot2, &
                              vx2, vy2, b, h
  implicit none

  double precision :: dx, dy, vmaxx, vmaxy, dt_cfl
  double precision, parameter :: cfl_number = 0.2d0

  dx = xl / dble(nx - 1)
  dy = yl / dble(ny - 1)

  ! Safe directional CFL restriction.
  vmaxx = maxval(abs(vx2))
  vmaxy = maxval(abs(vy2))

  dt_cfl = dt

  if (vmaxx > 0.d0) then
    dt_cfl = min(dt_cfl, cfl_number * dx / vmaxx)
  end if

  if (vmaxy > 0.d0) then
    dt_cfl = min(dt_cfl, cfl_number * dy / vmaxy)
  end if

  dt = dt_cfl
  ! Strang splitting
  call tvd_advect(h2,    vx2, vy2, dx, dy, dt)
  call tvd_advect(b2,    vx2, vy2, dx, dy, dt)
  call tvd_advect(etot2, vx2, vy2, dx, dy, dt)

  b = min(b, h)

contains

subroutine tvd_advect(field, velx, vely, dx, dy, dtloc)

  implicit none

  double precision, intent(inout) :: field(nx, ny)
  double precision, intent(in)    :: velx(nx, ny), vely(nx, ny)
  double precision, intent(in)    :: dx, dy, dtloc

  ! Strang splitting: x half-step, y full-step, x half-step.
  call tvd_x(field, velx, dx, 0.5d0 * dtloc)
  call tvd_y(field, vely, dy, dtloc)
  call tvd_x(field, velx, dx, 0.5d0 * dtloc)

end subroutine tvd_advect

subroutine tvd_x(field, velx, dx, dtloc)

  implicit none

  double precision, intent(inout) :: field(nx, ny)
  double precision, intent(in)    :: velx(nx, ny)
  double precision, intent(in)    :: dx, dtloc

  integer :: i, j
  double precision :: vface, flo, fhi, den, r, phi
  double precision, parameter :: eps = 1.d-14

  double precision, allocatable :: flux(:,:), field_new(:,:)

  ! flux(i,j) is the flux across the face between cells i and i+1.
  allocate(flux(nx-1, ny), field_new(nx, ny))

  do j = 1, ny
    do i = 1, nx-1

      ! Convert node-centred FastScape velocity to face velocity.
      vface = 0.5d0 * (velx(i,j) + velx(i+1,j))

      ! Godunov / first-order upwind flux.
      flo = max(vface, 0.d0) * field(i,j) + &
            min(vface, 0.d0) * field(i+1,j)

      den = field(i+1,j) - field(i,j)

      ! Use high order only where both required neighbouring slopes exist.
      phi = 0.d0

      if (i >= 2 .and. i <= nx-2 .and. abs(den) > eps) then

        if (vface >= 0.d0) then
          ! Upwind direction is from i-1 toward i.
          r = (field(i,j) - field(i-1,j)) / den
        else
          ! Upwind direction is from i+2 toward i+1.
          r = (field(i+2,j) - field(i+1,j)) / den
        end if

        phi = vanleer(r)

      end if

      ! Lax-Wendroff flux.
      fhi = 0.5d0 * vface * (field(i,j) + field(i+1,j)) - &
            0.5d0 * vface**2 * (dtloc / dx) * den

      ! Campforts/Toro-style limited flux.
      flux(i,j) = flo + phi * (fhi - flo)

    end do
  end do

  field_new = field

  ! Keep field(1,:) and field(nx,:) unchanged, as in Advect_Original.
  do j = 1, ny
    do i = 2, nx-1
      field_new(i,j) = field(i,j) - &
          (dtloc / dx) * (flux(i,j) - flux(i-1,j))
    end do
  end do

  field = field_new

  deallocate(flux, field_new)

end subroutine tvd_x

subroutine tvd_y(field, vely, dy, dtloc)

  implicit none

  double precision, intent(inout) :: field(nx, ny)
  double precision, intent(in)    :: vely(nx, ny)
  double precision, intent(in)    :: dy, dtloc

  integer :: i, j
  double precision :: vface, flo, fhi, den, r, phi
  double precision, parameter :: eps = 1.d-14

  double precision, allocatable :: flux(:,:), field_new(:,:)

  ! flux(i,j) is the flux across the face between cells j and j+1.
  allocate(flux(nx, ny-1), field_new(nx, ny))

  do i = 1, nx
    do j = 1, ny-1

      ! Convert node-centred FastScape velocity to face velocity.
      vface = 0.5d0 * (vely(i,j) + vely(i,j+1))

      ! Godunov / first-order upwind flux.
      flo = max(vface, 0.d0) * field(i,j) + &
            min(vface, 0.d0) * field(i,j+1)

      den = field(i,j+1) - field(i,j)

      phi = 0.d0

      if (j >= 2 .and. j <= ny-2 .and. abs(den) > eps) then

        if (vface >= 0.d0) then
          r = (field(i,j) - field(i,j-1)) / den
        else
          r = (field(i,j+2) - field(i,j+1)) / den
        end if

        phi = vanleer(r)

      end if

      fhi = 0.5d0 * vface * (field(i,j) + field(i,j+1)) - &
            0.5d0 * vface**2 * (dtloc / dy) * den

      flux(i,j) = flo + phi * (fhi - flo)

    end do
  end do

  field_new = field

  ! Keep field(:,1) and field(:,ny) unchanged, as in Advect_Original.
  do i = 1, nx
    do j = 2, ny-1
      field_new(i,j) = field(i,j) - &
          (dtloc / dy) * (flux(i,j) - flux(i,j-1))
    end do
  end do

  field = field_new

  deallocate(flux, field_new)

end subroutine tvd_y

double precision function vanleer(r)

  implicit none
  double precision, intent(in) :: r

  if (r <= 0.d0) then
    vanleer = 0.d0
  else
    vanleer = (2.d0 * r) / (1.d0 + r)
  end if

end function

end subroutine Advect_TVD
