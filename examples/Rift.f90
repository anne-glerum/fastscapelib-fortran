program Rift

  ! Provenance test: two-sided rift basin
  !
  ! The model represents a subsiding central rift basin surrounded by
  ! uplifting left and right rift shoulders.
  !
  ! Composition 1 occupies the left rift shoulder (x < 50 km).
  ! Composition 2 occupies the right rift shoulder (x >= 50 km).
  !
  ! Erosion from both shoulders supplies sediment to the central subsiding
  ! basin. Provenance tracking records the relative contribution of each
  ! shoulder to deposited rift-basin sediment.

  implicit none

  integer :: nx, ny, istep, nstep, i, j
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g, x, y, r2
  double precision, dimension(:), allocatable :: h, u, chi, kf, kd, ux, uy
  integer, dimension(:), allocatable :: comp

  ! initialize FastScape
  call FastScape_Init ()

  ! set grid size
  nx = 401
  ny = 401
  call FastScape_Set_NX_NY (nx,ny)
  ! Call advection scheme
  call FastScape_Set_Advection_Scheme(2) ! 1 - Original scheme; 2 - FVM-TVD scheme
  ! set number of composition
  call FastScape_Set_NComposition(2)

  ! allocate memory
  call FastScape_Setup ()

  ! set model dimensions
  xl = 100.d3
  yl = xl
  call FastScape_Set_XL_YL (xl,yl)

  ! Composition 1: 0 <= x < 50 km
  ! Composition 2: 50 km <= x <= 100 km
  allocate(comp(nx*ny))

  do j = 1, ny
    do i = 1, nx

      x = dble(i-1)*xl/dble(nx-1)

      if (x .lt. 50.d3) then
        comp(i+(j-1)*nx) = 1
      else
        comp(i+(j-1)*nx) = 2
      endif

    enddo
  enddo

  call FastScape_Set_Composition(comp)
  ! set time step
  dt = 1.d5
  call FastScape_Set_DT (dt)

  ! set sloping topography
    allocate(h(nx*ny))

  ! Small perturbation to initiate channels.
  call random_number(h)
  h = 5.d0*(h-0.5d0)

  do j = 1, ny
    y = dble(j-1)*yl/dble(ny-1)

    do i = 1, nx
      x = dble(i-1)*xl/dble(nx-1)

      ! Elevated rift shoulders; low central rift valley.
      h(i+(j-1)*nx) = h(i+(j-1)*nx) + &
        500.d0*(1.d0-exp(-((x-50.d3)/12.d3)**2))

      ! Central depocentre within the rift valley.
      h(i+(j-1)*nx) = h(i+(j-1)*nx) - &
        200.d0*exp(-((x-50.d3)/10.d3)**2 - &
                    ((y-50.d3)/20.d3)**2)

    enddo
  enddo

  ! Fixed top and bottom base-level boundaries for BC = 1010.
  h(1:nx) = 0.d0
  h(nx*(ny-1)+1:nx*ny) = 0.d0

  call FastScape_Init_H(h)

  call FastScape_Init_H(h)
  call FastScape_Init_H (h)

  ! set erosional parameters
  allocate (kf(nx*ny),kd(nx*ny))
  kf = 1.d-5
  kfsed = -1.d0
  m = 0.45d0
  n = 1.d0
  kd = 1.d-2
  kdsed = -1.d0
  g = 1.d0
  call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, g, -2.d0)

  ! set uplift rate (uniform while keeping boundaries at base level)
  allocate(u(nx*ny), ux(nx*ny), uy(nx*ny))

  do j = 1, ny
    do i = 1, nx

      x = dble(i-1)*xl/dble(nx-1)

      ! Rift tectonics:
      ! uplifted shoulders at 30 and 70 km,
      ! subsiding basin at 50 km.
      u(i+(j-1)*nx) = &
        1.d-3*exp(-((x-30.d3)/12.d3)**2) + &
        1.d-3*exp(-((x-70.d3)/12.d3)**2) - &
        1d-3*exp(-((x-50.d3)/10.d3)**2)

      ! Horizontal extension away from x = 50 km.
      ! x = 25 km -> ux < 0: moves left
      ! x = 50 km -> ux = 0
      ! x = 75 km -> ux > 0: moves right
      ux(i+(j-1)*nx) = -1.d-3 * &
        sin(2.d0*acos(-1.d0)*x/xl)

      ! No along-rift advection.
      uy(i+(j-1)*nx) = 0.d0

    enddo
  enddo

  ! Fixed bottom boundary.
  u(1:nx) = 0.d0
  ux(1:nx) = 0.d0
  uy(1:nx) = 0.d0

  ! Fixed right boundary.
  u(nx:nx*ny:nx) = 0.d0
  ux(nx:nx*ny:nx) = 0.d0
  uy(nx:nx*ny:nx) = 0.d0

  ! Fixed left boundary.
  u(1:nx*ny:nx) = 0.d0
  ux(1:nx*ny:nx) = 0.d0
  uy(1:nx*ny:nx) = 0.d0

  ! Fixed top boundary.
  u(nx*(ny-1)+1:nx*ny) = 0.d0
  ux(nx*(ny-1)+1:nx*ny) = 0.d0
  uy(nx*(ny-1)+1:nx*ny) = 0.d0

  call FastScape_Set_U(u)
  call FastScape_Set_V(ux,uy)

  ! set boundary conditions
  call FastScape_Set_BC (1111)

  ! set number of time steps and initialize counter istep
  nstep = 200
  call FastScape_Get_Step (istep)

  !allocate memory to extract chi
  allocate (chi(nx*ny))

  ! loop on time stepping
  do while (istep<nstep)
    ! execute step
    call FastScape_Execute_Step()
    ! get value of time step counter
    call FastScape_Get_Step (istep)
    ! extract solution
    call FastScape_Copy_Chi (chi)
    ! create VTK file
    call FastScape_VTK (chi, 2.d0)
    ! create prevenance VTk
    call Fastscape_Provenance_VTK(2.d0, istep)
    ! outputs h values
    call FastScape_Copy_h (h)
    print*,'step',istep
    print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
  enddo

  ! output timing
  call FastScape_Debug()

  ! end FastScape run
  call FastScape_Destroy ()

  deallocate (h,u,kf,kd,chi,comp)

end program Rift

