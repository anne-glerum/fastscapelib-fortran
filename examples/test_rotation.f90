program Rotation
! test advection with rotation. Initially build a steady state topography and then rotate without doing anything else
  implicit none

  integer :: nx, ny, istep, i, j, idx
  integer :: nbuild, nrotate, target_step
  integer :: nseed
  integer, allocatable :: seed(:)

  double precision :: xl, yl, dx, dy
  double precision :: dt, dt_rotate
  double precision :: kfsed, kdsed, m, n, g
  double precision :: omega, rotation_period
  double precision :: xc, yc
  double precision :: x, y
  double precision :: pi
  double precision :: rmax, vmax
  double precision, parameter :: cfl = 0.2d0

  double precision, dimension(:), allocatable :: h
  double precision, dimension(:), allocatable :: u
  double precision, dimension(:), allocatable :: ux, uy
  double precision, dimension(:), allocatable :: kf, kd

! Initialize FastScape

  call FastScape_Init()
  call FastScape_Set_Advection_Scheme(2) ! 1 - Original scheme; 2 - FVM-TVD scheme

  nx = 501
  ny = 501

  call FastScape_Set_NX_NY(nx,ny)
  call FastScape_Setup()

  xl = 500.d3
  yl = 500.d3

  call FastScape_Set_XL_YL(xl,yl)

  dx = xl/dble(nx-1)
  dy = yl/dble(ny-1)

! Initial timestep

  dt = 2.d4
  call FastScape_Set_DT(dt)

! Initial topography

  allocate(h(nx*ny))

  call random_seed(size=nseed)
  allocate(seed(nseed))

  seed = 12345

  call random_seed(put=seed)
  call random_number(h)

  deallocate(seed)

  h = h * 50.d0

  call FastScape_Init_H(h)

! Erosion parameters

  allocate(kf(nx*ny),kd(nx*ny))

  kf = 1.d-5
  kd = 1.d-3

  kfsed = -1.d0
  kdsed = -1.d0

  m = 0.45d0
  n = 1.d0
  g = 0.01d0

  call FastScape_Set_Erosional_Parameters( &
       kf,kfsed,m,n,kd,kdsed,g,g,-2.d0)

! Uplift and velocity arrays

  allocate(u(nx*ny))
  allocate(ux(nx*ny))
  allocate(uy(nx*ny))

  u  = 0.5d-3
  ux = 0.d0
  uy = 0.d0

! Boundary conditions

  u(1:nx)=0.d0
  u(nx:nx*ny:nx)=0.d0
  u(1:nx*ny:nx)=0.d0
  u(nx*(ny-1)+1:nx*ny)=0.d0

  call FastScape_Set_U(u)
  call FastScape_Set_V(ux,uy)

  call FastScape_Set_BC(1010)

! PHASE 1 : Build mountain

  nbuild = 500

  call FastScape_Get_Step(istep)

  do while (istep < nbuild)

     call FastScape_Execute_Step()
     call FastScape_Get_Step(istep)

     if(mod(istep,25)==0) then
        call FastScape_Copy_H(h)
        call FastScape_VTK(h,2.d0)
        print *,'Build step =',istep
     endif

  enddo

! PHASE 2 : Freeze erosion and rotate

  u = 0.d0
  call FastScape_Set_U(u)

  kf = 0.d0
  kd = 0.d0
  g  = 0.d0

  call FastScape_Set_Erosional_Parameters( &
       kf,kfsed,m,n,kd,kdsed,g,g,-2.d0)

! Rotation velocity field

  pi = acos(-1.d0)

  xc = xl/2.d0
  yc = yl/2.d0

  rotation_period = 80.d6

  omega = 2.d0*pi/rotation_period

  do j=1,ny

     y = (j-1)*dy

     do i=1,nx

        x = (i-1)*dx

        idx = i + (j-1)*nx

        ux(idx) = -omega*(y-yc)
        uy(idx) =  omega*(x-xc)

     enddo
  enddo

  call FastScape_Set_V(ux,uy)

! CFL timestep

  rmax = sqrt((xl/2.d0)**2 + (yl/2.d0)**2)

  vmax = omega*rmax

  dt_rotate = 0.9d0*cfl*min(dx,dy)/vmax

  call FastScape_Set_DT(dt_rotate)

  nrotate = nint(0.25d0*rotation_period/dt_rotate)

  target_step = istep + nrotate

! Rotate landscape

  do while (istep < target_step)

     call FastScape_Execute_Step()

     call FastScape_Get_Step(istep)

     if(mod(istep,25)==0) then
        call FastScape_Copy_H(h)
        call FastScape_VTK(h,2.d0)
        print *,'Rotate step =',istep
     endif

  enddo

! Finish

  call FastScape_Debug()

  call FastScape_Destroy()

  deallocate(h)
  deallocate(u)
  deallocate(ux)
  deallocate(uy)
  deallocate(kf)
  deallocate(kd)

end program Rotation