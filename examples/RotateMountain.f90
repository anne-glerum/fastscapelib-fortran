program RotateMountain

  ! Test for different advection schemes (Original vs FVM-TVD)
  ! Initially build a steady state topography till 10 Myr and then rotate it completely 
  ! by 360 degrees without any uplift or erosion 
  ! 1st phase - build topography with constant uplift of 0.5 mm/yr and linear erosion law (m=0.45, n=1)
  ! till it reaches equilibrium
  ! 2nd phase - stop uplift and erosin and rotate the topography by 360 degrees to check the impact of advection
  ! Advection scheme set by: FastScape_Set_Advection_Scheme(n)
  ! n = 1 - Original implicit scheme; n =2 - FVM-TVD scheme
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

  ! Call advection scheme
  call FastScape_Set_Advection_Scheme(2) ! 1 - Original scheme; 2 - FVM-TVD scheme

  ! Set grid size
  nx = 501
  ny = 501
  call FastScape_Set_NX_NY(nx,ny)

  ! Allocate memory
  call FastScape_Setup()

  ! Set model dimensions
  xl = 500.d3
  yl = 500.d3
  call FastScape_Set_XL_YL(xl,yl)
  
  ! Calculate cell size
  dx = xl/dble(nx-1)
  dy = yl/dble(ny-1)

  ! Set time step
  dt = 2.d4
  call FastScape_Set_DT(dt)

  ! Set same random topography every model run
  allocate(h(nx*ny))

  call random_seed(size=nseed)
  allocate(seed(nseed))

  seed = 12345

  call random_seed(put=seed)
  call random_number(h)

  deallocate(seed)

  h = h * 50.d0

  call FastScape_Init_H(h)

  ! Set erosion parameters
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

  dt_rotate = cfl*min(dx,dy)/vmax

  call FastScape_Set_DT(dt_rotate)

  nrotate = nint(rotation_period/dt_rotate)

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

  ! End FastScape run
  call FastScape_Debug()

  call FastScape_Destroy()

  deallocate(h)
  deallocate(u)
  deallocate(ux)
  deallocate(uy)
  deallocate(kf)
  deallocate(kd)

end program RotateMountain