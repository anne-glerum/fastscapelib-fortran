program AdvectMountain

  ! Test for different advection schemes (Original vs FVM-TVD)
  ! In this model a rectangular domain (1000 km x 300 km) is subjected to constant uplift
  ! of 0.5 mm/yr and constant high advection velocity of 5 mm/yr. 
  ! initial random topography
  ! non linear erosion law with n = 2.2 and m = 1. 
  ! deposition coefficient g = 0.01
  ! Advection scheme set by: FastScape_Set_Advection_Scheme(n)
  ! n = 1 - Original implicit scheme; n =2 - FVM-TVD scheme
  implicit none
  
  integer :: nx, ny, istep, nstep, i, j, idx
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g, dx, x
  double precision, dimension(:), allocatable :: h, u, chi, kf, kd, ux, uy

  ! Onitialize FastScape
  call FastScape_Init ()

  ! Call advection scheme
  call FastScape_Set_Advection_Scheme(2) ! 1 - Original scheme; 2 - FVM-TVD scheme

  ! Set grid size
  nx = 1001
  ny = 301
  call FastScape_Set_NX_NY (nx,ny)

  ! Allocate memory
  call FastScape_Setup ()

  ! Set model dimensions
  xl = 1000.d3
  yl = 300.d3
  call FastScape_Set_XL_YL (xl,yl)

  ! Set time step
  dt = 2.d4
  call FastScape_Set_DT (dt)

  ! Set random initial topography (0–50 m)
  allocate (h(nx*ny))
  call random_number (h)
  h = h * 50.d0
  call FastScape_Init_H (h)

  ! Set erosional parameters
  allocate (kf(nx*ny),kd(nx*ny))
  kf = 5.d-9
  kfsed = -1.d0
  m = 1.d0
  n = 2.2d0
  kd = 1.d-3
  kdsed = -1.d0
  g = 0.01d0
  call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, g, -2.d0)

  ! Uplift and velocity arrays
  allocate (u(nx*ny))
  allocate (ux(nx*ny), uy(nx*ny))

  ! Constant velocity fields
  u  = 0.5d-3 ! uplift
  ux = 5d-3.  ! x-direction advection
  uy = 0.d0.  ! y-direction advection

  ! Boundary conditions
  u(1:nx)=0.d0
  u(nx:nx*ny:nx)=0.d0
  u(1:nx*ny:nx)=0.d0
  u(nx*(ny-1)+1:nx*ny)=0.d0

  ux(1:nx)=0.d0
  ux(nx:nx*ny:nx)=0.d0
  ux(1:nx*ny:nx)=0.d0
  ux(nx*(ny-1)+1:nx*ny)=0.d0

  uy(1:nx)=0.d0
  uy(nx:nx*ny:nx)=0.d0
  uy(1:nx*ny:nx)=0.d0
  uy(nx*(ny-1)+1:nx*ny)=0.d0

  call FastScape_Set_U(u)
  call FastScape_Set_V(ux, uy)

  call FastScape_Set_BC (1010)

  ! Set number of time steps and initialize counter istep
  nstep = 1000
  call FastScape_Get_Step (istep)

  !Allocate memory to extract chi
  allocate (chi(nx*ny))

  ! Loop on time stepping
  do while (istep<nstep)
    ! execute step
    call FastScape_Execute_Step()
    ! get value of time step counter
    call FastScape_Get_Step (istep)
    ! extract solution
    call FastScape_Copy_Chi (chi)
    ! create VTK file
    if (mod (istep,5) == 0) then
       call FastScape_VTK (chi, 2.d0)
    end if
    ! outputs h values
    call FastScape_Copy_h (h)
    print*,'step',istep
    print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
  enddo

  ! End FastScape run
  call FastScape_Debug()

  call FastScape_Destroy ()

  deallocate (h,u,kf,kd,chi,ux,uy)

end program AdvectMountain

