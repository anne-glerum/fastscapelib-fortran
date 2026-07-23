program advection_test

  ! simple example of the use of the FastScapeLib
  ! test for advection schemes (Original vs FVM-TVD)
  implicit none
  
  integer :: nx, ny, istep, nstep, i, j, idx
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g, dx, x
  double precision, dimension(:), allocatable :: h, u, chi, kf, kd, ux, uy

  ! initialize FastScape
  call FastScape_Init ()
  call FastScape_Set_Advection_Scheme(2)
  ! set grid size
  nx = 1001
  ny = 301
  call FastScape_Set_NX_NY (nx,ny)

  ! allocate memory
  call FastScape_Setup ()

  ! set model dimensions
  xl = 1000.d3
  yl = 300.d3
  call FastScape_Set_XL_YL (xl,yl)

  ! set time step
  dt = 2.d4
  call FastScape_Set_DT (dt)

  ! set random initial topography (0–50 m)
  allocate (h(nx*ny))
  call random_number (h)
  h = h * 50.d0
  call FastScape_Init_H (h)

  ! set erosional parameters
  allocate (kf(nx*ny),kd(nx*ny))
  kf = 5.d-9
  kfsed = -1.d0
  m = 1.d0
  n = 2.2d0
  kd = 1.d-3
  kdsed = -1.d0
  g = 0.01d0
  call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, g, -2.d0)

  ! ---- spatially variable uplift ----

allocate (u(nx*ny))
allocate (ux(nx*ny), uy(nx*ny))

! constant fields
u  = 0.5d-3
ux = 5d-3
uy = 0.d0

! ---- boundary conditions (keep as before) ----
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


  ! set boundary conditions
  call FastScape_Set_BC (1010)

  ! set number of time steps and initialize counter istep
  nstep = 1000
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
    if (mod (istep,5) == 0) then
       call FastScape_VTK (chi, 2.d0)
   end if
    ! outputs h values
    call FastScape_Copy_h (h)
    print*,'step',istep
    print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
  enddo

  ! output timing
  call FastScape_Debug()

  ! end FastScape run
  call FastScape_Destroy ()

  deallocate (h,u,kf,kd,chi)

end program advection_test

