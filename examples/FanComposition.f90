program FanComposition

  ! Provenance test: alluvial fan sourced from a two-composition plateau
  !
  ! The model represents erosion of an uplifting plateau in the northern
  ! half of the domain (y >= 10 km). Sediment is transported southward toward
  ! the fixed base-level boundary at y = 0 km, where an alluvial fan develops.
  !
  ! The plateau is divided into two source lithologies:
  !   Composition 1: left half of the plateau/domain  (x < 5 km)
  !   Composition 2: right half of the plateau/domain (x >= 5 km)
  !
  ! Provenance routing records the volume of sediment supplied by each source
  ! composition, transported through the drainage network, deposited within
  ! the fan, and delivered to the base-level outlet.

  implicit none

  integer :: nx,ny,istep,nstep,nn,ibc, nseed
  double precision, dimension(:), allocatable :: h,x,y,kf,kd,b,u
  real :: time_in,time_out
  double precision :: kfsed,m,n,kdsed,g1,g2,expp
  double precision xl,yl,dt,pi,vex,x1,y1
  integer, dimension(:), allocatable :: comp
  integer i,j
  integer, allocatable :: seed(:)

  ! set model resolution
  nx = 101
  ny = 201
  nn = nx*ny

  pi=atan(1.d0)*4.d0

  ! initialize FastScape
  call FastScape_Init ()
  call FastScape_Set_NX_NY (nx,ny)
  ! set number of compositions
  call FastScape_Set_NComposition(2)
  call FastScape_Setup ()

  ! set model dimensions
  xl=10.d3
  yl=20.d3
  call FastScape_Set_XL_YL (xl,yl)

  ! construct nodal coordinate arrays x and y
  allocate (x(nx*ny),y(nx*ny))
  x = (/((xl*float(i-1)/(nx-1), i=1,nx),j=1,ny)/)
  y = (/((yl*float(j-1)/(ny-1), i=1,nx),j=1,ny)/)
  ! set the composition
  allocate(comp(nx*ny))

  do j = 1, ny
    do i = 1, nx
    ! compute the value of the x-coordinate
      x1 = dble(i-1)*xl/dble(nx-1)

      ! set composition 1 on the left side of the domain, and 2 on the right
      if (x1 .lt. 5.d3) then
        comp(i+(j-1)*nx) = 1
      else
        comp(i+(j-1)*nx) = 2
      endif

    enddo
  enddo

  call FastScape_Set_Composition(comp)
  ! set time step
  dt=2.d3
  call FastScape_Set_DT (dt)

  ! Uplift 
  allocate(u(nx*ny))

  do j = 1, ny
    do i = 1, nx
      ! compute value of y-coordinate
      y1 = dble(j-1)*yl/dble(ny-1)
      ! prescribe uplift on northern part of the domain
      if (y1 .ge. 10.d3) then
        u(i+(j-1)*nx) = 2.d-3
      else
        u(i+(j-1)*nx) = 0.d0
      endif

    enddo
  enddo

  ! Boundary conditions
  u(1:nx)=0.d0
  call FastScape_Set_U(u)

  ! we make the sediment slightly more easily erodible
  allocate (kf(nn),kd(nn))
  kf=1.d-4
  kfsed=1.5d-4
  m=0.4d0
  n=1.d0
  kd=1.d-2
  kdsed=1.5d-2
  g1=0.1d0
  g2=0.1d0
  expp=1.d0
  call FastScape_Set_Erosional_Parameters (kf,kfsed,m,n,kd,kdsed,g1,g2,expp)

  ! bottom side is fixed only
  ibc=1000
  call FastScape_Set_BC (ibc)

  ! initial topography is a 1000 m high plateau
  allocate(h(nn),b(nn))
  call random_seed(size=nseed)
  allocate(seed(nseed))
  seed = 12345
  call random_seed(put=seed)
  call random_number(h)
  deallocate(seed)
  h = 0.1d0*(h-0.5d0) ! Rescaling random topography

  do j = 1, ny
    do i = 1, nx
      ! compute value of y-coordinate
      y1 = dble(j-1)*yl/dble(ny-1)

      ! Constant regional southward slope of 0.01: elevation increases linearly from the southern outlet toward the north.
      h(i+(j-1)*nx) = h(i+(j-1)*nx) + 0.01d0*y1

      ! Smooth 1 km plateau above y = 10 km.
      h(i+(j-1)*nx) = h(i+(j-1)*nx) + &
        1000.d0*0.5d0*(1.d0+tanh((y1-10.d3)/200.d0))

    enddo
  enddo

  ! Fixed bottom boundary elevation.
  h(1:nx) = 0.d0

  call FastScape_Init_H(h)

  ! set number of time steps
  nstep = 500

  ! echo model setup
  call FastScape_View ()

  ! initializes time step
  call FastScape_Get_Step (istep)

  ! set vertical exaggeration
  vex = 3.d0

  ! start of time loop
  call cpu_time (time_in)
  do while (istep.lt.nstep)

    ! execute FastScape step
    call FastScape_Execute_Step ()
    call FastScape_Get_Step (istep)

    ! output vtk with sediment thickness information
    call FastScape_Copy_H (h)
    call FastScape_Copy_Basement (b)
    call FastScape_VTK (h-b, vex)
    call Fastscape_Provenance_VTK(2.d0, istep) ! 2.d0 is the vertical exaggeration

  enddo

  ! display timing information
  call FastScape_Debug()
  call cpu_time (time_out)
  print*,'Total run time',time_out-time_in

  ! exits FastScape
  call FastScape_Destroy ()

  ! deallocate memory
  deallocate(h,x,y,kf,kd,b,u,comp)

end program FanComposition

