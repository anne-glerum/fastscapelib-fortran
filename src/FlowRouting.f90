!--------------------------------------------------------------------------------------------

subroutine FlowRouting ()

  use FastScapeContext

  implicit none

  double precision :: dx, dy

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

  ! finds receiver

  call find_receiver (h, nx, ny, dx, dy, rec, length, &
    bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic)

  ! finds donors

  call find_donor (rec, nn, ndon, don)

  ! find stack

  call find_stack (rec, don, ndon, nn, catch0, stack, catch)

  ! removes local minima
  call LocalMinima (stack,rec,bounds_bc,ndon,don,h,length,nx,ny,dx,dy)

  ! computes receiver and stack information for mult-direction flow
  call find_mult_rec (h,rec,stack,hwater,mrec,mnrec,mwrec,mlrec,mstack,nx,ny,dx,dy,p,p_mfd_exp, &
    bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic)

  ! compute lake depth
  lake_depth = hwater - h

  return

end subroutine FlowRouting

!--------------------------------------------------------------------------------------------

subroutine FlowRoutingSingleFlowDirection ()

  use FastScapeContext

  implicit none

  integer :: i, ijk, ijr
  double precision :: dx, dy,deltah

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

  ! finds receiver

  call find_receiver (h, nx, ny, dx, dy, rec, length, &
    bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic)

  ! finds donors

  call find_donor (rec, nn, ndon, don)

  ! find stack

  call find_stack (rec, don, ndon, nn, catch0, stack, catch)

  ! removes local minima
  call LocalMinima (stack,rec,bounds_bc,ndon,don,h,length,nx,ny,dx,dy)

  ! find hwater

  hwater = h

  ! fill the local minima with a nearly planar surface

  deltah = 1.d-8
  do i=1,nn
    ijk = stack(i)
    ijr = rec(ijk)
    if (ijr.ne.0) then
      if (hwater(ijr).gt.hwater(ijk)) then
        hwater(ijk) = hwater(ijr) + deltah
      endif
    endif
  enddo

  ! compute lake depth
  lake_depth = hwater - h

  return

end subroutine FlowRoutingSingleFlowDirection

!--------------------------------------------------------------------------------------------

subroutine FlowAccumulation ()

  use FastScapeContext

  implicit none

  integer :: ij, ijk, k
  double precision :: dx,dy

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  a=dx*dy*precip
  do ij=1,nn
    ijk=mstack(ij)
    do k =1,mnrec(ijk)
      a(mrec(k,ijk))=a(mrec(k,ijk))+a(ijk)*mwrec(k,ijk)
    enddo
  enddo

  return

end subroutine FlowAccumulation

!--------------------------------------------------------------------------------------------

subroutine ProvenanceRouting ()

  use FastScapeContext

  implicit none

  integer :: ij, ijk, ijr, k, c
  double precision :: dx, dy, cellarea
  double precision :: erosion_supply, deposition
  double precision :: total_flux, deposited_fraction

  if (ncomp .eq. 0) return

  dx = xl/(nx-1)
  dy = yl/(ny-1)
  cellarea = dx*dy

  ! These arrays represent the current timestep only.
  prov_flux = 0.d0
  donor_count = 0.d0

  ! mstack processes all donors before their receivers.
  do ij = 1, nn

    ijk = mstack(ij)
    c = composition(ijk)

    ! Local sediment supply from erosion.
    ! erate > 0 means erosion, in m/yr.
    ! erosion_supply is in m3/yr.
    erosion_supply = max(erate(ijk),0.d0)*cellarea

    prov_flux(c,ijk) = prov_flux(c,ijk) + erosion_supply

    ! Count cells that are actively contributing sediment.
    if (erosion_supply .gt. 0.d0) then
      donor_count(c,ijk) = donor_count(c,ijk) + 1.d0
    endif

    ! Local deposition.
    ! erate < 0 means deposition. Remove that amount
    ! from the arriving provenance mixture.
    deposition = max(-erate(ijk),0.d0)*cellarea
    total_flux = sum(prov_flux(:,ijk))

    if (deposition .gt. 0.d0 .and. total_flux .gt. 0.d0) then

      deposited_fraction = min(deposition,total_flux)/total_flux

      ! Save cumulative deposition, m3.
      prov_deposited(:,ijk) = prov_deposited(:,ijk) + &
        deposited_fraction*prov_flux(:,ijk)*dt

      ! Material remaining in transport.
      prov_flux(:,ijk) = (1.d0-deposited_fraction)*prov_flux(:,ijk)

    endif

    ! Either collect at a base-level node or route
    ! to all MFD receivers.
    if (bounds_bc(ijk)) then

      ! Flux rate m3/yr converted to volume m3.
      prov_delivered(:,ijk) = prov_delivered(:,ijk) + &
        prov_flux(:,ijk)*dt

    else

      do k = 1, mnrec(ijk)
        ijr = mrec(k,ijk)

        prov_flux(:,ijr) = prov_flux(:,ijr) + &
          mwrec(k,ijk)*prov_flux(:,ijk)

        donor_count(:,ijr) = donor_count(:,ijr) + &
          mwrec(k,ijk)*donor_count(:,ijk)
      enddo

    endif

  enddo

end subroutine ProvenanceRouting

!--------------------------------------------------------------------------------------------

subroutine FlowAccumulationSingleFlowDirection ()

  use FastScapeContext

  implicit none

  integer :: ij, ijk
  double precision :: dx,dy

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  a=dx*dy*precip
  do ij=nn,1,-1
    ijk=stack(ij)
    a(rec(ijk))=a(rec(ijk))+a(ijk)
  enddo

  return

end subroutine FlowAccumulationSingleFlowDirection

!--------------------------------------------------------------------------------------------

subroutine ProvenanceRoutingSingleFlowDirection ()

  use FastScapeContext

  implicit none

  integer :: ij, ijk, ijr, c
  double precision :: dx, dy, cellarea
  double precision :: erosion_supply, deposition
  double precision :: total_flux, deposited_fraction

  if (ncomp .eq. 0) return

  dx = xl/(nx-1)
  dy = yl/(ny-1)
  cellarea = dx*dy

  prov_flux = 0.d0
  donor_count = 0.d0

  ! For SFD, stack(nn:1:-1) is donor -> receiver order.
  do ij = nn, 1, -1

    ijk = stack(ij)
    c = composition(ijk)

    erosion_supply = max(erate(ijk),0.d0)*cellarea
    prov_flux(c,ijk) = prov_flux(c,ijk) + erosion_supply

    if (erosion_supply .gt. 0.d0) then
      donor_count(c,ijk) = donor_count(c,ijk) + 1.d0
    endif

    deposition = max(-erate(ijk),0.d0)*cellarea
    total_flux = sum(prov_flux(:,ijk))

    if (deposition .gt. 0.d0 .and. total_flux .gt. 0.d0) then
      deposited_fraction = min(deposition,total_flux)/total_flux
      prov_deposited(:,ijk) = prov_deposited(:,ijk) + &
        deposited_fraction*prov_flux(:,ijk)*dt
      prov_flux(:,ijk) = (1.d0-deposited_fraction)*prov_flux(:,ijk)
    endif

    if (bounds_bc(ijk)) then
      prov_delivered(:,ijk) = prov_delivered(:,ijk) + &
        prov_flux(:,ijk)*dt
    else
      ijr = rec(ijk)
      prov_flux(:,ijr) = prov_flux(:,ijr) + prov_flux(:,ijk)
      donor_count(:,ijr) = donor_count(:,ijr) + donor_count(:,ijk)
    endif

  enddo

end subroutine ProvenanceRoutingSingleFlowDirection

!--------------------------------------------------------------------------------------------

subroutine find_mult_rec (h,rec0,stack0,water,rec,nrec,wrec,lrec,stack,nx,ny,dx,dy,p,p_mfd_exp, &
  bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic)

  ! subroutine to find multiple receiver information
  ! in input:
  ! h is topography
  ! rec0 is single receiver information
  ! stack0 is stack (from bottom to top) obtained by using single receiver information
  ! water is the surface of the lakes (or topography where there is no lake)
  ! nx, ny resolution in x- and y-directions
  ! dx, dy grid spacing in x- and y-directions
  ! p is exponent to which the slope is put to share the water/sediment among receivers
  ! bounds: boundaries type (boundary conditions)
  ! in output:
  ! rec: multiple receiver information
  ! nrec: number of receivers for each node
  ! wrec: weight for each receiver
  ! lrec: distance to each receiver
  ! stack: stoack order for multiple receivers (from top to bottom)

  integer, intent(in) :: bounds_i1, bounds_i2, bounds_j1, bounds_j2
  logical, intent(in) :: bounds_xcyclic, bounds_ycyclic
  integer nx,ny
  double precision h(nx*ny),wrec(8,nx*ny),lrec(8,nx*ny),dx,dy,p,water(nx*ny),p_mfd_exp(nx*ny)
  integer rec(8,nx*ny),nrec(nx*ny),stack(nx*ny),rec0(nx*ny),stack0(nx*ny)

  integer :: nn,i,j,ii,jj,iii,jjj,ijk,k,ijr,nparse,nstack,ijn,ij
  double  precision :: slopemax,sumweight,deltah,slope
  integer, dimension(:), allocatable :: ndon,vis,parse
  integer, dimension(:,:), allocatable :: don
  double precision, dimension(:), allocatable :: h0

  nn=nx*ny

  allocate (h0(nn))

  h0=h

  ! fill the local minima with a nearly planar surface

  deltah = 1.d-8
  do i=1,nn
    ijk=stack0(i)
    ijr=rec0(ijk)
    if (ijr.ne.0) then
      if (h0(ijr).gt.h0(ijk)) then
        h0(ijk)=h0(ijr)+deltah
      endif
    endif
  enddo

  water = h0

  nrec=0
  wrec=0.d0
  lrec=0.d0
  rec=0.d0

  ! loop on all nodes
  do j=bounds_j1,bounds_j2
    do i=bounds_i1,bounds_i2
      ij = (j-1)*nx + i
      slopemax = 0.
      do jj=-1,1
        jjj= j + jj
        if (jjj.lt.1.and.bounds_ycyclic) jjj=jjj+ny
        jjj=max(jjj,1)
        if (jjj.gt.ny.and.bounds_ycyclic) jjj=jjj-ny
        jjj=min(jjj,ny)
        do ii=-1,1
          iii = i + ii
          if (iii.lt.1.and.bounds_xcyclic) iii=iii+nx
          iii=max(iii,1)
          if (iii.gt.nx.and.bounds_xcyclic) iii=iii-nx
          iii=min(iii,nx)
          ijk = (jjj-1)*nx + iii
          if (h0(ij).gt.h0(ijk)) then
            nrec(ij)=nrec(ij)+1
            rec(nrec(ij),ij) = ijk
            lrec(nrec(ij),ij) = sqrt((ii*dx)**2 + (jj*dy)**2)
            wrec(nrec(ij),ij) = (h0(ij) - h0(ijk))/lrec(nrec(ij),ij)
          endif
        enddo
      enddo
    enddo
  enddo

  do ij =1,nn
    if (p<0.d0) then
      slope = 0.d0
      if (nrec(ij).ne.0) slope = real(sum(wrec(1:nrec(ij),ij))/nrec(ij))
      p_mfd_exp(ij) = 0.5 + 0.6*slope
    endif
    do k=1,nrec(ij)
      wrec(k,ij) = wrec(k,ij)**p_mfd_exp(ij)
    enddo
    sumweight = sum(wrec(1:nrec(ij),ij))
    wrec(1:nrec(ij),ij) = wrec(1:nrec(ij),ij)/sumweight
  enddo

  allocate (ndon(nn),don(8,nn))

  ndon=0

  do ij=1,nn
    do k=1,nrec(ij)
      ijk = rec(k,ij)
      ndon(ijk)=ndon(ijk)+1
      don(ndon(ijk),ijk) = ij
    enddo
  enddo

  allocate (vis(nn),parse(nn))

  nparse=0
  nstack=0
  stack=0
  vis=0
  parse=0

  ! we go through the nodes
  do ij=1,nn
    ! when we find a "summit" (ie a node that has no donors)
    ! we parse it (put it in a stack called parse)
    if (ndon(ij).eq.0) then
      nparse=nparse+1
      parse(nparse)=ij
    endif
    ! we go through the parsing stack
    do while (nparse.gt.0)
      ijn=parse(nparse)
      nparse=nparse-1
      ! we add the node to the stack
      nstack=nstack+1
      stack(nstack)=ijn
      ! for each of its receivers we increment a counter called vis
      do  ijk=1,nrec(ijn)
        ijr=rec(ijk,ijn)
        vis(ijr)=vis(ijr)+1
        ! if the counter is equal to the number of donors for that node we add it to the parsing stack
        if (vis(ijr).eq.ndon(ijr)) then
          nparse=nparse+1
          parse(nparse)=ijr
        endif
      enddo
    enddo
  enddo
  if (nstack.ne.nn) stop 'error in stack'

  deallocate (ndon,don,vis,parse,h0)

  return

end subroutine find_mult_rec

!--------------------------------------------------------------------------------------------

subroutine find_receiver (h, nx, ny, dx, dy, rec, length, &
  bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic)

  implicit none

  integer, intent(in) :: bounds_i1, bounds_i2, bounds_j1, bounds_j2
  logical, intent(in) :: bounds_xcyclic, bounds_ycyclic
  integer, intent(in) :: nx, ny
  double precision, dimension(nx*ny), intent(in) :: h
  double precision, intent(in) :: dx, dy
  integer, dimension(nx*ny), intent(out) :: rec
  double precision, dimension(nx*ny), intent(out) :: length

  integer :: nn, ij, i, j, ii, jj, iii, jjj, ijk
  double precision :: l, smax, slope

  nn = nx*ny

  ! resets receiver and distance between node and its receiver

  do ij=1,nn
    rec(ij)=ij
    length(ij)=0.d0
  enddo

  ! finds receiver using steepest descent/neighbour method

  do j=bounds_j1,bounds_j2
    do i=bounds_i1,bounds_i2
      ij=i+(j-1)*nx
      smax=tiny(smax)
      do jj=-1,1
        do ii=-1,1
          iii=i+ii
          if (iii.lt.1.and.bounds_xcyclic) iii=iii+nx
          iii=max(iii,1)
          if (iii.gt.nx.and.bounds_xcyclic) iii=iii-nx
          iii=min(iii,nx)
          jjj=j+jj
          if (jjj.lt.1.and.bounds_ycyclic) jjj=jjj+ny
          jjj=max(jjj,1)
          if (jjj.gt.ny.and.bounds_ycyclic) jjj=jjj-ny
          jjj=min(jjj,ny)
          ijk=iii+(jjj-1)*nx
          if (ijk.ne.ij) then
            l=sqrt((dx*ii)**2+(dy*jj)**2)
            slope=(h(ij)-h(ijk))/l
            if (slope.gt.smax) then
              smax=slope
              rec(ij)=ijk
              length(ij)=l
            endif
          endif
        enddo
      enddo
    enddo
  enddo

  return

end subroutine find_receiver

!--------------------------------------------------------------------------------------------

subroutine find_donor (rec, nn, ndon, don)

  implicit none

  integer, intent(in) :: nn
  integer, dimension(nn), intent(in) :: rec
  integer, dimension(nn), intent(out) :: ndon
  integer, dimension(8,nn), intent(out) :: don

  integer ij, ijk

  ! inverts receiver array to compute donor arrays
  ndon=0
  do ij=1,nn
    if (rec(ij).ne.ij) then
      ijk=rec(ij)
      ndon(ijk)=ndon(ijk)+1
      don(ndon(ijk),ijk)=ij
    endif
  enddo

  return

end subroutine find_donor

!--------------------------------------------------------------------------------------------

subroutine find_stack (rec, don, ndon, nn, catch0, stack, catch)

  implicit none

  integer, intent(in) :: nn
  integer, intent(in), dimension(nn) :: rec, ndon
  integer, intent(in), dimension(8,nn) :: don
  double precision, intent(in), dimension(nn) :: catch0
  integer, intent(out), dimension(nn) :: stack
  double precision, dimension(nn), intent(out) :: catch

  integer :: ij, nstack

  ! computes stack by recursion
  nstack=0
  catch=catch0
  do ij=1,nn
    if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack_recursively (ij,don,ndon,nn,stack,nstack,catch)
    endif
  enddo

  return

end subroutine find_stack

!--------------------------------------------------------------------------------------------

recursive subroutine find_stack_recursively (ij,don,ndon,nn,stack,nstack,catch)

! recursive routine to go through all nodes following donor information

  implicit none

  integer k,ij,ijk,nn,nstack
  integer don(8,nn),ndon(nn),stack(nn)
  double precision catch(nn)

  do k=1,ndon(ij)
    ijk=don(k,ij)
    nstack=nstack+1
    stack(nstack)=ijk
    catch(ijk)=catch(ij)
    call find_stack_recursively (ijk,don,ndon,nn,stack,nstack,catch)
  enddo

  return

end subroutine find_stack_recursively

!--------------------------------------------------------------------------------------------
