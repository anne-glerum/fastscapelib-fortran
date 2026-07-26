subroutine Advect()

  use FastScapeContext, only: advection_scheme, ADVECTION_ORIGINAL, ADVECTION_TVD
  implicit none

  select case (advection_scheme)

  case (ADVECTION_TVD)
    call Advect_TVD()

  case (ADVECTION_ORIGINAL)
    call Advect_Original()

  case default
    call Advect_Original()

  end select

end subroutine Advect