module openunits
	use iso_fortran_env, only: int32

    implicit none

    private

    public :: &
    & flashInput_unit,&
    & intrcn32_unit, qPar150_unit, rPar150_unit

    ! These are module variables that can be used by any program unit 
    ! that uses this module
    integer(kind=int32) :: &
    & flashInput_unit, intrcn32_unit, qPar150_unit, rPar150_unit


  end module openunits