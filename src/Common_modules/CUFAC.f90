module CUFAC
    use iso_fortran_env, only: real64, int32
    
    implicit none
    
    private

    public :: NKK,NGG,Pxx,Txx
  
 
    ! These are module variables that can be used by any program unit that uses 
    ! this module
    integer(kind=int32) :: NKK, NGG
    real(kind=real64) :: Txx
    real(kind=real64), dimension(10,10) :: Pxx
    
  
  end module CUFAC