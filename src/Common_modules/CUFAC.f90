module CUFAC
    implicit none
    private
    public :: NKK,NGG,Pxx,Txx
  
 
    ! These are module variables that can be used by any program unit that uses this module
    integer :: NKK, NGG
    real(8) :: Txx
    real(8), dimension(10,10) :: Pxx
    
  
  end module CUFAC