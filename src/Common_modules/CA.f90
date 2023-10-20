module CA
    implicit none
    private
    public :: XC, GE, GC
  
 
    ! These are module variables that can be used by any program unit that uses this module
    real(8), dimension(5) :: XC
    real(8), dimension(5,2) :: GE,GC
    
  
  end module CA