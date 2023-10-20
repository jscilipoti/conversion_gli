module CVAP
    implicit none
    private
    public :: NOVAP,NDUM,IDUM,PRAT
  
 
    ! These are module variables that can be used by any program unit that uses this module
    integer :: NOVAP, NDUM 
    real(8), dimension(4) ::  IDUM
    real(8), dimension(10) :: PRAT
    
  
  end module CVAP