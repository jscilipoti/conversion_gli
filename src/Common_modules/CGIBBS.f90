module CGIBBS
    implicit none
    private
    public :: NF, MAXZ, GNUL, Z, A, XVL, SFAS,&
    &GAM, AL, DA, XM
  
 
    ! These are module variables that can be used by any program unit that uses this module
    integer :: NF
    real(8) :: MAXZ, GNUL
    real(8), dimension(4) :: SFAS
    real(8), dimension(10) :: Z, A, AL
    real(8), dimension(10,4) :: XVL,XM
    real(8), dimension(10,10) ::  GAM,DA
    
  
  end module CGIBBS