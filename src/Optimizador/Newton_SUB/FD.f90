double precision FUNCTION FD (X,P)
!*****************************************************
!     Derivada de la Funci�n Objetivo
!*****************************************************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   !   REAL :: FD
      DIMENSION X(1),XP(1)
      
      XP(1)=X(1)+P
             
      FD=(f_main(XP,1)-f_main(X,1))/P
      
      END FUNCTION FD