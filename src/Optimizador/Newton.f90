     double precision FUNCTION NEWTON (X)
           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER :: NumIterations, N
    !  REAL :: OldApprox, FDeriv_Old, NewApprox, Epsilon, F_Value, PER
    ! REAL :: f_main, FD  
      DIMENSION X(1),XOLD(1),XNEW(1)
      CHARACTER (1) :: A

      ! Get termination values Epsilon and NumIterations and 
      ! initial approximation

      NumIterations=1000
      Epsilon=1E-9
      PRINT *
      ! Initialize F_Value and N; print headings and initial values
        
      F_Value = f_main(X,1)
      N = 0
      PRINT *, "  N      X(N)       f_main(X(N))"
      PRINT *, "============================="
      PRINT 10, 0, X(1), F_Value, 0.0
10    FORMAT(1X, I3, F11.5, 2E14.5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !do i = 0,100
      !    x(1)=i*0.01
      !    write (333,*) x(1),f_main(x,1)
      !enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Iterate using Newton's method while ABS(F_Value) is greater
! than or equal to Epsilon and N has not reached NumIterations

      !DO
        PER = X(1)*1E-9
        ! If a termination condition met, stop generating approximations

        ! Otherwise continue with the following
        N = N + 1
        FDeriv_Old = FD(X,PER)

        ! Terminate if the derivative is 0 at some approximation
        IF (FDeriv_Old == 0) THEN
            PRINT *, "Newton's method fails -- derivative = 0"
        END IF
        XOLD(1)=X(1)
        ! Generate a new approximation
        X(1) = XOLD(1) - (F_Value / FDeriv_Old)
       ! if (x(1)>1) x(1) = xold(1) + per
        if (x(1)<0) x(1) = 0
        F_Value = f_main(X,1)
        F_ValueOld = f_main(XOLD,1)
        PRINT 10, N, X(1), F_Value, FDeriv_Old

      !END DO
!C
!C-----M�todo de la secante
!C
      DO
        IF ((ABS(X(1)-XOLD(1)) < Epsilon) .OR. (N > NumIterations)) EXIT
        
        N = N + 1

        XNEW(1) = X(1) - ((X(1)-XOLD(1))/(F_Value-F_ValueOld))*F_Value
        if( xnew(1) < 0 ) xnew(1) = 0
        
        XOLD(1) = X(1)
        X(1) = XNEW (1)
        F_Value = f_main(X,1)
        F_ValueOld = f_main(XOLD,1)
        PRINT 10, N, X(1), F_Value, F_ValueOld
        
      ENDDO
      
      
      
      
      NEWTON = F_Value

END FUNCTION NEWTON