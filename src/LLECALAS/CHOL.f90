subroutine CHOL(N, A)                                              
    IMPLICIT REAL*8(A-H, O-Z)                                          
    dimension A(2, 2)                                                  
    do 50 I = 1, N                                                       
    I1 = I-1                                                            
    if (I1.EQ.0) goto 30                                               
    do 20 J = I, N                                                       
    do 20 K = 1, I1                                                      
 20 A(I, J) = A(I, J)-A(I, K)*A(J, K)                                       
 30 if (A(I, I).LT.1.D-14) A(I, I) = 1.D-14                                
    A(I, I) = DSQRT(A(I, I))                                              
    if (I.EQ.N) goto 100                                               
    J1 = I+1                                                            
    do 50 J = J1, N                                                      
 50 A(J, I) = A(I, J)/A(I, I)                                              
100 return                                                            
    end                                                              
! subroutine GAUSL SOLVES N LINEAR ALGEBRAIC EQUATIONS BY GAUSS        
! ELIMINATION WITH ROW PIVOTING                                        
! TO SOLVE THE PROBLEM QX = U, WHERE Q IS A NXN MATRIX AND U IS NXNS, 
! ONE PLACES Q IN THE FIRST N COLUMNS OF A AND U IS PLACED IN THE      
! FOLLOWING NS COLUMNS.                                                
! THE PROGRAM RETURNS X = Q**(-1)*U AT THE PREVIOUS POSITION OF U.       
! *                                                                    
! ND IS THE ROW dimension AND NCOL IS THE COLUMN dimension OF A.       
! BOTH MUST BE TRANSFERRED TO THE subroutine.                          
! *****************                                                    
!                                                                 