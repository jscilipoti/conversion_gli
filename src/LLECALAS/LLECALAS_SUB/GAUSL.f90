subroutine GAUSL(ND, NCOL, N, NS, A)                                  
                                                                     
    IMPLICIT REAL*8 (A-H, O-Z)                                         
    dimension A(ND, NCOL)                                              
    N1 = N+1                                                            
    NT = N+NS                                                           
    if(N .EQ. 1) GO TO 50                                            
! START ELIMINATION                                                
!                                                                 
!                                                                 
    do 10 I = 2, N                                                       
    IP = I-1                                                            
    I1 = IP                                                             
    X = DABS(A(I1, I1))                                                  
    do 11 J = I, N                                                       
    if(DABS(A(J, I1)) .LT. X) GO TO 11                                
    X = DABS(A(J, I1))                                                   
    IP = J                                                              
 11 continue                                                          
    if(IP .EQ. I1) GO TO 13                                          
!                                                                 
! ROW INTERCHANGE                                                   
!                                                                 
    do 12 J = I1, NT                                                     
    X = A(I1, J)                                                         
    A(I1, J) = A(IP, J)                                                   
 12 A(IP, J) = X                                                         
 13 do 10 J = I, N                                                       
    X = A(J, I1)/A(I1, I1)                                                
    do 10 K = I, NT                                                      
 10 A(J, K) = A(J, K) - X*A(I1, K)                                         
!                                                                 
! ELIMINATION FINISHED, NOW BACKSUBSTITUTION                       
!                                                                 
 50 do 20 IP = 1, N                                                      
    I = N1-IP                                                           
    do 20 K = N1, NT                                                     
    A(I, K) = A(I, K)/A(I, I)                                            
    if(I .EQ. 1) GO TO 20                                            
    I1 = I-1                                                            
    do 25 J = 1, I1                                                      
 25 A(J, K) = A(J, K) - A(I, K)*A(J, I)                                   
 20 continue                                                          
    return                                                            
    end