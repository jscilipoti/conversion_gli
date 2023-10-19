subroutine PARAM2                                                 
    IMPLICIT REAL*8(A-H, O-Z)                                          
    common/CUFAC/NK, NG, P(10, 10), T                                     
    common/CPAR/TAU(10, 10), S(10, 10), F(10)                             
    common/CQT/QT(10, 10), Q(10), R(10)                                  
    do 30 I = 1, NG                                                      
      do 30 J = 1, NG                                                      
 30       TAU(I, J) = DEXP(-P(I, J)/T)                                          
 40 continue                                                          
    do 50 I = 1, NK                                                      
      do 50 K = 1, NG                                                      
          S(K, I) = 0.D0                                                       
    do 50 M = 1, NG                                                      
 50   S(K, I) = S(K, I)+QT(M, I)*TAU(M, K)                                    
    do 60 I = 1, NK                                                      
      F(I) = 1.D0                                                         
    do 60 J = 1, NG                                                      
 60   F(I) = F(I)+QT(J, I)*DLOG(S(J, I))                                    
    return                                                            
    end                                                              
    