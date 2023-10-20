subroutine PARAM2                                                 
    Use CUFAC
   IMPLICIT REAL*8(A-H, O-Z)                                          
    !common/CUFAC/NKK, NGG, Pxx(10, 10), Txx                                     
    common/CPAR/TAU(10, 10), S(10, 10), F(10)                             
    common/CQT/QT(10, 10), Q(10), R(10)                                  
    do 30 I = 1, NGG                                                      
      do 30 J = 1, NGG                                                      
 30       TAU(I, J) = DEXP(-Pxx(I, J)/Txx)                                          
 40 continue                                                          
    do 50 I = 1, NKK                                                      
      do 50 K = 1, NGG                                                      
          S(K, I) = 0.D0                                                       
    do 50 M = 1, NGG                                                      
 50   S(K, I) = S(K, I)+QT(M, I)*TAU(M, K)                                    
    do 60 I = 1, NKK                                                      
      F(I) = 1.D0                                                         
    do 60 J = 1, NGG                                                      
 60   F(I) = F(I)+QT(J, I)*DLOG(S(J, I))                                    
    return                                                            
    end                                                              
    