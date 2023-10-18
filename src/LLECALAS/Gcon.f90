subroutine GCON(NK, X, ACT, DACT, ICVEX)                              
    IMPLICIT REAL*8(A-H, O-Z)                                          
    dimension X(3), DG(2), DDG(2, 2), ACT(3), DACT(10, 10)                  
    ICVEX = 1                                                           
    do 1 I = 1, NK                                                       
  1 if (ACT(I).LT.1.D-10) ACT(I) = 1.D-10                                
    do 5 I = 1, NK                                                       
    do 5 J = 1, NK                                                       
  5 DACT(I, J) = DACT(I, J)/ACT(I)                                        
    if (NK.EQ.3) goto 9                                                
    DDG(2, 2) = DACT(2, 2)-DACT(1, 2)-DACT(2, 1)+DACT(1, 1)                  
    goto 30                                                           
9     do 20 I = 2, NK                                                      
    II = I-1                                                            
    do 20 J = 2, NK                                                      
    JJ = J-1                                                            
 20 DDG(II, JJ) = DACT(I, J)-DACT(1, J)-DACT(I, 1)+DACT(1, 1)                
    if (X(1).LE.1.D-12.OR.X(2).LE.1.D-12) goto 30                      
    DET = DDG(1, 1)*DDG(2, 2)-DDG(2, 1)*DDG(2, 1)                           
    if (DET.LE.0.D0.OR.DDG(1, 1).LE.0.D0.OR.DDG(2, 2).LE.0.D0) ICVEX = -1  
    goto 100                                                          
 30 continue                                                          
    if (DDG(2, 2).LE.0.D0) ICVEX = -1                                     
100 continue                                                          
    return                                                            
    end