subroutine FUNC(N, M, NDIF, X, SSQ)                                   
    IMPLICIT REAL*8(A-H, O-Z)                                          
    common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
    common/CMARQ/GRAD(2), XJTJ(2, 2)                                    
    common/CUFAC/NK, NG, P(10, 10), T                                     
    common/CA/XC(5), GE(5, 2), GC(5, 2)                                   
    dimension X(2), F(2)                                               
    JD = 0                                                              
    if (NDIF.EQ.1) JD = 4                                                
    P(1, 2) = X(1)*300.D0                                                
    P(2, 1) = X(2)*300.D0                                                
    call PARAM2                                                       
    SSQ = 0.                                                            
    if (NDIF.EQ.0) goto 11                                             
    do 10 I = 1, 2                                                       
    GRAD(I) = 0.                                                        
    do 10 J = 1, 2                                                       
 10 XJTJ(I, J) = 0.                                                      
 11 continue                                                          
    do 21 L = 1, 5                                                       
    Y1(1) = XC(L)                                                       
    Y1(2) = 1.D0-XC(L)                                                  
    call unifac(JD, Y1, ACT1, DACT1, PACT)                                
    do 17 I = 1, 2                                                       
    F(I) = ACT1(I)-GE(L, I)                                              
    GC(L, I) = ACT1(I)                                                   
 17 SSQ = SSQ+F(I)*F(I)                                                 
    if (JD.EQ.0) goto 21                                               
    do 19 I = 1, 2                                                       
    do 20 J = 1, 2                                                       
    GRAD(J) = GRAD(J)+F(I)*PACT(I, J)                                    
    do 20 K = 1, 2                                                       
 20 XJTJ(J, K) = XJTJ(J, K)+PACT(I, J)*PACT(I, K)                           
 19 continue                                                          
 21 continue                                                          
    return                                                            
    end