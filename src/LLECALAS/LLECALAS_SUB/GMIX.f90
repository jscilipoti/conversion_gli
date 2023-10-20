subroutine GMIX(NARG, NDIF, FUN, GRAD, XMAT, YVAL)                     
   use CVAP
   use CGIBBS 
   IMPLICIT REAL*8(A-H, O-Z)                                          
    !common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
    common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
    !common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
    dimension YVAL(30), GRAD(30), X(10), TM(10), FG(10), XMAT(30, 30)       
    NG = NF-1                                                           
    N = NARG/NG                                                         
    JD = 1                                                              
    if (NDIF.EQ.2) JD = 2                                                
    if (NDIF.EQ.2) call CHECK(N, YVAL)                                  
    if (NF.NE.NG) goto 20                                              
    NARG = NARG-N                                                       
    NG = NG-1                                                           
    do 10 I = 1, N                                                       
    do 10 J = 1, NG                                                      
 10 YVAL(I+(J-1)*N) = XVL(I, J)                                          
 20 FUN = -GNUL                                                         
    do 50 I = 1, N                                                       
    XVL(I, NF) = 1.                                                      
    do 30 J = 1, NG                                                      
    XVL(I, J) = YVAL(I+(J-1)*N)                                          
 30 XVL(I, NF) = XVL(I, NF)-XVL(I, J)                                      
    do 40 J = 1, NF                                                      
    if (XVL(I, J).GT.0.) goto 40                                        
    FUN = 0.                                                            
    goto 1000                                                         
 40 continue                                                          
 50 continue                                                          
    do 200 J = 1, NF                                                     
    SFAS(J) = 0.                                                        
    do 60 I = 1, N                                                       
    X(I) = XVL(I, J)*Z(I)                                                
 60 SFAS(J) = SFAS(J)+X(I)                                              
    do 65 I = 1, N                                                       
    AL(I) = X(I)/SFAS(J)                                                
 65 XM(I, J) = AL(I)                                                     
    call unifac(JD, AL, FG, DA, PACT)                                     
    IDUM(J) = NDUM                                                      
    do 70 I = 1, N                                                       
    TM(I) = DLOG(XVL(I, J)/SFAS(J))+FG(I)                                
 70 FUN = FUN+X(I)*TM(I)                                                
    if (NDIF.EQ.0) goto 200                                            
    do 80 I = 1, N                                                       
    S = Z(I)*TM(I)                                                      
    if (J.EQ.NF) goto 75                                               
    GRAD(I+(J-1)*N) = S                                                 
    goto 80                                                           
 75 do 76 K = 1, NG                                                      
    NK = I+(K-1)*N                                                      
 76 GRAD(NK) = GRAD(NK)-S                                               
 80 continue                                                          
    if (NDIF.EQ.1) goto 200                                            
    do 100 I = 1, N                                                      
    ST = Z(I)/SFAS(J)                                                   
    do 100 L = 1, N                                                      
    S = ST*(DA(I, L)-1.)*Z(L)                                            
    if (L.EQ.I)S = S+Z(I)/XVL(I, J)                                       
    if (J.EQ.NF) goto 90                                               
    XMAT(I+(J-1)*N, L+(J-1)*N) = S                                       
    goto 95                                                           
 90 do 92 K = 1, NG                                                      
    do 92 M = 1, K                                                       
    NK = I+(K-1)*N                                                      
    NM = L+(M-1)*N                                                      
    if (K.NE.M) XMAT(NK, NM) = S                                          
    if (K.EQ.M) XMAT(NK, NM) = XMAT(NK, NM)+S                              
 92 continue                                                          
 95 continue                                                          
100 continue                                                          
200 continue                                                          
1000 return                                                            
    end