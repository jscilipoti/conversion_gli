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
    subroutine GAMINF(N1, N2, GAM)                                      
    IMPLICIT REAL*8(A-H, O-Z)                                          
    common/CUFAC/NK, NG, P(10, 10), T                                     
    common/CPAR/TAU(10, 10), S(10, 10), F(10)                             
    common/CQT/QT(10, 10), Q(10), R(10)                                  
!-------------------------------------------------------------------------------
    common/asoc/nktt, igamt(20, 12), nytt(20, 12)  
    common/nga/nga, mass(12)
    common/grupas1/rkass(6, 12, 6, 12), enass(6, 12, 6, 12), deloh(6, 12, 6, 12) !Alfonsina
!-------------------------------------------------------------------------------
    dimension xohi0(10), xnohi0(10), tgt(10), xgamk(20), x(2)
  common/ioh2/rngoh(12, 12)
  common/ioh2sis/rngoht(3, 2)

!-------------------------------------------------------------------------------
    dk = 1.381e-23
    deloh = 0.0
    xnoh = 0.0
    xnoh1 = 0.0
  xoh = 0.0
    xgam = 0.0
    xgamt = 0.0
    do 7777 i = 1, 10
    xnohi0(i) = 0.0
    tgt(i) = 0.0
7777 continue
    do 8888 k = 1, 20
    xgamk(k) = 0.0
8888 continue

    Q1 = Q(N2)/Q(N1)                                                    
    R1 = R(N2)/R(N1)                                                    
    QR = R1/Q1                                                          
    GAM = F(N2)+Q(N2)*(1.-DLOG(Q1))-R1+DLOG(R1)-5.D0*Q(N2)*(1.-QR+DLOG(R1)-DLOG(Q1))                                                      
    do 10 I = 1, NG                                                      
 10 GAM = GAM-S(I, N2)/S(I, N1)*QT(I, N1)-QT(I, N2)*DLOG(S(I, N1))           

    return                                                            
    end