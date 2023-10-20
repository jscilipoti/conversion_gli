subroutine STIG(Y, S)                                              
    use CUFAC
    use CVAP
    use CGIBBS
   IMPLICIT REAL*8(A-H, O-Z)                                          
    !common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
    !common/CUFAC/NKK, NGG, Pxx(10, 10), Txx                                      
    common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
    !common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
    dimension Y(10), V(10), YGEM(10)                                    
    common/nga/nga, mass(12)
    JPUR = 0                                                            
    AMAX = 0.                                                           
    do 10 I = 1, NKK                                                       
    if (A(I).LT.AMAX) goto 10                                          
    JPUR = I                                                            
    AMAX = A(I)                                                         
 10 continue                                                          
    RMAX = 1.D5                                                         
    NGN = NKK                                                             
    if (NF.GT.1) NGN = NKK+NF                                              
    NEG = 0                                                             
    do 100 KK = 1, NGN                                                   
    JM = KK                                                             
    if (JPUR.NE.0) JM = JPUR                                             
    if (JM.LE.NKK) goto 30                                               
    do 20 I = 1, NKK                                                       
 20 Y(I) = Z(I)*(2+XVL(I, JM-NKK)/SFAS(JM-NKK))/3                            
    goto 40                                                           
 30 SUM = 0.                                                            
    do 35 I = 1, NKK                                                       
       GG = A(I)-GAM(JM, I)                                                 
       if (GG.LT.-50.D0) GG = -50.D0                                        
       Y(I) = DEXP(GG)                                                     
 35    SUM = SUM+Y(I)                                                      
 40 NA = 3                                                              
    do 43 K = 1, NA                                                      
    do 36 I = 1, NKK                                                       
 36 Y(I) = Y(I)/SUM                                                     
    call unifac(1, Y, AL, DA, PACT)                                       
    if (K.EQ.NA) goto 44                                               
    do 41 I = 1, NKK                                                       
 41 Y(I) = DEXP(A(I)-AL(I))                                             
 42 SUM = 0.                                                            
    do 43 I = 1, NKK                                                       
 43 SUM = SUM+Y(I)                                                      
 44 continue                                                          
    YV1 = 0.                                                            
    do 50 J = 1, NF                                                      
 50 V(J) = 0.                                                           
    do 60 I = 1, NKK                                                       
    GD = DLOG(Y(I))+AL(I)-A(I)                                          
    YV1 = YV1+Y(I)*GD                                                   
    do 60 J = 1, NF                                                      
    K = J                                                               
    VV = XVL(I, K)*Z(I)/SFAS(K)                                          
    D = GD*(Y(I)-VV)                                                    
 60 V(J) = V(J)+D                                                       
    YV2 = V(1)                                                          
    do 70 J = 1, NF                                                      
    if (V(J).LT.YV2) YV2 = V(J)                                          
 70 continue                                                          
    RT1 = YV1                                                           
    if (YV2.GT.0.) RT1 = RT1-YV2/2                                       
    if (NEG.EQ.0.AND.YV1.GT.0.) goto 80                                
    RT1 = YV1                                                           
    if (NEG.EQ.0) RMAX = 0.                                              
    NEG = 1                                                             
 80 if (RT1.GT.RMAX) goto 100                                          
    S = YV1                                                             
    RMAX = RT1                                                          
    CC = DEXP(-YV1)                                                     
    do 90 I = 1, NKK                                                       
 90 YGEM(I) = Y(I)*CC                                                   
    if (JPUR.NE.0) goto 110                                            
100 continue                                                          
110 do 120 I = 1, NKK                                                      
120 Y(I) = YGEM(I)                                                      
    return                                                            
    end