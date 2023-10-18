subroutine TMSSJ(ND, N, IPR, NMAX, XLAM, GLIM, F   , X, GD, G, W, ifunc)     
    IMPLICIT REAL*8(A-H, O-Z)                                          
    common/COUT/IOUT                                                  
    integer::ifunc
    dimension X(ND), GD(ND), G(ND, ND), W(ND, 5)                           
    NTM = 0                                                             
    DEL = 1.D-7                                                         
150 NTM = NTM+1                                                         
    GNM = 0.                                                            
    
    if(ifunc == 1)then
      call STABIL(N, 2, F, GD, G, X)                                       
    elseif(ifunc == 2)then
      call GMIX(N, 2, F, GD, G, X)  
    endif
    do 5 I = 1, N                                                        
  5 GNM = GNM+GD(I)**2                                                  
    if(GNM.LT. GLIM) GO TO 200                                       
    BETA = 0.                                                           
    do 10 I = 1, N                                                       
    W(I, 4) = X(I)                                                       
    W(I, 5) = GD(I)                                                      
    do 10 J = 1, I                                                       
    S = DABS(G(I, J))                                                    
    if(I.EQ.J) S = S*N                                                 
    if(S.GT. BETA) BETA = S                                            
 10 W(I, 1) = 0.                                                         
    BETA = DSQRT(BETA/N)                                                
    call SPLIT(ND, N, IDI, BETA, DEL, G, W)                                 
    XLM = 0.                                                            
    NTS = 0                                                             
350 NTS = NTS+1                                                         
    call LINE(ND, N, XLM, GD, G, W)                                        
    SMAX = 0.                                                           
    GP = 0.                                                             
    DP = 0.                                                             
    do 205 I = 1, N                                                      
    S = W(I, 3)                                                          
    if(DABS(S).GT. SMAX) SMAX = DABS(S)                                
    GP = GP+S*GD(I)                                                     
205 DP = DP+S*S*W(I, 1)                                                  
    FAC = 1.                                                            
    if(SMAX.GT. XLAM) FAC = XLAM/SMAX                                  
    DER = FAC*GP                                                        
    ALFA = 1.  
    cont = 0
210 FF = ALFA*FAC                                                       
    do 215 I = 1, N                                                      
215 X(I) = W(I, 4)+FF*W(I, 3)                                             

     if(ifunc == 1)then
      call STABIL(N, 1, FNEW, GD, G, X)                           
    elseif(ifunc == 2)then
      call GMIX(N, 1, FNEW, GD, G, X)   
    endif

    cont = cont + 1
    if (cont > 10000)then 
    goto 230
    endif

                                           
    if (FNEW.NE.0.) goto 220                                           
    ALFA = ALFA/2.                                                      
    goto 210                                                          
220 DELS = FNEW-F                                                       
     if(FNEW.NE. 0.  .AND. DELS.LT. 1.D-10) GO TO 230                
    if(NTS.GT. 1) GO TO 125                                          
    DE2 = (DELS-ALFA*DER)/ALFA**2/2                                     
    GT = -DER/DE2                                                       
    if(GT.GT. ALFA) ALFA = ALFA/3.                                     
    if(GT.LE. ALFA) ALFA = GT                                          
    GO TO 210                                                         
230 PRED = FF*GP-FF**2/2*(GP+DP)                                        
    CALPRE = DELS/PRED                                                  
    F = FNEW                                                            
    do 345 I = 1, N                                                      
345 W(I, 4) = X(I)                                                       
    if(NTS.GE.3) GO TO 125                                           
    if(ALFA.EQ. 1.  .AND. CALPRE.GT. .8) GO TO 350                   
125 do 126 I = 1, N                                                      
126 X(I) = W(I, 4)                                                       
    if(IPR.NE.0) write(6, 130) NTM, F, CALPRE, GNM, ALFA                  
    if (IOUT.NE.6.AND.IPR.NE.0) write(IOUT, 130) NTM, F, CALPRE, GNM, ALFA  
130 format(1X, I2, ':  F = ', D16.8, '  CALPRE = ', F8.3, '  GRAD = ', D16.8, '  ALFA = ', D10.2)                                                
    if(IPR.GT.1) write(6, 131) (X(I), I = 1, N)                           
    if (IOUT.NE.6.AND.IPR.GT.1) write(IOUT, 131) (X(I), I = 1, N)           
131 format(' X-VECTOR', 10F11.5, (/, 9X, 10F11.5))                        
    if(IPR.GT. 1) write(6, 132)                                       
    if (IOUT.NE.6.AND.IPR.GT.1) write(IOUT, 132)                        
132 format(' ')                                                       
    if(NTM.LT. NMAX) GO TO 150                                       
200 if (IPR.GT.0) write(6, 133) NTM, GNM, F                               
    if (IOUT.NE.6.AND.IPR.GT.0) write(IOUT, 133) NTM, GNM, F              
133 format(/, '  NUMBER OF ITERATIONS = ', I2, ', NORM OF GRADIENT = ', D12.5, ', F = ', D14.6)                                             
    if (IPR.GT.0) write(6, 131) (X(I), I = 1, N)                            
    if (IOUT.NE.6.AND.IPR.GT.0) write(IOUT, 131) (X(I), I = 1, N)           
    return                                                            
    end