subroutine MARQ(FUNC, N, M, X, XLAMB, FAC, EPSG, MAXF)                   
    Use CIPR
    IMPLICIT REAL*8(A-H, O-Z)                                          
    common/COUT/IOUT                                                  
    common/CMARQ/GRAD(2), XJTJ(2, 2)                                    
    !common/CIPR/IPR                                                   
    dimension X(2), Y(2), XNY(2), A(2, 2), DX(2)                           
    IEVAL = 0                                                           
    ISTOP = 0                                                           
    IER = 0                                                             
    XMAXL = XLAMB*1.D+4                                                 
    ITER = 0                                                            
    if (IOUT.NE.6.AND.IPR.EQ.1) write(IOUT, 603)                        
    if (IPR.EQ.1) write(6, 603)                                         
    call FUNC(N, M, 1, X, SRES)                                           
    IEVAL = IEVAL+1                                                     
    SSQ = SRES                                                          
    if (IPR.EQ.1.AND.IOUT.NE.6) write(IOUT, 601) ITER, SSQ               
    if (IPR.EQ.1) write(6, 601) ITER, SSQ                                
 10 continue                                                          
    if (IEVAL.NE.1) call FUNC(N, M, 1, X, SRES)                            
    GNORM = 0.D0                                                        
    do 140 I = 1, N                                                      
140 GNORM = GNORM+GRAD(I)**2                                            
    GNORM = DSQRT(GNORM)                                                
    if (GNORM.LT.EPSG) ISTOP = ISTOP+1                                   
    if (ISTOP.GT.0) goto 1000                                          
    ITER = ITER+1                                                       
 49 continue                                                          
    if (IEVAL.GT.MAXF) goto 998                                        
    do 41 I = 1, N                                                       
    do 40 J = 1, N                                                       
 40 A(I, J) = XJTJ(I, J)                                                  
 41 A(I, I) = A(I, I)+XLAMB                                               
    call CHOL(N, A)                                                    
    Y(1) = -GRAD(1)/A(1, 1)                                              
    do 81 I = 2, N                                                       
    SUM = 0.D0                                                          
    II = I-1                                                            
    do 80 J = 1, II                                                      
 80 SUM = SUM+A(I, J)*Y(J)                                               
 81 Y(I) = (-GRAD(I)-SUM)/A(I, I)                                        
    DX(N) = Y(N)/A(N, N)                                                 
    do 85 I = 2, N                                                       
    II = N-I+1                                                          
    SUM = 0.D0                                                          
    III = II+1                                                          
    do 84 J = III, N                                                     
 84 SUM = SUM+A(J, II)*DX(J)                                             
 85 DX(II) = (Y(II)-SUM)/A(II, II)                                       
    do 90 I = 1, N                                                       
 90 XNY(I) = X(I)+DX(I)                                                 
    call FUNC(N, M, 0, XNY, SRES)                                         
    IEVAL = IEVAL+1                                                     
    SSQNY = SRES                                                        
    SQ1 = 0.D0                                                          
    SQ2 = 0.D0                                                          
    do 110 I = 1, N                                                      
    Y(I) = XNY(I)*300.D0                                                
    SQ1 = SQ1+DX(I)**2                                                  
    SQ2 = SQ2-DX(I)*GRAD(I)                                             
110 continue                                                          
    CCAL = SSQ-SSQNY                                                    
    CPRE = SQ2+XLAMB*SQ1                                                
    CALPRE = CCAL/(CPRE+1.D-14)                                         
    if (IPR.EQ.1) write(6, 601) ITER, SSQNY, Y(1), Y(2), GNORM, XLAMB, CALPRE 
    if (IOUT.NE.6.AND.IPR.EQ.1) write(IOUT, 601) ITER, SSQNY, Y(1), Y(2), GNORM, XLAMB, CALPRE                                                
    if (SSQ-SSQNY.GT..75*CPRE) XLAMB = XLAMB/FAC                         
    if (SSQ-SSQNY.LT..25*CPRE) XLAMB = XLAMB*FAC                         
    if (XLAMB.GT.XMAXL) goto 999                                       
    if (SSQNY-SSQ) 120, 120, 119                                         
119 continue                                                          
    if (SSQNY-SSQ.GT.DABS(SSQ)*.5D0) XLAMB = XLAMB*FAC                   
    goto 49                                                           
120 continue                                                          
    if (SSQ-SSQNY.GT.DABS(SSQ)*.8D0) XLAMB = XLAMB/FAC                   
    do 130 I = 1, N                                                      
130 X(I) = XNY(I)                                                       
    SSQ = SSQNY                                                         
    goto 10                                                           
998 IER = 1                                                             
    goto 1000                                                         
999 IER = 2                                                             
    goto 1000                                                         
601 format(1X, I3, 2X, D12.4, 2F10.3, 5X, 2D12.4, 3X, F10.2)                  
602 format(//, ' ISTOP = ', I2, 5X, 'IER = ', I2, 5X, 'IEVAL = ', I5, //)       
603 format(///, '  ** ITERATIONCOURSE, UNIQUAC-PARAMETERS FROM UNIFAC **'/)                                                              
1000 continue                                                          
    if (IOUT.NE.6.AND.IPR.EQ.1) write(IOUT, 602) ISTOP, IER, IEVAL        
    if (IPR.EQ.1) write(6, 602)ISTOP, IER, IEVAL                        
    return                                                            
    end