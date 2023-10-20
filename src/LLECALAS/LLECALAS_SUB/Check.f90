subroutine CHECK(N, YVAL)
   use CGIBBS
   use CIPR                                          
    IMPLICIT REAL*8(A-H, O-Z)                                          
    !common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
    !common/CIPR/IPR                                                   
    common/COUT/IOUT                                                  
    dimension YVAL(30)                                                
    JMAX = NF                                                           
    SMAX = SFAS(NF)                                                     
    NT = NF-1                                                           
    do 5 J = 1, NT                                                       
    if (SFAS(J).LT.SMAX) goto 5                                        
    SMAX = SFAS(J)                                                      
    JMAX = J                                                            
  5 continue                                                          
    if (JMAX.EQ.NF) goto 20                                            
    do 10 I = 1, N                                                       
    DS = 1.                                                             
    do 15 J = 1, NT                                                      
 15 DS = DS-YVAL(I+(J-1)*N)                                             
 10 YVAL(I+(JMAX-1)*N) = DS                                             
 20 if (NF.EQ.2) goto 100                                              
    do 21 I = 1, N                                                       
    XVL(I, NF) = 1.                                                      
    do 21 J = 1, NT                                                      
    XVL(I, J) = YVAL(I+(J-1)*N)                                          
 21 XVL(I, NF) = XVL(I, NF)-XVL(I, J)                                      
    do 23 J = 1, NF                                                      
    SFAS(J) = 0.                                                        
    do 22 I = 1, N                                                       
    AL(I) = XVL(I, J)*Z(I)                                               
 22 SFAS(J) = SFAS(J)+AL(I)                                             
    do 23 I = 1, N                                                       
 23 XM(I, J) = AL(I)/SFAS(J)                                             
    do 30 I = 1, NT                                                      
    JN = I+1                                                            
    do 30 J = JN, NF                                                     
    if (DABS(XM(MAXZ, I)-XM(MAXZ, J)).LT.5.D-3) goto 40                  
 30 continue                                                          
    goto 100                                                          
 40 IV = I                                                              
    JV = J                                                              
    DMAX = 0.                                                           
    do 50 I = 1, N                                                       
    R1 = DABS(XM(I, IV)-XM(I, JV))                                        
    if (R1.GT.DMAX) DMAX = R1                                            
 50 continue                                                          
    if (DMAX.GT.2.5D-2) goto 100                                       
    write(6, 601)                                                      
    if (IOUT.NE.6) write(IOUT, 601)                                     
    NF = NF-1                                                           
    NT = NT-1                                                           
    do 60 I = 1, N                                                       
    XVL(I, IV) = XVL(I, IV)+XVL(I, JV)                                     
 60 XVL(I, JV) = XVL(I, NF+1)                                             
100 continue                                                          
    if (NF.LT.3) goto 250                                              
    MINF = 0                                                            
    do 200 I = 1, NF                                                     
    if (SFAS(I).LT.1.D-12) MINF = I                                      
200 continue                                                          
    if (MINF.EQ.0) goto 250                                            
    write(6, 602) MINF, SFAS(MINF)                                      
    if (IOUT.NE.6) write(IOUT, 602) MINF, SFAS(MINF)                     
    do 220 I = 1, NF                                                     
    if (I.EQ.MINF) goto 220                                            
    do 210 J = 1, N                                                      
210 XVL(J, I) = XVL(J, I)+SFAS(I)*XVL(J, MINF)                             
220 continue                                                          
    if (MINF.EQ.NF) goto 250                                           
    NNF = NF-1                                                          
    do 230 I = 1, NNF                                                    
    if (I.LT.MINF) goto 230                                            
    do 240 J = 1, N                                                      
240 XVL(J, I) = XVL(J, I+1)                                               
230 continue                                                          
    NF = NF-1                                                           
250 continue                                                          
601 format(//, ' * * * TWO PHASES ARE IDENTICAL * * *'//)              
602 format(/, ' * THE AMOUNT OF PHASE', I2, ' IS', D12.4, '. THE PHASE IS ELIMINATED *', /)                                                   
    return                                                            
    end