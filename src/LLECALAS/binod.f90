subroutine BINOD                                                  
    IMPLICIT REAL*8 (A-H, O-Z)                                         
    common/CY/Y13, Y21, STEP                                            
    common/COUT/IOUT                                                  
    dimension YMAT(70, 6), NITE3(70)                                    
    common/CBISO/NIC1, NIC2, IC1(120), IC2(120)                          
    dimension Y(4), YOLD(4), DY(4), DYOLD(4), DMAT(4, 4)                   
    dimension K(6)                                                    
!-------------------------------------------------------------------------------
    dimension xandrea1(10), xandrea2(10), act1(10), act2(10), dact(10, 10), pact(2, 2), acttot(70, 6)
    common/ig/ig
!-------------------------------------------------------------------------------      
    DATA K/1, 3, 2, 4, 5, 6/                                               
    write(6, 601)                                                      
    if (IOUT.NE.6) write(IOUT, 601)                                     
601 format(///, 10X, ' ** BINODAL CURVE, CONCENTRATIONS IN MOLE PERCENT **', //, 18X, 'COMP. 1', 13X, 'COMP. 2', 13X, 'COMP. 3')                 
!-------------------------------------------------------------------------------
    write(6, 1601)
    if(iout.ne.6) write(iout, 1601)
1601 format(18x, 'gamma1', 14x, 'gamma2', 14x, 'gamma3', /)
!-------------------------------------------------------------------------------     
    do 1 I = 1, 4                                                        
    do 1 J = 1, 4                                                        
1     DMAT(I, J) = 0.D0                                                    
    IRUND = 200                                                         
    NIC1 = 0                                                            
    NIC2 = 0                                                            
    ICOND = 0                                                           
    NOLD = 2                                                            
    N = 0                                                               
    Y(1) = 1.D0-Y13/100.D0                                              
    Y(2) = 0.D0                                                         
    Y(3) = Y21/100.D0                                                   
    Y(4) = 0.D0                                                         
 12 call SOLVE(Y, DY, NOLD, NEW, NITER, N, NT)                              
    if (NITER.LE.NT) goto 16                                           
    if (N.GT.0) goto 19                                                
    ICOND = -10                                                         
    write(6, 902)                                                      
    if (IOUT.NE.6) write(IOUT, 902)                                     
902 format(/, ' THE BASE LINE CALCULATION DID NOT CONVERGE IN 10 ITERATIONS'/)                                                           
    goto 3                                                            
 19 if (IHALF.LT.5) goto 20                                            
    ICOND = -3                                                          
    goto 3                                                            
 20 IHALF = IHALF+1                                                     
    ST = ST/2.D0                                                        
    goto 17                                                           
 16 if (DABS(Y(1)-Y(3))+DABS(Y(2)-Y(4)).GT.1.D-8) goto 21              
    if (N.GT.0) goto 19                                                
    write(6, 903)                                                      
    if (IOUT.NE.6) write(IOUT, 903)                                     
903 format(/, ' THE CONCENTRATIONS ON THE BASE LINE ARE IDENTICAL IN THE TWO PHASES'/)                                                   
    goto 3                                                            
 21 N = N+1                                                             
    NITE3(N) = NITER                                                    
    IHALF = 0                                                           
    do 2 I = 1, 4                                                        
2     YMAT(N, I) = Y(I)                                                    
    if (ICOND.EQ.2.AND.Y(1).LT.1.D-10) goto 3                          
    DYMAX = DABS(DY(NEW))                                               
    do 4 I = 1, 4                                                        
4     DY(I) = DY(I)/DYMAX                                                 
    if (N.EQ.1)goto 5                                                  
    STAP = DABS(Y(NEW)-YOLD(NEW))                                       
    if (DY(NEW)*DYOLD(NEW).GT.0.D0)goto 6                              
    do 7 I = 1, 4                                                        
7     DY(I) = -DY(I)                                                      
6     if (NEW.EQ.NOLD)goto 8                                             
    RR = DY(NEW)/DYOLD(NEW)                                             
    do 9 I = 1, 4                                                        
9     DYOLD(I) = DYOLD(I)*RR                                              
8     do 10 I = 1, 4                                                       
    Z = (YOLD(I)-Y(I))/STAP                                             
    DMAT(I, 3) = (3.D0*Z+2.D0*DY(I)+DYOLD(I))/STAP                       
10    DMAT(I, 4) = (2.D0*Z+DY(I)+DYOLD(I))/STAP**2                         
5     ST = RUND(Y(NEW), DY(NEW), STEP, IRUND)                                
    do 18 I = 1, 4                                                       
    DMAT(I, 1) = Y(I)                                                    
    DMAT(I, 2) = DY(I)                                                   
    YOLD(I) = Y(I)                                                      
18    DYOLD(I) = DY(I)                                                    
17    do 11 I = 1, 4                                                       
    Y(I) = DMAT(I, 4)                                                    
    do 11 J = 1, 3                                                       
11    Y(I) = ST*Y(I)+DMAT(I, 4-J)                                          
    if (IHALF.GT.0)goto 12                                             
    call TERM(Y, DMAT, ICOND, NEW)                                       
    NOLD = NEW                                                          
    if (ICOND.EQ.0.OR.ICOND.EQ.2) goto 12                              
3     NGIT = N                                                            
    if (N.EQ.0) return                                                 
    if (ICOND.NE.1)goto 13                                             
    N = N+1                                                             
    do 14 I = 1, 4                                                       
14    YMAT(N, I) = Y(I)                                                    
    NITE3(N) = 0                                                        
13    if (N.EQ.0)return                                                  
    do 15 I = 1, N                                                       
    YMAT(I, 5) = 1.D0-YMAT(I, 1)-YMAT(I, 2)                                
15    YMAT(I, 6) = 1.D0-YMAT(I, 3)-YMAT(I, 4)                                
    do 30 I = 1, N                                                       

    xandrea1(1) = ymat(i, k(1))
    xandrea1(2) = ymat(i, k(3))
    xandrea1(3) = ymat(i, k(5))
    xandrea2(1) = ymat(i, k(2))
    xandrea2(2) = ymat(i, k(4))
    xandrea2(3) = ymat(i, k(6))
    call unifac(0, xandrea1, act1, dact, pact)
    call unifac(0, xandrea2, act2, dact, pact)
    acttot(i, 1) = dexp(act1(1))
    acttot(i, 2) = dexp(act2(1))
    acttot(i, 3) = dexp(act1(2))
    acttot(i, 4) = dexp(act2(2))
    acttot(i, 5) = dexp(act1(3))
    acttot(i, 6) = dexp(act2(3))

    write(6, 600) I, NITE3(I), (YMAT(I, K(J)), J = 1, 6)   
  write(3, 611) I, NITE3(I), (YMAT(I, K(J)), J = 1, 6)   
  
    write(8, 47) YMAT (i, 1), YMAT (i, 2), YMAT (i, 5)    !Alfonsina
   write(8, 47) YMAT (i, 3), YMAT (i, 4), YMAT (i, 6)  !Alfonsina
   write(8, *) 
   
 47	format (8X, F12.6 , 8X, F12.6 , 8X, F12.6) !Alfonsina                  

    if(ig.eq.1) write(6, 1600) (acttot(i, j), j = 1, 6)
 30 continue

600 format(5X, 2I3, 2X, 3(2P2F9.4, 2X))     
611 format(5X, 2I3, 2X, 3(2P2F9.4, 2X))                                   
    if (NIC1.GT.0) write(6, 900) IC1(1), IC1(NIC1)                       
1600 format(12x, 6(f9.3, 1x))
900 format(/, ' FALSE SOLUTION IN PHASE 1 FOR CALCULATED TIE LINES', I3, ' TO', I3, /)                                                       
    if (NIC2.GT.0) write(6, 901) IC2(1), IC2(NIC2)                       
901 format(/, ' FALSE SOLUTION IN PHASE 2 FOR CALCULATED TIE LINES', I3, ' TO', I3, /)                                                       
    if (IOUT.EQ.6) goto 50                                             
    do 35 I = 1, N                                                       
    write(IOUT, 600) I, NITE3(I), (YMAT(I, K(J)), J = 1, 6)                   
  


    if(ig.eq.1) write(iout, 1600) (acttot(i, j), j = 1, 6)
 35 continue

    if (NIC1.GT.0) write(IOUT, 900) IC1(1), IC1(NIC1)                    
    if (NIC2.GT.0) write(IOUT, 901) IC2(1), IC2(NIC2)                    
 50 continue                                                          
    return                                                            
    end