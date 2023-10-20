subroutine SOLBIN                                                 
   use CY 
   IMPLICIT REAL*8(A-H, O-Z)                                          
    common/CACT/X1(10), X2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
    !common/CY/Y13, Y21, STEP                                            
    common/COUT/IOUT                                                  
    dimension DMAT(2, 3)                                               
    common/nga/nga, mass(12)

    NITER = 0                                                           
    X1(1) = 1.D0-Y13/100.D0                                             
    X2(1) = Y21/100.D0                                                  
 10 NITER = NITER+1                                                     
    if (NITER.GT.10) goto 50                                           
    if (X1(1).LT.0.D0) X1(1) = 0.D0                                      
    if (X2(1).LT.0.D0) X2(1) = 0.D0                                      
    X1(2) = 1.D0-X1(1)                                                  
    X2(2) = 1.D0-X2(1)                                                  
    call unifac(3, X1, ACT1, DACT1, PACT)                                 
    call unifac(3, X2, ACT2, DACT2, PACT)                                 
    do 20 I = 1, 2                                                       
    DMAT(I, 1) = DACT1(I, 1)-DACT1(I, 2)                                   
    DMAT(I, 2) = DACT2(I, 2)-DACT2(I, 1)                                   
 20 DMAT(I, 3) = ACT1(I)-ACT2(I)                                         
    call GAUSL(2, 3, 2, 1, DMAT)                                          
    RES = DMAT(1, 3)**2+DMAT(2, 3)**2                                     
    X1(1) = X1(1)-DMAT(1, 3)                                             
    X2(1) = X2(1)-DMAT(2, 3)                                             
    if (RES.GT.1.D-20) goto 10                                         
 50 continue                                                          
    write(6, 603)                                                      
    if (IOUT.NE.6) write(IOUT, 603)                                     
    if (IOUT.NE.6) write(IOUT, 604) X1(1), X2(1), X1(2), X2(2)             
    write(6, 604) X1(1), X2(1), X1(2), X2(2)                              
603 format(///, 5X, '** BINARY SOLUBILITIES IN MOLE FRACTIONS **', //, 11X, 'COMPONENT 1', 15X, 'COMPONENT 2', /)                               
604 format(2(2X, 2P2D12.2)//)                                          
    call GCON(2, X1, ACT1, DACT1, ICVEX)                                  
    if (IOUT.NE.6.AND.ICVEX.EQ.-1) write(IOUT, 601)                     
    if (ICVEX.EQ.-1) write(6, 601)                                      
    call GCON(2, X2, ACT2, DACT2, ICVEX)                                  
    if (IOUT.NE.6.AND.ICVEX.EQ.-1) write(IOUT, 602)                     
    if (ICVEX.EQ.-1) write(6, 602)                                      
601 format(' FALSE SOLUTION IN PHASE 1')                              
602 format(' FALSE SOLUTION IN PHASE 2')                              
    return                                                            
    end