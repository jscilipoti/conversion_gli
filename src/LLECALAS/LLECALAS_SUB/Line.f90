subroutine LINE(ND, N, XLAM, GD, G, W)                                 
    IMPLICIT REAL*8(A-H, O-Z)                                          
    dimension GD(ND), G(ND, ND), W(ND, 5)                                 
 65 W(1, 2) = -GD(1)/(G(1, 1)+XLAM)                                       
    if(N.EQ.1) GO TO 75                                              
    do 70 I = 2, N                                                       
    I1 = I-1                                                            
    S = -GD(I)                                                          
    do 80 J = 1, I1                                                      
 80 S = S-G(I, J)*W(J, 2)                                                 
 70 W(I, 2) = S/(G(I, I)+XLAM)                                            
 75 W(N, 3) = W(N, 2)/(G(N, N)+XLAM)                                       
    if(N.EQ.1) GO TO 85                                              
    do 90 II = 2, N                                                      
    I = N+1-II                                                          
    I1 = I+1                                                            
    S = W(I, 2)                                                          
    do 100  J = I1, N                                                    
100 S = S-G(J, I)*W(J, 3)                                                 
 90 W(I, 3) = S/(G(I, I)+XLAM)                                            
 85 return                                                            
    end