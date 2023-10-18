      FUNCTION RUND(Y, DY, S, IRUND)                                       
      IMPLICIT REAL*8(A-H, O-Z)                                          
      X = Y+S*DY+1.D-8*DY**2                                              
      IX = IRUND*X                                                        
      Z = DFLOAT(IX)/IRUND-Y                                              
      RUND = DABS(Z)                                                      
      return                                                            
      end                                                              
      subroutine TERM(Y, DMAT, ICOND, NEW)                                 
      IMPLICIT REAL*8(A-H, O-Z)                                          
      dimension Y(4), DMAT(4, 4), A(4)                                     
      if (Y(1).LT.1.D-14.OR.Y(3).LT.1.D-14) goto 1                       
      if (Y(1)+Y(2).GT.1.D0.OR.Y(3)+Y(4).GT.1.D0) goto 2                 
      if (Y(1)+Y(2)-.01D0.LT.Y(3)+Y(4).AND.Y(1)-.01D0.LT.Y(3))goto 3     
      return                                                            
1     ICOND = 2                                                           
      DS = DMAT(1, 1)/(DMAT(1, 1)-Y(1))                                     
      do 5 I = 1, 4                                                        
5     Y(I) = DMAT(I, 1)+DS*(Y(I)-DMAT(I, 1))                                
      Y(1) = 0.D0                                                         
      NEW = 1                                                             
      return                                                            
2     ICOND = -2                                                          
      return                                                            
3     ICOND = 1                                                           
      ND = 2+NEW                                                          
      if (ND.GT.4)ND = ND-4                                                
      do 6 I = 1, 4                                                        
6     A(I) = DMAT(NEW, I)-DMAT(ND, I)                                       
      DS = 0.D0                                                           
      NITER = 0                                                           
7     NITER = NITER+1                                                     
      if (NITER.LE.10)goto 8                                             
      ICOND = -1                                                          
      return                                                            
8     F = ((A(4)*DS+A(3))*DS+A(2))*DS+A(1)                                
      DF = (3.D0*A(4)*DS+2.D0*A(3))*DS+A(2)                               
      DF = -F/DF                                                          
      DS = DS+DF                                                          
      if (DABS(DF).GT.1.D-6)goto 7                                       
      do 9 I = 1, 4                                                        
      Y(I) = DMAT(I, 4)                                                    
      do 9 J = 1, 3                                                        
9     Y(I) = Y(I)*DS+DMAT(I, 4-J)                                          
      return                                                            
      end