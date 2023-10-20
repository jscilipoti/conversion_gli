subroutine unifac(NDIF, X, ACT, DACT, PACT)                           
    use CUFAC
    IMPLICIT REAL*8(A-H, O-Z)                                          

    common/asoc/nktt, igamt(20, 12), nytt(20, 12)   
    common/nga/nga, mass(12)
    common/grupas1/rkass(6, 12, 6, 12), enass(6, 12, 6, 12),&
    deloh(6, 12, 6, 12)!Alfonsin

    common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
    !common/CUFAC/NKK, NGG, Pxx(10, 10), Txx                                     
    common/CPAR/TAU(10, 10), S(10, 10), F(10)                             
    common/CQT/QT(10, 10), Q(10), R(10)                                  
    dimension X(10), GAM(10), ACT(10), DACT(10, 10), THETA(10), PHI(10),&
    RI(10), QI(10), ETA(10), QIL(10), RIL(10), QID(10), ETAL(10), TETAR(10) 
    dimension U(10, 10), V(10, 10), PACT(2, 2), DTAU(2, 2, 2)                 

    dimension goh(10), xgamk(20), dxohdx(10), dxxdx(10, 10), dasdx1(10, 10),&
    dasdx2(10, 10), dasdx(10, 10)
     common/ioh2/rngoh(12, 12)

    dimension dif(12, 12), dif1(10, 12, 12) !Alfonsina
    common/zzzas/xoh(6, 12), xohi0(12, 6, 12), xoh_old(6, 12), xohi(6, 12),&
    xohi_old(6, 12), xohi0_old(12, 6, 12)  !Alfonsina
     dimension m_lambda(nga*2, nga*2), m_lambda1(nga*2, nga*2) !Alfonsina
     dimension psin(12) !Alfonsina
     dimension indx(20)
    double precision  m_lambda, m_lambda1, xoh, xohi0, xoh_old, xohi0_old  
     double precision del, del1, dif, dif1, d1, psin, xgam, xnoh1 !Alfonsina
     integer order !Alfonsina
    double precision sum1, sum2, sum3, sum4, SUMA1J, sumaj !Alfonsina
    dimension xnohi0(12, 12), tgt(12), dnohdx(12, 12), actas(12) !Alfonsina
     dimension xnoh1(12), xnoh(12), das1(3), das3(3), dxkdni(12, 6, 12),&
    dxkdnic(12, 6, 12), dgasdx(12)  !Alfonsina
    dimension dgasdxij (12, 12), drhodx(12), drhodni(12, 6, 12)
  

    actas(:) = 0.0
    dk = 1.381e-23
    deloh = 0.0
    xnoh = 0.0
    xnoh1 = 0.0
  xoh = 0.0
    xgam = 0.0
    do 7777 i = 1, 10
  xohi0 = 0
    xnohi0 = 0.0
    tgt(i) = 0.0
7777 continue

    THETS = 0.                                                          
    PHS = 0.                                                            
    do 10 I = 1, NKK                                                      
    THETA(I) = X(I)*Q(I)                                                
    PHI(I) = R(I)*X(I)                                                  
    THETS = THETS+THETA(I)                                              
 10 PHS = PHS+PHI(I)                                                    
    do 20 I = 1, NKK                                                      
    RI(I) = R(I)/PHS                                                    
    RIL(I) = DLOG(RI(I))                                                
    QI(I) = Q(I)/THETS                                                  
 20 QIL(I) = DLOG(QI(I))                                                

    do 33 i = 1, NKK
    goh(i) = 0.
    tgt(i) = 0.0
    xnohi0 = 0.0
    xgam = 0.0
   
    do j = 1, NKK
      tgt(j) = 0.0
      xgam = 0.0
  end do
 33 continue
    do i = 1, NKK
  if(nga.gt.0) then
  
      do k = 1, nktt
          tgt(i) = tgt(i)+nytt(k, i)
      end do

      do j = 1, nga  
          xnohi0(i, j) = rngoh(i, j)/R(i)  
      end do
      xgam = xgam+R(i)*x(i)
      end if  
  end do

    xnoh1 = 0d0
  do ja = 1, nga
      do i = 1, NKK
      xnoh1(ja) = xnoh1(ja)+rngoh(i, ja)*x(i)
    end do
    end do
  
    do ja = 1, nga
      xnoh(ja) = xnoh1(ja)/xgam
  end do
    do 40 I = 1, NGG                                                      
    ETA(I) = 0.                                                         
    do 45 J = 1, NKK                                                      
 45 ETA(I) = ETA(I)+S(I, J)*X(J)                                         
 40 ETAL(I) = DLOG(ETA(I))                                              
    do 55 I = 1, NGG                                                      
    TETAR(I) = 0.                                                       
    do 55 J = 1, NKK                                                      
 55 TETAR(I) = TETAR(I)+QT(I, J)*X(J)                                    
    do 60 I = 1, NKK                                                      
    QID(I) = 1.-RI(I)/QI(I)                                             
    XX = F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                             
    XX = XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                              
    ACT(I) = XX                                                         
    do 661 J = 1, NGG                                                     
    U(J, I) = S(J, I)/ETA(J)                                              
    V(J, I) = U(J, I)*TETAR(J)                                            
661 ACT(I) = ACT(I)-V(J, I)-QT(J, I)*ETAL(J)                              
!**************************calculo de las fuerzas de asociacion*****************
    if(nga.ne.0) then
      do J = 1, NGA 
            if (MASS(J).EQ.0) GO TO 201
              do m = 1, NGA
                if (MASS(m).EQ.0) GO TO 101
                      do L = 1, MASS(J)
                              do K = 1, MASS(m)
                                if (ENASS(K, m, L, J).EQ.0) then
                                       continue
                                  else
    deloh(k, m, l, j) = (DEXP(ENASS(K, m, L, J)/Txx) - 1 )*RKASS(K, m, L, J)
                                endif 
                              enddo
                      enddo

101            continue
              enddo
201          continue
      enddo


  end if  
!***********************************calculo de xoh******************************
! Calculo de la fracci�n no asociada Paper:Ind.Eng.Chem.Res, 2004, 43, 203-208   
! Inicializacion
    if(nga.ne.0) then
    xoh = 1.0d0 
    del = 1.d0
! Iteraciones con tolerancia de 10-9
    do while (del>1.0d-10)
    xoh_old = xoh
  do m = 1, nga
      do j = 1, 2
    sum1 = 0.D0
  do k = 1, nga
  sum2 = 0.D0
  do l = 1, 2

  sum2 = sum2+xoh_old(l, k)*deloh(l, k, j, m)
  end do
  sum1 = sum1+sum2*xnoh1(k)
    end do
    xoh(j, m) = 1.D0/(1.D0+sum1/xgam)          
  dif(j, m) = dabs((xoh(j, m)-xoh_old(j, m))/xoh(j, m))
  end do
  end do
  del = maxval(dif)
    end do
    end if
!cc Fin del C�lculo de la  fracci�n no asociada 
!*******calculo de xohi0(esta implementaci�n permite que una molecula tenga mas 
! de un grupo asociativo 14/07/06***********************************************
! Inicializacion
   
    if(nga.ne.0) then
    xohi0 = 1.D0
    del1 = 1.D0
! Iteraciones con tolerancia de 10-12
    do while (del1>1.0d-10)
    xohi0_old = xohi0
  do m = 1, nga
  if	(rngoh(i, m).gt.0d0) then
      do j = 1, 2
    sum3 = 0.D0
  do k = 1, nga
  sum4 = 0.D0
  do l = 1, 2 
  sum4 = sum4+ xohi0_old(i, l, k)*deloh(l, k, j, m)*xnohi0(i, k)
  end do
  sum3 = sum3+sum4
    end do
    xohi0(i, j, m) = 1.D0/(1.D0+sum3)    
   dif1(i, j, m) = dabs((xohi0(i, j, m)-xohi0_old(i, j, m))/xohi0(i, j, m))
  end do
  else
  end if
  end do
  del1 = maxval(dif1)
    end do
    end if

!******************************fin del calculo de xohi0*************************
!C�lculo del gama de asociaci�n ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! actas(M) = LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I  

     if(nga.ne.0) then	

    SUMAJ = 0.D0
    do J = 1, NGA 
    if (MASS(J).NE.0) then      
    do K = 1, MASS(J)
  if(XOH(K, J).gt.1d-13)then
    SUMAJ = SUMAJ + RNGOH(i, j)*(dlog(XOH(K, J)/XOHI0(I, K, J))&
    +0.5D0*(XOHi0(i, K, J)-1))+0.5D0*R(i)*xnoh(j)*(1-xoh(k, j))
  end if
    enddo
    else
    continue
    endif 
    enddo
    actas(I) = SUMAJ
  end if
    act(i) = act(i)+actas(i)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 60 continue
    NDUM = 0                                                            
    if (NOVAP.EQ.0) goto 69                                            
    SS = 0                                                              
    do 61 I = 1, NKK                                                      
 61 SS = SS+X(I)*(PRAT(I)-ACT(I))                                       
    if (SS.GT.0.) goto 69                                              
    NDUM = 1                                                            
    do 62 I = 1, NKK                                                      
    ACT(I) = PRAT(I)                                                    
    do 62 J = 1, NKK                                                      
 62 DACT(I, J) = 0.                                                      
    goto 100                                                          
 69 continue                                                          
    if (NDIF.EQ.4) goto 90                                             
    if (NDIF.LT.2) goto 100                                            
    do 70 I = 1, NKK                                                      
    do 70 J = I, NKK                                                      
    XX = Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))         
    do 75 K = 1, NGG                                                      
 75 XX = XX+U(K, I)*(V(K, J)-QT(K, J))-U(K, J)*QT(K, I)                      

!********************************calculo de dxkdni Alfonsina**************************************
!cCalcula los elementos de la matriz deltapq para el c�lculo de la derivada de la fracci�n 
! no asociada respecto a la fracci�n molar del componente
    psin = 0.0d0
  if(nga.ne.0) then
  m_lambda1 = 0.0d0
  m_lambda = 0.0d0
    z = 0; y = 0
  do n = 1, 2
  do m = 1, nga
       z = z+1
  do l = 1, 2
  do k = 1, nga
      y = y+1
    m_lambda(z, y) = xnoh(k)*deloh(l, k, n, m)*xoh(n, m)**2 
     if (z.eq.y)  then
    m_lambda(z, y) = m_lambda(z, y)+ 1.0d0
    end if
  end do
  end do
  y = 0
  end do
  end do
    order = nga*2
  end if 
! Calculo de  los elementos de la matriz [Yp] 
 !para el calculo de la derivada de la fraccion 
! no asociada respecto a la fracci�n molar del componente
    if(nga.ne.0) then
    do k = 1, nga

  do ll = 1, 2
  do m = 1, nga
  drhodni(j, ll, m) = (((rngoh(j, m)*xgam-xnoh1(m)*R(j)))/xgam**2)	  
  end do
  end do

    z = 0     
  do ll = 1, 2
  do m = 1, nga
  sum3 = 0.0d0
    do l = 1, nga
  sum4 = 0.0d0


  do kk = 1, 2
  sum4 = sum4+ (xoh(kk, l)*deloh(kk, l, ll, m))*drhodni(j, kk, l)
  end do
  sum3 = sum3+sum4
  end do
  z = z+1
  psin(z) = -(xoh(ll, m)**2)*sum3
  end do
  end do


    N = order
  NP = order
    m_lambda1 = m_lambda
    call  ludcmp(m_lambda1, N, NP, indx, d1)
     call lubksb(m_lambda1, N, NP, indx, psin)
! colectando las derivadas en su correspondiente sub�ndice
    z = 0
  do m = 1, 2
  do l = 1, nga
  z = z+1
  dxkdni(k, m, l) = psin(z)
  end do
  end do 
  end do


  do l = 1, nga
  do m = 1, 2
  do kk = 1, nga
  if (rngoh(i, kk).ne.0) then
    dxkdnic(j, m, l) = dxkdni(kk, m, l)   
  end if
  end do
    end do
  end do

! fin del c�lculo de la derivada de la fracci�n no asociada respecto a la 
! fracci�n molar del componente
! dgasdx(M) = derivada LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I
    if(nga.ne.0) then
 
    do l = 1, NGA 
     SUMA1J = 0.D0

    if (MASS(l).NE.0) then      
    do K = 1, MASS(l)
  if(XOH(K, l).gt.1d-13)then

    SUMA1J = SUMA1J + RNGOH(i, l)*1.D0/XOH(K, l)*dxkdnic(j, k, l)+0.5D0*&
   r(i)*(drhodni(j, k, l)-xnoh(l)*dxkdnic(j, k, l)-drhodni(j, k, l)*&
   XOH(K, l))


  end if
    enddo
    else
    continue
    endif 
  xx = xx+ SUMA1J
   end do   
  end if
      end if
    DACT(I, J) = XX                                                      
 70 DACT(J, I) = XX    
     if (NDIF.LT.3) goto 100                                           
    do 80 I = 1, NKK                                                      
    GAM(I) = DEXP(ACT(I))                                               
 80 ACT(I) = GAM(I)*X(I)                                                
    do 85 I = 1, NKK                                                      
    do 85 J = 1, NKK                                                      
    DACT(I, J) = ACT(I)*(DACT(I, J)-1.D0)                                 
    if (J.EQ.I)DACT(I, J) = DACT(I, J)+GAM(I)                              
 85 continue                                                          
    goto 100                                                          
 90 continue                                                          
    do 91 I = 1, 2                                                       
    do 91 K = 1, 2                                                       
    DTAU(I, K, K) = 0.                                                    
    do 91 L = 1, 2                                                       
    if (L.EQ.K) goto 91                                                
    H1 = TETAR(L)-QT(L, I)*ETA(L)/S(L, I)                                 
    H2 = QT(K, I)-S(L, I)*TETAR(K)/ETA(L)                                 
    DTAU(I, K, L) = -H1*H2/ETA(L)                                         
 91 continue                                                          
    do 92 I = 1, NKK                                                      
    PACT(I, 1) = -DTAU(I, 1, 2)*TAU(1, 2)/Txx*300.D0                          
 92 PACT(I, 2) = -DTAU(I, 2, 1)*TAU(2, 1)/Txx*300.D0                          
100 return                                                            
    end 