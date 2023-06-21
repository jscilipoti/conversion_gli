   subroutine llecalas(Tf, Pf, Zf) 
! Temperatura, presion, composicion, nro de componentes
! *******************************************************************  
! *                                                                 *  
! *  PROGRAM   L L E C A L A S                                      *  
! *  (asociacion incorporada para el calculo de flash y             *  
! *  curva binodal (icalc 0 y 1)                                    *
! *                                                                 *
! *                       ASOCIACION CRUZADA					           *		  
! *                      VERSION GENERALIZADA                       *        
! *                          Junio 2023                             *
! *  MODIFICADA POR:                                                *
! *                    ALFONSINA  ESTER ANDREATTA                   *
! *                    JOSE ANTONIO SCILIPOTI                       *  
! *                    JUAN PABLO ROVEZZI                           *
! *       Revisada en Octubre del 2007 en el chequeo de estabilidad *
! *                                                                 *  
! *        BASADA EN LAS SIMPLIFICACIONES DE LOS PAPERS:            * 	
! *  Michelsen, et al. (Fluid Phase Equilibria, 180(2001)165-174 )  *		
! *  Tan, et al.  (Ind. Eng. Chem. Res, 2004, 43, 203-208).			  *	   	
! *																                 *	   	
! *        Esto permitio  que todos los casos particulares          *         
! *       de asociacion se puedan simplificar a un unico calculo.   * 
! *                                                                 *
! *  Valido para un maximo numero grupo asociativo de 12.           *
! *  Con la implementacion en el calculo de la fraccion no asociada *   
! *  en el componente puro por  el metodo iterativo aqui            *
! *  implementado se permite queuna molecula tenga mas de un grupo  *
! *  asociativo 14/07/06.                                           * 
! *  El calculo se limita a que el numero maximo de sitios sea      * 
! *  dos (por razones matematicas).                                 *
! *                                                                 *                                                   
! *******************************************************************  
! *                                           DATE: 24/3 - 1982 /TJ *
! *******************************************************************
   
      use InputData
      use flashout
      IMPLICIT REAL*8(A-H, O-Z)      !C QUITAR                                    
      EXTERNAL STABIL, GMIX, FUNC                                         
      common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
      common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4),&
      & GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
      common/CUFAC/N, NG, P(10, 10), T                                      
      common/CY/Y13, Y21, STEP                                            
      common/CA/XC(5), GE(5, 2), GC(5, 2)                                   
      common/CIPR/IPR                                                   
      common/CQT/QT(10, 10), Q(10), R(10)                                  
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10),&
      & DACT2(10, 10), PACT(2, 2)                                                     
      common/CMODEL/MODEL                                               
      common/COUT/IOUT                                                  
      common/nga/nga, mass(12)
      common/ig/ig
      dimension DLX(10), YVAL(30), Y(10), GRAD(30), XMAT(30, 30), WORK(30, 5)  
      dimension X(2)         
      integer::MODEL, IPR, IOUT, NOVAP, ig            
      character(len = 36):: name1 
      integer:: parameters 
      real*8::Zf(10)
      dimension xmj(10), actgam(10), de(10, 10), pe(2, 2)
   
      ! Variables leidas desde main se pasan a los commons de llecalas
      model = modelo
      IPR = iprm
      IOUT = ioutm
      NOVAP = novapm
      IG = igm
      !  icalc:   0-' **** FLASH CALCULATION ****'                            
      !           1-' **** BINODAL CURVE CALCULATION ****'
      !           2-' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** '
      !  model:   0-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'     
      !           1-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'
      !  ipr:     1-' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND 
      !              UNIQUAC, RESPECTIVELY **'
      !  iout:    1-'open 'lleasoccuzada.OUT''
      !  novap:   0-'VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS'
      !  ig:      0-'write compositions'
      !           1-'write compositions and activities'
      !  ipareq:  1-'liquid-liquid parameters table (UNIFAC)'
      !           2-'liquid-vapor parameters table (UNIFAC)'
      !           3-'infinite dilution parameters table (UNIFAC)'
      !           4-'GC-EOS parameters'
  
   
    !  call ab_ban1(model) ahora se llama desde el main
	   if(IOUT.EQ.1) OPEN (UNIT = 1, FILE = 'lleasoccuzada.OUT',&
      & FORM = 'FORMATTED')
      OPEN (UNIT = 3, FILE = 'output.OUT', FORM = 'FORMATTED')
	   output = 3                                                  
      if (IOUT.EQ.0) IOUT = 6                                                                                
      if (IOUT.EQ.6) then 
         goto 5 !C Corregir
      else
         write(IOUT, 608)                                                   
         write(IOUT, 610)                                                   
         if (ICALC.EQ.0) write(IOUT, 620)                                    
         if (ICALC.EQ.1) write(IOUT, 621)                                    
         if (ICALC.EQ.2) write(IOUT, 622)                                    
         if (NOVAP.NE.0) write(IOUT, 629)                                    
         if (MODEL.EQ.0) write(IOUT, 624)                                    
         if (MODEL.EQ.1) write(IOUT, 625)                                    
         write(IOUT, 623) NTEXT  
      endif                                                                                
    5 continue                                                          
    !  CALL PARIN2     hago que se llame desde el main                                                                                                   
      T1 = 0.                                                             
      NN = 0                                                              
      T = Tf
      PP = Pf
      if (T.EQ.0.) goto 10000 
      !if (PP.EQ.0..OR.NOVAP.EQ.0) then !reemplazado por la siguiente estructura
      if (PP/=  0..and.NOVAP/=  0) then                                 
        do 1 I = 1, N                                                        
    1       PRAT(I) = DLOG(PP)-ANT(1, I)+ANT(2, I)/(T-273.15+ANT(3, I))                                                    
      endif       
      Z(:) = Zf(:)
      ZSUM = 0.                                                           
      ZMAX = 0.                                                           
      do 15 I = 1, N                                                       
        ZSUM = ZSUM+Z(I)                                                    
        if (Z(I).LT.ZMAX) cycle !goto 15                                          
        ZMAX = Z(I)                                                         
        MAXZ = I                                                            
   15 continue                                                          
      if (T.EQ.T1) goto 30                                               
      CALL PARAM2                                                       
      if (ICALC.NE.1) goto 16                                            
      if (N.NE.2.AND.N.NE.3) write(6, 616)                                
      if (IOUT.NE.6.AND.N.NE.2.AND.N.NE.3) write(IOUT, 616)               
      Y13 = Z(1)                                                          
      Y21 = Z(2)                                                          
      !write(6, 633) T                                                    
      if (IOUT.NE.6) write(IOUT, 633) T                                   
      if (N.EQ.3) goto 12                                                
      CALL SOLBIN                                                       
      goto 10000                                                        
   12 STEP = Z(3)/100.D0                                                  
      if (STEP.EQ.0.) STEP = .02D0                                         
      CALL BINOD                                                        
      goto 10000                                                        
   16 continue                                                          
      if (ICALC.NE.2) goto 19                                            
      if (N.NE.2) write(6, 616)                                           
      if (IOUT.NE.6.AND.N.NE.2) write(IOUT, 616)                          
      XC(1) = 0.                                                          
      XC(2) = .2D0                                                        
      XC(3) = .5D0                                                        
      XC(4) = .8D0                                                        
      XC(5) = 1.D0                                                        
      do 17 K = 1, 5                                                       
        Y(1) = XC(K)                                                        
        Y(2) = 1.D0-XC(K)                                                   
        CALL unifac(1, Y, ACT1, DACT1, PACT)                                  
        GE(K, 1) = ACT1(1)                                                   
   17   GE(K, 2) = ACT1(2)                                                   
      READ(2, *) R(1), Q(1)                                               
      READ(2, *) R(2), Q(2)                                                                                    
      write(6, 627)                                                      
      do 14 I = 1, 2                                                       
   14   write(6, 626) I, R(I), Q(I)                                          
      if (IOUT.EQ.6) goto 13                                             
      write(IOUT, 627)                                                   
      do 11 I = 1, 2                                                       
   11   write(IOUT, 626) I, R(I), Q(I)                                       
   13 continue                                                          
      X(1) = Z(1)/300.D0                                                  
      X(2) = Z(2)/300.D0                                                  
      do 18 I = 1, 2                                                       
        do 18 J = 1, 2                                                       
            QT(I, J) = 0.                                                        
   18       P(I, J) = 0.                                                         
      QT(1, 1) = Q(1)                                                      
      QT(2, 2) = Q(2)                                                      
      NK = 2                                                              
      NG = 2                                                              
      XLAMB = 1.                                                          
      CALL MARQ(FUNC, 2, 10, X, XLAMB, 3.D0, 1.D-7, 99)                        
      write(6, 633) T                                                    
      if (IOUT.NE.6) write(IOUT, 633) T                                   
      write(6, 617) P(1, 2), P(2, 1)       
      !(///, ' ** UNIQUAC PARAMETERS FROM UNIFAC **',
      !//, 5X, 'A12/R = ', F12.3, ' K , A21/R = ', F12.3, ' K', ///)                                  
      if (IPR.EQ.1) write(6, 618)                                         
      do 21 L = 1, 5                                                       
        do 21 I = 1, 2                                                       
            GE(L, I) = DEXP(GE(L, I))                                             
   21       GC(L, I) = DEXP(GC(L, I))                                             
      if (IPR.EQ.1) write(6, 619) ((GE(L, I), L = 1, 5), I = 1, 2)                 
      if (IPR.EQ.1) write(6, 619) ((GC(L, I), L = 1, 5), I = 1, 2)                 
      if (IOUT.EQ.6) goto 22                                             
      write(IOUT, 617) P(1, 2), P(2, 1)                                     
      if (IPR.EQ.1) write(IOUT, 618)                                      
      if (IPR.EQ.1) write(IOUT, 619) ((GE(L, I), L = 1, 5), I = 1, 2)              
      if (IPR.EQ.1) write(IOUT, 619) ((GC(L, I), L = 1, 5), I = 1, 2)              
   22 continue                                                          
      goto 10000                                                        
   19 continue                                                          
      do 20 I = 1, N                                                       
        do 20 J = 1, N                                                       
            GAM(I, J) = 0.D0                                                     
            if (J.EQ.I) goto 20                                                
            CALL GAMINF(I, J, G)                                                
            GAM(I, J) = G                                                        
   20 continue                                                          
   30 T1 = T                                                              
      NN = NN+1                                                                                                       
      do 35 I = 1, N                                                       
   35 Z(I) = Z(I)/ZSUM                                                                                 
      if (IOUT.NE.6) write(IOUT, 602) NN                                  
      if (IOUT.NE.6) write(IOUT, 605) T, PP, ZSUM, (Z(I), I = 1, N)              
      CALL unifac(1, Z, AL, DA, PACT)                                       
      SFAS(1) = 1.                                                        
      GNUL = 0.                                                           
      do 40 I = 1, N                                                       
        XVL(I, 1) = 1.                                                       
        Z(I) = Z(I)+1.D-20                                                  
        DLX(I) = DLOG(Z(I))                                                 
        A(I) = AL(I)+DLX(I)                                                 
   40   GNUL = GNUL+Z(I)*AL(I)                                              
      NF = 1                                                              
   50 CALL STIG(Y, S)                                                    
      if (S.GT.-1.D-7) goto 70                                                                                              
      if (IOUT.NE.6) write(IOUT, 603)                                     
      do 60 I = 1, N                                                       
        YVAL(I) = 1.D-5*Y(I)/Z(I)                                           
   60 continue                                                          
      goto 100                                                          
   70 do 75 I = 1, N                                                       
   75   YVAL(I) = DLOG(Y(I))                                                
      XLAM = 1.                                                                              
      if (IOUT.NE.6.AND.NF.EQ.1.AND.IPR.GT.0) write(IOUT, 606)            
      if (IOUT.NE.6.AND.NF.GT.1.AND.IPR.GT.0) write(IOUT, 609) NF         
      CALL TMSSJ(30, N, IPR, 15, XLAM, 1.D-12, FUN, YVAL, GRAD, XMAT, WORK, 1)     
      if (FUN.LT.-1.D-7) goto 80                                         
      write(output, *) 1
        write(output, 2613) (Z(j), J = 1, N)
        write(output, 2613) (AL(j), j = 1, N)        
      write(output, *) "SYSTEM IS STABLE"                                                   

	 write(7, 46) T, (xM(l, 1), l = 1, N)   
	 write(7, 46) T, (xM(l, 2), l = 1, N)
	 write(7, *)                    



      if (IOUT.NE.6) write(IOUT, 604)                                                                                         
   80 if (IOUT.NE.6) write(IOUT, 603)                                     
      do 90 I = 1, N                                                       
   90   YVAL(I) = 1.D-5*DEXP(YVAL(I))/Z(I)                                  
  100 NF = NF+1                                                           
  104 do 105 I = 1, N                                                      
        if (YVAL(I).GT.1.D0) goto 106                                      
  105 continue                                                          
      goto 109                                                          
  106 do 107 I = 1, N                                                      
  107   YVAL(I) = YVAL(I)/10.                                               
      goto 104                                                          
  109 continue                                                          
      SFAS(NF) = 1.                                                       
      XLAM = .2                                                           
      if (NF.EQ.2) XLAM = .5                                               
      M = (NF-1)*N                                                                                          
      if (IOUT.NE.6.AND.IPR.GT.0) write(IOUT, 607) NF                     
      CALL TMSSJ(30, M, IPR, 60, XLAM, 1.D-16, FUN, YVAL, GRAD, XMAT, WORK, 2)     
      NT = NF*N                                                           
      NB = NT-N                                                           
      do 110 I = 1, NB                                                     
  110   YVAL(NT+1-I) = YVAL(NB+1-I)                                                                                     
      NVAP = 0                                                            
      do 111 J = 1, NF                                                     
        if (IDUM(J).EQ.1) NVAP = J                                           
  111 continue                                                                                         
      if (IOUT.NE.6.AND.NVAP.EQ.0) write(IOUT, 630)                       
      if (IOUT.NE.6.AND.NVAP.NE.0) write(IOUT, 631) NVAP                                                   
      if (IOUT.NE.6) write(IOUT, 614) NF                                  
      if (IOUT.NE.6) write(IOUT, 611)(J, SFAS(J), J = 1, NF)                   
      if (IOUT.NE.6) write(IOUT, 612) (J, J = 1, NF)                          
      SUM = 0.                                                            
      do 115 I = 1, N                                                      
        DLX(I) = XVL(I, NF)*Z(I)/SFAS(NF)                                    
  115   SUM = SUM+DLX(I)                                                    
      SUM = DLOG(SUM)                                                     
      CALL unifac(1, DLX, A, DA, PACT)                                      
      do 120 I = 1, N                                                      
        DLX(I) = DLOG(DLX(I))                                               
  120   A(I) = A(I)+DLX(I)-SUM                                              
!-------------------------------------------------------------------------------
      do 1130 j = 1, nf
        do 1131 i = 1, n
 1131       xmj(i) = xm(i, j)
        call unifac(1, xmj, actgam, de, pe)
        do 1132 i = 1, n
 1132       agam(i, j) = actgam(i)
 1130 continue
      write (output, *) NF
      compfases(:, :) = xm(:, :)
      do i = 1, NF !escribe resultados para el output que ser� le�do por excel
        write(output, 2613) (XM(j, i), J = 1, N)
        write(output, 2613) (agam(j, i), j = 1, N)  
      enddo
    
      if (IOUT.EQ.6) goto 132                                            
      do 131 I = 1, N                                                      
      write(IOUT, 613) I, (XM(I, J), J = 1, NF)    
  131 write(iout, 1613) i, (agam(i, j), j = 1, nf)
   46	format (2X, F12.2, 8X, F12.6, 8X, F12.6 , 8X, F12.6, 8X, F12.6, &                                                 
    8X, F12.6, 8X, F12.6, 8X, F12.6, 8X, F12.6, 8X, F12.6, 8X, F12.6, 8X, F12.6)
!-------------------------------------------------------------------------------
  132 continue                                                          
                                                       


  501 format(36A2)                                                      
  502 format(8F10.2)                                                    
  503 format(20I3)                                                      
  602 format(///, ' * * * FLASH NUMBER', I3, ' * * *', //)                  
  603 format(/, ' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')            
  604 format(/, ' * SYSTEM IS STABLE *', /)                               
  605 format(' TEMPERATURE = ', F10.4, ' K, PRESSURE = ', F7.3, ' ATM, FEED = '&
     , F10.2, ' MOLES', /, ' FEED COMPOSITION (MOLE PERCENT):', /, 1X, 15(2PF7&
     .3))                                                              
  606 format(//, ' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')       
  607 format(/, ' PHASE SPLIT CALCULATION, ', I2, ' PHASES:')               
  608 format(1H1)                                                       
  609 format(//, ' DIFFERENTIAL STABILITY TEST FOR', I2, '-PHASE SYSTEM')  
  610 format(///)                                                       
  611 format(/, '  PHASE FRACTIONS (PERCENT):', 4(5X, I3, 2PF7.3, 5X))       
  612 format(/, '  COMPOSITION  ', 10X, 4(8X, I3, 9X))                       
  613 format('   X(', I2, ')            ', 5(8X, F12.8))                    
!-------------------------------------------------------------------------------
 1613 format('  ln(G', i2, ')            ', 5(8x, f12.8))
 2613 format(5(2x, f12.8))

!-------------------------------------------------------------------------------
  614 format(//, ' RESULT OF', I2, '-PHASE CALCULATION:')                  
  616 format(//, ' * WRONG INPUT SPECIFICATION *', //)                    
  617 format(///, ' ** UNIQUAC PARAMETERS FROM UNIFAC **', //, 5X, 'A12/R = ',&
      F12.3, ' K , A21/R = ', F12.3, ' K', ///)                         
  618 format(//,' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC,&
       RESPECTIVELY **'//)                                       
  619 format(10F12.5)                                                   
  620 format(' **** FLASH CALCULATION ****')                            
  621 format(' **** BINODAL CURVE CALCULATION ****', //)                 
  622 format(' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** ', //)                                                              
  623 format(1X, 'COMPONENTS : ', 40A2, //)                                
  624 format(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'//)     
  625 format(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'//)    
  626 format(I5, 2F15.4)                                                 
  627 format(//, ' SPECIFIED UNIQUAC R AND Q', /)                         
  628 format(/, ' IOUT = ', I2, /' IfIOUT = 0: OUTPUT ONLY ON UNIT 6', /,&
      ' IfIOUT = 1: OUTPUT ON BOTH UNIT 6 AND 1')                      
  629 format(/, ' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS', //)        
  630 format(' NO VAPOR PHASE')                                         
  631 format(' PHASE', I2, ' IS A VAPOR PHASE')                           
  633 format(///, '   TEMPERATURE = ', F10.2, ' DEG K')                     
10000 continue
      if(IOUT.EQ.1) CLOSE (UNIT = 1)
      close (unit = 3)

                                                                    
      endsubroutine llecalas   
                                                                  
      subroutine unifac(NDIF, X, ACT, DACT, PACT)                           
      IMPLICIT REAL*8(A-H, O-Z)                                          

      common/asoc/nktt, igamt(20, 12), nytt(20, 12)   
      common/nga/nga, mass(12)
      common/grupas1/rkass(6, 12, 6, 12), enass(6, 12, 6, 12),&
      deloh(6, 12, 6, 12)!Alfonsin

      common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
      common/CUFAC/NK, NG, P(10, 10), T                                     
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
      do 10 I = 1, NK                                                      
      THETA(I) = X(I)*Q(I)                                                
      PHI(I) = R(I)*X(I)                                                  
      THETS = THETS+THETA(I)                                              
   10 PHS = PHS+PHI(I)                                                    
      do 20 I = 1, NK                                                      
      RI(I) = R(I)/PHS                                                    
      RIL(I) = DLOG(RI(I))                                                
      QI(I) = Q(I)/THETS                                                  
   20 QIL(I) = DLOG(QI(I))                                                

      do 33 i = 1, nk
      goh(i) = 0.
      tgt(i) = 0.0
      xnohi0 = 0.0
      xgam = 0.0
     
      do j = 1, nk
		tgt(j) = 0.0
		xgam = 0.0
	end do
   33 continue
      do i = 1, nk
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
	    do i = 1, nk
		xnoh1(ja) = xnoh1(ja)+rngoh(i, ja)*x(i)
	  end do
      end do
	
      do ja = 1, nga
		xnoh(ja) = xnoh1(ja)/xgam
	end do
      do 40 I = 1, NG                                                      
      ETA(I) = 0.                                                         
      do 45 J = 1, NK                                                      
   45 ETA(I) = ETA(I)+S(I, J)*X(J)                                         
   40 ETAL(I) = DLOG(ETA(I))                                              
      do 55 I = 1, NG                                                      
      TETAR(I) = 0.                                                       
      do 55 J = 1, NK                                                      
   55 TETAR(I) = TETAR(I)+QT(I, J)*X(J)                                    
      do 60 I = 1, NK                                                      
      QID(I) = 1.-RI(I)/QI(I)                                             
      XX = F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                             
      XX = XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                              
      ACT(I) = XX                                                         
      do 661 J = 1, NG                                                     
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
      deloh(k, m, l, j) = (DEXP(ENASS(K, m, L, J)/T) - 1 )*RKASS(K, m, L, J)
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
      SUMAJ = SUMAJ + RNGOH(i, j)*(dlog(XOH(K, J)/XOHI0(I, K, J))+0.5D0*(XOHi0(i, K, J)-1))+0.5D0*R(i)*xnoh(j)*(1-xoh(k, j))
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
      do 61 I = 1, NK                                                      
   61 SS = SS+X(I)*(PRAT(I)-ACT(I))                                       
      if (SS.GT.0.) goto 69                                              
      NDUM = 1                                                            
      do 62 I = 1, NK                                                      
      ACT(I) = PRAT(I)                                                    
      do 62 J = 1, NK                                                      
   62 DACT(I, J) = 0.                                                      
      goto 100                                                          
   69 continue                                                          
      if (NDIF.EQ.4) goto 90                                             
      if (NDIF.LT.2) goto 100                                            
      do 70 I = 1, NK                                                      
      do 70 J = I, NK                                                      
      XX = Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))         
      do 75 K = 1, NG                                                      
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
      do 80 I = 1, NK                                                      
      GAM(I) = DEXP(ACT(I))                                               
   80 ACT(I) = GAM(I)*X(I)                                                
      do 85 I = 1, NK                                                      
      do 85 J = 1, NK                                                      
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
      do 92 I = 1, NK                                                      
      PACT(I, 1) = -DTAU(I, 1, 2)*TAU(1, 2)/T*300.D0                          
   92 PACT(I, 2) = -DTAU(I, 2, 1)*TAU(2, 1)/T*300.D0                          
  100 return                                                            
      end 
      subroutine PARIN2           
      use InputData                                      
      IMPLICIT REAL*8(A-H, O-Z)    
    !PARAMETER(NMS = 2) !, NCOM = 3, NSCM = 10)                                      
      common/asoc/nktt, igamt(20, 12), nytt(20, 12)  !, p4, p5
      common/grupas1/rkass(6, 12, 6, 12), enass(6, 12, 6, 12), deloh(6, 12, 6, 12)!Alfon
      common/nga/nga, mass(12) 
		common/ioh2/rngoh(12, 12)

	common/ioh2sis/rngoht(3, 2)
      common/CUFAC/NK, NG, P(10, 10), T                                     
      common/CQT/QT(10, 10), Q(10), R(10)                                  
      common/CMODEL/MODEL                                               
      common/COUT/IOUT                                                  
      dimension RT(10, 10), A(100, 100), NGM(10) !, MAINSG(57)                   
      dimension MS(10, 10, 2), NY(10, 20), JH(150), IH(20)  

    dimension NPUNTA(NMG), NGRUPA(NMG), NPUNT(NMG), NGRUP(NMG) !xnoh1(12), &
 
    INTEGER CS(NMG), TS(NMG, NMS), NPUNTMG(NMG), NMAINGR(NMG), J, K, CANIL(10) !, A, B, rngoh(NC, NSCM), GA, &
             !
        real*8::par             
      logical VL
      external mainsgfunc
                         
      REAL*4 RR(150), QQ(150)                                              
      REAL*4    A1(32), A2(32), A3(32), A4(32), A5(32), A6(32), A7(32), A8(32), &
     A9(32), A10(32), A11(32), A12(32), A13(32), A14(32), A15(32), A16(32), A17(32)&
     , A18(32), A19(32), A20(32), A21(32), A22(32), A23(32), A24(32), A25(32)&
     , A26(32), A27(32), A28(32), A29(32), A30(32), A31(32), A32(32)        
  
      DATA RR/.9011, .6744, .4469, .2195, 1.3454, 1.1167, .8886, 1.1173, .5313, .3652, &
      1.2663, 1.0396, .8121, 1., 3.2499, 3.2491, .92, .8952, 1.6724, 1.4457, &
     .998, 3.168, 1.3013, 1.528, 1.9031, 1.6764, 1.145, .9183, .6908, .9183, 1.4654, &
     1.238, 1.006, 2.2564, 2.0606, 1.8016, 2.87, 2.6401, 3.39, 1.1562, 1.870132&
     , 1.6434, 1.06, 2.0086, 1.7818, 1.5544, 1.4199, 2.4088, 4.0013, 2.9993, 2.83&
     , 2.667, 3.3092, 2.4317, 3.0856, 4.0358, 2.8266, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0/                    
      DATA QQ/.848, .54, .228, 0., 1.176, .867, .676, .988, .4, .12, .968, .66, .348&
     , 1.2, 3.128, 3.124, 1.4, .68, 1.488, 1.18, .948, 2.484, 1.224, 1.532, 1.728, 1.42&
     , 1.088, .78, .468, 1.1, 1.264, .952, .724, 1.988, 1.684, 1.448, 2.41, 2.184&
     , 2.91, .844, 1.724, 1.416, .816, 1.868, 1.56, 1.248, 1.104, 2.248, 3.568, 2.113&
     , 1.833, 1.553, 2.86, 2.192, 2.736, 3.2, 2.472, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
     0.0, 0.0, 0.0/                        
      DATA A1/0., 292.3, 156.5, 104.4, 328.2, -136.7, -131.9, 342.4, -159.8, 66.56&
     , 146.1, 14.78, 1744., -320.1, 1571., 73.8, 27.9, 21.23, 89.97, -59.06, 29.08&
     , 175.8, 94.34, 193.6, 108.5, 81.49, -128.8, 147.3, -11.91, 14.91, 67.84, 36.42/                                                              
      DATA A2/74.54, 0., -94.78, -269.7, 470.7, -135.7, -135.7, 220.6, 1., 306.1, &
     517., 1., -48.52, 485.6, 76.44, -24.36, -52.71, -185.1, -293.7, 1., 34.78, 1.&
     , 375.4, 5*1., 176.7, 132.1, 42.73, 60.82/                              
      DATA A3/-114.8, 340.7, 0., -146.8, -9.21, -223., -252., 372.8, -473.2, -78.31&
     , -75.3, -10.44, 75.49, 114.8, 52.13, 4.68, 1., 288.5, -4.7, 777.8, 56.41, -218.9&
     , 113.6, 7.18, 247.3, -50.71, -255.3, 1., -80.48, -17.78, 59.16, 29.77/
      DATA A4/-115.7, 4102., 167., 0., 1.27, -162.6, -273.6, 203.7, -470.4, -73.87&
     , 223.2, -184.9, 147.3, -170., 65.69, 122.9, 1., 33.61, 134.7, -47.13, -53.29&
     , -15.41, -97.05, -127.1, 453.4, -30.28, -124.6, 3*1., 26.59, 55.97      /
      DATA A5/644.6, 724.4, 703.9, 4000., 0., -281.1, -268.8, -122.4, -63.15, 216.&
     , -431.3, 444.7, 118.4, 180.6, 137.1, 455.1, 669.2, 418.4, 713.5, 1989., 2011.&
     , 529., 483.8, 332.6, -289.3, -99.56, -319.2, 837.9, 4*1.              /
      DATA A6/329.6, 1731., 511.5, 136.6, 937.3, 2*0., 247., -547., 401.7, 643.4, &
     -94.64, 728.7, -76.64, -218.1, 351.5, -186.1, -465.7, -260.3, 3*1., 264.7, 9*1./                                                              
      DATA A7/310.7, 1731., 577.3, 906.8, 991.3, 2*0., 104.9, -547.2, -127.6, 231.4&
     , 732.3, 349.1, -152.8, -218.1, 351.5, -401.6, -465.7, 512.2, 3*1., 264.7, 9*1./                                                             
      DATA A8/1300., 896., 859.4, 5695., 28.73, -61.29, 5.89, 0., -595.9, 634.8, 623.7&
     , 211.6, 652.3, 385.9, 212.8, 770., 740.4, 793.2, 1205., 390.7, 63.48, -239.8&
     , 13.32, 439.9, -424.3, 1., 203., 1153., -311., -262.6, 1.11, 1.       /
      DATA A9/2255., 1., 1649., 292.6, -195.5, -153.2, -153.2, 344.5, 0., -568., 3*1.&
     , -337.3, 4*1., 1616., 2*1., -860.3, 1., -230.4, 523., 1., -222.7, 5*1.  /
      DATA A10/472.6, 343.7, 593.7, 916.7, 67.07, -47.41, 353.8, -171.8, -825.7, &
     0., 128., 48.93, -101.3, 58.84, 52.38, 483.9, 550.6, 342.2, 550., 190.5, -349.2&
     , 857.7, 377., 211.6, 82.77, 2*1., 417.4, 4*1.          /              
      DATA A11/158.1, -214.7, 362.3, 1218., 1409., -344.1, -338.6, -349.9, 1., -37.36, &
      0., -311.6, 1051., 1090., 1., -47.51, 16*1./                       
      DATA A12/383., 1., 31.14, 715.6, -140.3, 299.3, -241.8, 66.95, 1., 120.3, 1724.&
     , 0., -115.7, -46.13, 2*1., 808.8, 203.1, 70.14, 5*1., -75.23, 1., -201.9, 123.2, 1., &
     -281.9, 2*1./                                             
      DATA A13/139.4, 1647., 461.8, 339.1, -104., 244.4, -57.98, -465.7, 1., 1247.&
     , .75, 1919., 0., 1417., 1402., 337.1, 437.7, 370.4, 438.1, 1349., 1., 681.4, &
     152.4, 1., -1707., 2*1., 639.7, 4*1.          /                        
      DATA A14/972.4, -577.5, 6., 5688., 195.6, 19.57, 487.1, -6.32, -898.3, 258.70&
     , -245.8, 57.7, -117.6, 0., 461.3, 1., -132.9, 176.5, 129.5, -246.3, 2.41, 3*1., &
     29.86, 7*1./                      
      DATA A15/662.1, 289.3, 32.14, 213.1, 262.5, 1970., 1970., 64.42, 1., 5.202, 2*1.&
     , -96.62, -235.7, 0., 225.4, -197.7, -20.93, 113.9, 3*1., -94.49, 9*1. /
      DATA A16/42.14, 99.61, -18.81, -114.1, 62.05, -166.4, -166.4, 315.9, 1., 1000.&
     , 751.8, 1., 19.77, 1., 301.1, 0., -21.35, -157.1, 11.8, 13*1.       /   
      DATA A17/-243.9, 337.1, 2*1., 272.2, 128.6, 507.8, 370.7, 1., -301., 1., -347.9, &
      1670., 108.9, 137.8, 110.5, 0., 1., 17.97, 13*1.         /           
      DATA A18/7.5, 4583., -231.9, -12.14, -61.57, 2*1544., 356.8, 1., 12.01, 1., &
     -249.3, 48.15, -209.7, -154.3, 249.2, 1., 0., 51.9, 1., -15.62, -216.3, 4*1., &
     -114.7, 5*1. /                                                     
      DATA A19/-5.55, 5831., 3000., -141.3, -41.75, 224.6, -207., 502.9, 4894., &
     -10.88, 1., 61.59, 43.83, 54.57, 47.67, 62.42, 56.33, -30.1, 0., -255.4, -54.86, &
     8455., -34.68, 514.6, 8*1. /                                       
      DATA A20/924.8, 1., -878.1, -107.3, -597.1, 2*1., -97.27, 1., 902.6, 2*1., &
     874.3, 629., 4*1., 475.8, 0., -465.2, 1., 794.4, 1., -241.7, 1., -906.5, 5*1. /
      DATA A21/696.8, 405.9, 29.13, 1208., -189.3, 2*1., 198.3, 1., 430.6, 3*1., &
      -149.2, 3*1., 70.04, 492., 346.2, 0., 5*1., -169.7, 5*1.        /          
      DATA A22/902.2, 1., 1.64, 689.6, -348.2, 1., 1., -109.8, -851.6, 1010., 2*1., &
      942.2, 4*1., -75.5, 1302., 2*1., 0., 1., 175.8, 164.4, 1., -944.9, 5*1./    
      DATA A23/556.7, 425.7, -1.77, 3629., -30.7, 150.8, 150.8, 1538.6, 1., 400., &
      2*1., 446.3, 1., 95.18, 3*1., 490.9, -154.5, 2*1., 0., 1., 481.3, 7*1.      /
      DATA A24/575.7, 1., -11.19, -175.6, -159., 2*1., 32.92, -16.13, -328.6, 8*1., &
      534.7, 2*1., 179.9, 1., 0., -246., 7*1.           /                   
      DATA A25/527.5, 1., 358.9, 337.7, 536.6, 2*1., -269.2, -538.6, 211.6, 1., &
      -278.2, 572.7, 343.1, 5*1., 124.8, 1., 125.3, 139.8, 963., 0., 7*1.    /      
      DATA A26/269.2, 1., 363.5, 1023., 53.37, 20*1., 0., 6*1.          /      
      DATA A27/-300., 1., -578.2, -390.7, 183.3, 2*1., -873.6, -637.3, 2*1., -208.4, &
      5*1., 18.98, 1., -387.7, 134.3, 924.5, 4*1., 0., 5*1.     /            
      DATA A28/-63.6, 3*1., -44.44, 2*1., 1429., 1., 148., 1., -13.91, -2.16, 14*1., 0., &
      4*1./                                                        
      DATA A29/928.3, 500.7, 364.2, 4*1., -364.2, 20*1., 0., 3*1.    /         
      DATA A30/331., 115.4, -58.1, 4*1., -117.4, 3*1., 173.8, 17*1., 0., 2*1.   /
      DATA A31/561.4, 784.4, 21.97, 238., 3*1., 18.41, 22*1., 0., 1.       /    
      DATA A32/956.5, 265.4, 84.16, 132.2, 27*1., 0./                        
      do 5 I = 1, 32                                                       
      A(I, 1) = A1(I)                                                      
      A(I, 2) = A2(I)                                                      
      A(I, 3) = A3(I)                                                      
      A(I, 4) = A4(I)                                                      
      A(I, 5) = A5(I)                                                      
      A(I, 6) = A6(I)                                                      
      A(I, 7) = A7(I)                                                      
      A(I, 8) = A8(I)                                                      
      A(I, 9) = A9(I)                                                      
      A(I, 10) = A10(I)                                                    
      A(I, 11) = A11(I)                                                    
      A(I, 12) = A12(I)                                                    
      A(I, 13) = A13(I)                                                    
      A(I, 14) = A14(I)                                                    
      A(I, 15) = A15(I)                                                    
      A(I, 16) = A16(I)                                                    
      A(I, 17) = A17(I)                                                    
      A(I, 18) = A18(I)                                                    
      A(I, 19) = A19(I)                                                    
      A(I, 20) = A20(I)                                                    
      A(I, 21) = A21(I)                                                    
      A(I, 22) = A22(I)                                                    
      A(I, 23) = A23(I)                                                    
      A(I, 24) = A24(I)                                                    
      A(I, 25) = A25(I)                                                    
      A(I, 26) = A26(I)                                                    
      A(I, 27) = A27(I)                                                    
      A(I, 28) = A28(I)                                                    
      A(I, 29) = A29(I)                                                    
      A(I, 30) = A30(I)                                                    
      A(I, 31) = A31(I)                                                    
    5 A(I, 32) = A32(I)                                                    
      if (IOUT.EQ.0) IOUT = 6                                              
                                                  
      READ(2, *) NK                                                      
      NC = NK                                                           
      do 15 I = 1, 10                                                      
        do 15 J = 1, NK                                                      
            QT(I, J) = 0.D0                                                      
   15 RT(I, J) = 0.D0                                                      
      if (MODEL.NE.1) goto 19                                            
      NG = NK                                                             
      do 16 I = 1, NK                                                      
                   
   16 READ(2, *) RT(I, I), QT(I, I), (P(I, J), J = 1, NK)                         
   19 continue                                                          
      if (MODEL.EQ.1) goto 21                                            
                                    

!Lectura de par�metros R, Q, int
    read(2, *) IOWNRQ, IOWNP                                            

    if(IOWNRQ/=  0)then !Lectura de R y Q
        do I = 1, IOWNRQ                                                   
            read(2, *) K, RR(K), QQ(K)  !K = n�m de grupo                                     
        enddo
    endif            

10  continue
  
    if(IOWNP/=  0)then !Lectura de par�metros de int
        do I = 1, IOWNP                                                   
            read(2, *) J, K, par !A(J, K) !j y k = n�m de grupos
            j = mainsgfunc(j, ipareq)
            k = mainsgfunc(k, ipareq)
            a(j, k) = par
        enddo
    endif
                                                  
14  continue                                                          
!-------------------------------------------------------------------------------
!ccccccccccccccccccccccccccAlfonsinaccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!rngoh = 0.
	!LECTURA DE LA COMPOSICION GRUPAL ASOCIATIVA
           do ja = 1, nga
		read(2, *) (rngoh(i, ja), i = 1, nc) 
	end do     
	


      	!LECTURA DEL NUMERO DE SITIOS Y PARAMETROS ASOCIATIVOS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC
  


	if (NGA.GT.0) READ(2, *)(MASS(I), I = 1, NGA)
    MS(:, :, :) = 0                                                       
    JH(:) = 0                                                           
    do 50 I = 1, NK !Lectura de compuestos
      write(*, *) "----"
            READ(2, *) (MS(I, J, 1), MS(I, J, 2), j = 1, size(ms(1, :, 1)))    
        do 50 j = 1, size(ms(1, :, 1))
            igamt(j, i) = ms(i, j, 1)
            nytt(j, i) = ms(i, j, 2)
            write(*, *)ms(i, j, 1)
            write(*, *)ms(i, j, 2)
50  continue	

    NUM = 0
    NMGR = 0
    NGA = 0
	do J = 1, NC
	  do I = 1, 10 
	      if (MS(J, I, 1).EQ.0)CYCLE
	      IREC2 = MS(J, I, 1)+(ipareq-1)*150      
      	    READ(14, 500, REC = IREC2)MGR, RRT, ICS, ITS1, ITS2 !500  format(2x, I2, 37X, D15.8, 120X, 3I2)
      	    if(model == 2) then
              if(ICS.NE.0) then !GRUPOS ASOCIATIVOS
                NG = MS(J, I, 1)
                CALL BUSCARAS (NG, NGRUPA, NMG, VL)
                if (VL)then
                    if (MS(J, I, 1).EQ.10) then
                        RNGOH(J, NPUNTA(NG)) = CANIL(J)
                    else
                        RNGOH(J, NPUNTA(NG)) = MS(J, I, 2)
                    endif                  
                else
                    NGA = NGA+1
                    NPUNTA(NG) = NGA
                    NGRUPA(NGA) = NG
                    CS(NGA) = ICS
                    MASS(nga) = ICS
                    TS(NGA, 1) = ITS1
                    TS(NGA, 2) = ITS2
                    if (MS(J, I, 1).EQ.10) then
                        RNGOH(J, NGA) = CANIL(J)
                    else
                        RNGOH(J, NGA) = MS(J, I, 2)
                    endif
                endif
              endif
            endif        	    
      	    if(NPUNT(MS(J, I, 1)).EQ.0) then !Ordena todos los grupos
                NUM = NUM + 1
                NPUNT(MS(J, I, 1)) = NUM
                NGRUP(NUM) = MS(J, I, 1)
               ! RR(NUM) = RRT
                CALL BUSCARAS (MGR, NMAINGR, NMG, VL)
                if (.NOT.VL) then !Crea vectores NPUNTMG y NMAINGR
                    NMGR = NMGR+1
                    NPUNTMG(MGR) = NMGR
                    NMAINGR(NMGR) = MGR
                endif
            endif 
        enddo
    enddo


!.....................................................................................
!......Genera las matrices ENASS Y RKASS con los par�metros de energ�a y volumen 
!......de asociaci�n, respectivamente, seg�n los grupos asociativos de los componentes
!......del sistema que se est� corriendo. (si el modelo elegido es A-UNIFAC)
!.....................................................................................     
    enass(:, :, :, :) = 0.0
    rkass(:, :, :, :) = 0.0
    if(model == 2)then
      do J = 1, nga
        do K = 1, nga
            do BB = 1, CS(J)
                do AA = 1, CS(K)
                    IREC1 = MAINSGfunc(NGRUPA(K), IPAREQ)+(ipareq-1)*70
                    if (TS(J, BB).EQ.1.OR.TS(K, AA).EQ.1)then !Si alguno de ambos sitios es del tipo 1
                       CALL LEEPAR (J, IREC1, IPAREQ, NGRUPA, ENASST, RKASST)
                       ENASS(AA, K, BB, J) = ENASST
                       RKASS(AA, K, BB, J) = RKASST                        
                    ELSEIF (TS(J, BB).NE.TS(K, AA)) then
                       CALL LEEPAR (J, IREC1, IPAREQ, NGRUPA, ENASST, RKASST)
                       ENASS(AA, K, BB, J) = ENASST
                       RKASS(AA, K, BB, J) = RKASST
                    else 
                       ENASS(AA, K, BB, J) = 0.0
                       RKASS(AA, K, BB, J) = 0.0     
                    endif   
                enddo
            enddo
        enddo
      enddo
    endif
      IC = 1                                                              
      do 71 I = 1, NK                                                      
        do 70 J = 1, 10                                                      
            if (MS(I, J, 1).EQ.0) goto 71                                        
            IH(IC) = MS(I, J, 1)                                                  
            if (IC.EQ.1) goto 69                                               
            if (IH(IC).EQ.IH(IC-1)) goto 70                                    
            if (IH(IC).GT.IH(IC-1)) goto 69                                    
            if (IC.GT.2) goto 55                                               
            IHH = IH(1)                                                         
            IH(1) = IH(2)                                                       
            IH(2) = IHH                                                         
            goto 69                                                           
   55       I1 = IC-1                                                           
            do 65 I2 = 1, I1                                                     
                if (IH(IC).GT.IH(I2)) goto 65                                      
                if (IH(IC).EQ.IH(I2)) goto 70                                      
                I4 = IC-I2                                                          
                do 61 I3 = 1, I4                                                     
   61               IH(IC+1-I3) = IH(IC-I3)                                             
                IH(I2) = MS(I, J, 1)                                                  
   65       continue                                                          
   69       IC = IC+1                                                           
            if (IC.GT.20) write(6, 607)                                         
            if (IOUT.NE.6.AND.IC.GT.20) write(IOUT, 607)                        
   70   continue                                                          
   71 continue                                                          
      IC = IC-1                                                           
!-------------------------------------------------------------------------------
      nktt = ic
!-------------------------------------------------------------------------------
      do 73 I = 1, IC                                                      
   73   JH(IH(I)) = I                                                       
      do 72 I = 1, 10                                                      
        do 72 J = 1, 20                                                      
   72       NY(I, J) = 0                                                         
      do 75 I = 1, NK                                                      
        do 74 J = 1, 10                                                      
            if (MS(I, J, 1).EQ.0) goto 75                                        
            N1 = MS(I, J, 1)                                                      
            N2 = MS(I, J, 2)                                                      
            if (N1.EQ.0) goto 75                                               
            N3 = JH(N1)                                                         
   74   NY(I, N3) = N2                                                       
   75 continue                                                          
      I = 0                                                               
      NGMGL = 0                                                           
      do 80 K = 1, IC                                                      
        NSG = IH(K)                                                         
        NGMNY = MAINSGfunc(NSG, ipareq)                                                 
        if (NGMNY.NE.NGMGL) I = I+1                                          
        NGM(I) = NGMNY                                                      
        NGMGL = NGMNY                                                       
      do 80 J = 1, NK                                                      
      RT(I, J) = RT(I, J)+NY(J, K)*RR(NSG)                                   
   80 QT(I, J) = QT(I, J)+NY(J, K)*QQ(NSG)                                   
      NG = I                                                                                                             
      if (IOUT.EQ.6) goto 85                                                                                             
   85 continue                                                          
      do 20 I = 1, NG                                                      
      do 20 J = 1, NG                                                      
      NI = NGM(I)                                                         
      NJ = NGM(J)                                                         
   20 P(I, J) = A(NI, NJ)                                                                                                     
      do 95 K = 1, IC                                                      
   95 NN = IH(K)                                                                                                          
      if (IOUT.EQ.6) goto 99                                             
      write(IOUT, 612)                                                   
      do 96 K = 1, IC                                                      
      NN = IH(K)                                                          
   96 write(IOUT, 613) NN, RR(NN), QQ(NN)                                  
      write(IOUT, 699)                                                   
   99 continue                                                          
   21 continue                                                                                            
      if (IOUT.EQ.6) goto 26                                             
      write(IOUT, 604)                                                   
      do 27 I = 1, NG                                                      
   27 write(IOUT, 603) (P(I, J), J = 1, NG)                                   


!ccccccccccccccccc Escritura de los par�metros de asociaci�n ALFONSINAccccccccccccccccc
       if(NGA.GT.0) then

	write(1, 218) (I, I = 1, NC)
 218	format(/, X, '"ASSOC GROUP COMPOSITION" ', /, 23X, 'COMPONENTES', /, ' GRUPO 	  #  SIT ASOC  ', 2I5, /) 

	do ja = 1, NGA
	write(1, 219) ja, MASS(ja)  , (rngoh(i, ja), i = 1, nc) 
 219	format(3X, I3, 9X, I3, 6X, 2f5.1)
	enddo

	write(1, 220)
 220	format(/, X, 'PARAMETROS DE ENERGIA DE ASOCIACION (Kelvin)  ', /)  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
	do J = 1, NGA
			if (MASS(J).EQ.0) GO TO 202
		do L = 1, NGA
				if (MASS(L).EQ.0) GO TO 104
			do I = 1, MASS(J)
				do K = 1, MASS(L)
						write(1, 221) I, J, K, L, ENASS(I, J, K, L)
  221	format(X, ' ENASS( ', I3, I3, I3, I3, ' ) = ', F10.4)
				enddo
			enddo
  104			continue
			enddo
  202			continue
	enddo


	write(1, 222)
  222	format (/, X, 'PARAMETROS DE VOLUMEN DE ASOCIACI�N (cm3/mol) ', /)

	do J = 1, NGA
			if (MASS(J).EQ.0) GO TO 301
		do L = 1, NGA
				if (MASS(L).EQ.0) GO TO 5011
			do I = 1, MASS(J)
				do K = 1, MASS(L)
					write(1, 223) I, J, K, L, RKASS(I, J, K, L)
  223	format(X, ' RKASS( ', I3, I3, I3, I3, ' ) = ', F10.4)
				enddo
			enddo
 5011			continue
			enddo
  301		   continue
	enddo
	endif 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC
!-------------------------------------------------------------------------------
      write(IOUT, 699)                                                   
      if (MODEL.EQ.0) write(IOUT, 605)                                    
      if (MODEL.EQ.1) write(IOUT, 627)                                    
   26 continue                                                          
      do 30 I = 1, NK                                                      
      Q(I) = 0.D0                                                         
      R(I) = 0.D0                                                         
      do 30 K = 1, NG                                                      
      Q(I) = Q(I)+QT(K, I)                                                 
   30 R(I) = R(I)+RT(K, I)                                                 
      !do 40 I = 1, NK                                                      
   !40 write(6, 606) I, R(I), Q(I)                                          
      if (IOUT.EQ.6) goto 42                                             
      do 41 I = 1, NK                                                      
   41 write(IOUT, 606) I, R(I), Q(I)                                       
   42 continue      
  500 format(2x, I2, 37X, D15.8, 120X, 3I2)                                                     
  501 format(20I3)                                                      
  502 format(8F10.2)                                                    
  503 format(I3, 2F10.2)                                                 
  504 format(2I3, F10.2)                                                 
  603 format(1X, 10F12.3)                                                
  604 format('  INTERACTION PARAMETERS', /)                              
  605 format(' UNIFAC MOLECULAR R AND Q', /)                             
  606 format(I5, 2F15.4)                                                 
  607 format(' ** WARNING: NUMBER OF SUB GROUPS MUST NOT EXCEED 20 **') 
  608 format(//, ' SUB GROUPS :', 20I3)                                   
  609 format(' MAIN GROUPS:', 20I3)                                      
  610 format(' COMPONENT')                                              
  611 format(6X, I2, 5X, 20I3)                                             
  612 format(' GROUP R- AND Q-VALUES', /)                                
  613 format(1X, I3, 2F10.4)                                              
  627 format(' SPECIFIED UNIQUAC R AND Q', /)                            
  699 format(//)                                                        
 1603 format('  ASSOCIATION PARAMETERS', //, 10X, 'K(OH)   :', F7.3, /, 10X, 'E(OH)/k :', F7.1, ' K-1')
      return                                                            
    endsubroutine PARIN2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
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
      subroutine STIG(Y, S)                                              
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
      common/CUFAC/N, NG, P(10, 10), T                                      
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
      common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), AA(10), DA(10, 10), XM(10, 4)                                       
      dimension Y(10), V(10), YGEM(10)                                    
      common/nga/nga, mass(12)
      JPUR = 0                                                            
      AMAX = 0.                                                           
      do 10 I = 1, N                                                       
      if (A(I).LT.AMAX) goto 10                                          
      JPUR = I                                                            
      AMAX = A(I)                                                         
   10 continue                                                          
      RMAX = 1.D5                                                         
      NGN = N                                                             
      if (NF.GT.1) NGN = N+NF                                              
      NEG = 0                                                             
      do 100 KK = 1, NGN                                                   
      JM = KK                                                             
      if (JPUR.NE.0) JM = JPUR                                             
      if (JM.LE.N) goto 30                                               
      do 20 I = 1, N                                                       
   20 Y(I) = Z(I)*(2+XVL(I, JM-N)/SFAS(JM-N))/3                            
      goto 40                                                           
   30 SUM = 0.                                                            
      do 35 I = 1, N                                                       
         GG = A(I)-GAM(JM, I)                                                 
         if (GG.LT.-50.D0) GG = -50.D0                                        
         Y(I) = DEXP(GG)                                                     
   35    SUM = SUM+Y(I)                                                      
   40 NA = 3                                                              
      do 43 K = 1, NA                                                      
      do 36 I = 1, N                                                       
   36 Y(I) = Y(I)/SUM                                                     
      CALL unifac(1, Y, AA, DA, PACT)                                       
      if (K.EQ.NA) goto 44                                               
      do 41 I = 1, N                                                       
   41 Y(I) = DEXP(A(I)-AA(I))                                             
   42 SUM = 0.                                                            
      do 43 I = 1, N                                                       
   43 SUM = SUM+Y(I)                                                      
   44 continue                                                          
      YV1 = 0.                                                            
      do 50 J = 1, NF                                                      
   50 V(J) = 0.                                                           
      do 60 I = 1, N                                                       
      GD = DLOG(Y(I))+AA(I)-A(I)                                          
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
      do 90 I = 1, N                                                       
   90 YGEM(I) = Y(I)*CC                                                   
      if (JPUR.NE.0) goto 110                                            
  100 continue                                                          
  110 do 120 I = 1, N                                                      
  120 Y(I) = YGEM(I)                                                      
      return                                                            
      end                                                              
      subroutine GMIX(NARG, NDIF, FUN, GRAD, XMAT, YVAL)                     
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
      common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), XX(10), DA(10, 10), XM(10, 4)                                       
      dimension YVAL(30), GRAD(30), X(10), TM(10), FG(10), XMAT(30, 30)       
      NG = NF-1                                                           
      N = NARG/NG                                                         
      JD = 1                                                              
      if (NDIF.EQ.2) JD = 2                                                
      if (NDIF.EQ.2) CALL CHECK(N, YVAL)                                  
      if (NF.NE.NG) goto 20                                              
      NARG = NARG-N                                                       
      NG = NG-1                                                           
      do 10 I = 1, N                                                       
      do 10 J = 1, NG                                                      
   10 YVAL(I+(J-1)*N) = XVL(I, J)                                          
   20 FUN = -GNUL                                                         
      do 50 I = 1, N                                                       
      XVL(I, NF) = 1.                                                      
      do 30 J = 1, NG                                                      
      XVL(I, J) = YVAL(I+(J-1)*N)                                          
   30 XVL(I, NF) = XVL(I, NF)-XVL(I, J)                                      
      do 40 J = 1, NF                                                      
      if (XVL(I, J).GT.0.) goto 40                                        
      FUN = 0.                                                            
      goto 1000                                                         
   40 continue                                                          
   50 continue                                                          
      do 200 J = 1, NF                                                     
      SFAS(J) = 0.                                                        
      do 60 I = 1, N                                                       
      X(I) = XVL(I, J)*Z(I)                                                
   60 SFAS(J) = SFAS(J)+X(I)                                              
      do 65 I = 1, N                                                       
      XX(I) = X(I)/SFAS(J)                                                
   65 XM(I, J) = XX(I)                                                     
      CALL unifac(JD, XX, FG, DA, PACT)                                     
      IDUM(J) = NDUM                                                      
      do 70 I = 1, N                                                       
      TM(I) = DLOG(XVL(I, J)/SFAS(J))+FG(I)                                
   70 FUN = FUN+X(I)*TM(I)                                                
      if (NDIF.EQ.0) goto 200                                            
      do 80 I = 1, N                                                       
      S = Z(I)*TM(I)                                                      
      if (J.EQ.NF) goto 75                                               
      GRAD(I+(J-1)*N) = S                                                 
      goto 80                                                           
   75 do 76 K = 1, NG                                                      
      NK = I+(K-1)*N                                                      
   76 GRAD(NK) = GRAD(NK)-S                                               
   80 continue                                                          
      if (NDIF.EQ.1) goto 200                                            
      do 100 I = 1, N                                                      
      ST = Z(I)/SFAS(J)                                                   
      do 100 L = 1, N                                                      
      S = ST*(DA(I, L)-1.)*Z(L)                                            
      if (L.EQ.I)S = S+Z(I)/XVL(I, J)                                       
      if (J.EQ.NF) goto 90                                               
      XMAT(I+(J-1)*N, L+(J-1)*N) = S                                       
      goto 95                                                           
   90 do 92 K = 1, NG                                                      
      do 92 M = 1, K                                                       
      NK = I+(K-1)*N                                                      
      NM = L+(M-1)*N                                                      
      if (K.NE.M) XMAT(NK, NM) = S                                          
      if (K.EQ.M) XMAT(NK, NM) = XMAT(NK, NM)+S                              
   92 continue                                                          
   95 continue                                                          
  100 continue                                                          
  200 continue                                                          
 1000 return                                                            
      end                                                              
      subroutine STABIL(N, NDIF, FUN, GRAD, XMAT, Y)                         
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
      common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
      dimension GRAD(30), XMAT(30, 30), Y(30), YEX(10), XEX(10)              
      common/nga/nga, mass(12)

      SUM = 0.                                                            
      do 10 I = 1, N                                                       
      YEX(I) = 0.                                                         
      if (Y(I).GT.-40.) YEX(I) = DEXP(Y(I))                                
   10 SUM = SUM+YEX(I)                                                    
      do 15 I = 1, N                                                       
   15 XEX(I) = YEX(I)/SUM                                                 
      JD = 1                                                              
      if (NDIF.EQ.2) JD = 2                                                
      CALL unifac(JD, XEX, AL, DA, PACT)                                    
      FUN = 1.                                                            
      do 20 I = 1, N                                                       
      S = Y(I)+AL(I)-A(I)                                                 
      if (NDIF.EQ.0) goto 20                                             
      GRAD(I) = YEX(I)*S                                                  
   20 FUN = FUN+YEX(I)*(S-1)                                              
      if (NDIF.LT.2) goto 50                                             
      do 30 I = 1, N                                                       
      S = XEX(I)                                                          
      do 40 J = 1, I                                                       
      XMAT(I, J) = S*YEX(J)*DA(I, J)                                        
   40 XMAT(J, I) = XMAT(I, J)                                               
   30 XMAT(I, I) = XMAT(I, I)+YEX(I)+GRAD(I)                                

  
  
   50 return                                                            
      end                                                              
      subroutine CHECK(N, YVAL)                                          
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4), GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
      common/CIPR/IPR                                                   
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
      subroutine SPLIT(ND, N, IDI, BETA, DEL, G, W)                           
      IMPLICIT REAL*8(A-H, O-Z)                                          
      dimension G(ND, ND), W(ND, 5)                                        
      IDI = 0                                                             
      do 20 I = 1, N                                                       
      I1 = I-1                                                            
      I2 = I+1                                                            
      ID = 0                                                              
      TV = G(I, I)                                                         
      SV = TV                                                             
      if(I1.EQ.0) GO TO 25                                             
      do 30 J = 1, I1                                                      
   30 SV = SV-G(I, J)**2                                                   
   25 if(SV.LT. DEL*DEL) ID = 1                                          
      SVR = DSQRT(DABS(SV))                                               
      if(SVR.LT. DEL) SVR = DEL                                          
      XM = 0.                                                             
      if(I.EQ.N) GO TO 35                                              
      do 40 J = I2, N                                                      
      S = G(J, I)                                                          
      if(I1.EQ.0) GO TO 45                                             
      do 50 K = 1, I1                                                      
   50 S = S-G(I, K)*G(J, K)                                                 
   45 S = S/SVR                                                           
      if(DABS(S).GT. XM) XM = DABS(S)                                    
   40 G(J, I) = S                                                          
   35 if(XM.LT. BETA) GO TO 55                                         
      ID = 1                                                              
      XM = XM/BETA                                                        
      do 60 J = I, N                                                       
   60 G(J, I) = G(J, I)/XM                                                  
      SVR = SVR*XM                                                        
   55 if(ID.EQ.1) W(I, 1) = SVR**2-SV                                     
      G(I, I) = SVR                                                        
      IDI = IDI+ID                                                        
   20 continue                                                          
      return                                                            
      end                                                              
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
        CALL STABIL(N, 2, F, GD, G, X)                                       
      elseif(ifunc == 2)then
        CALL GMIX(N, 2, F, GD, G, X)  
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
      CALL SPLIT(ND, N, IDI, BETA, DEL, G, W)                                 
      XLM = 0.                                                            
      NTS = 0                                                             
  350 NTS = NTS+1                                                         
      CALL LINE(ND, N, XLM, GD, G, W)                                        
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
        CALL STABIL(N, 1, FNEW, GD, G, X)                           
      elseif(ifunc == 2)then
        CALL GMIX(N, 1, FNEW, GD, G, X)   
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
   12 CALL SOLVE(Y, DY, NOLD, NEW, NITER, N, NT)                              
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
      CALL TERM(Y, DMAT, ICOND, NEW)                                       
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
      subroutine SOLVE(Y, DY, NOLD, NEW, NITER, N, NT)                        
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
      common/CBISO/NIC1, NIC2, IC1(120), IC2(120)                          
      dimension Y(4), DY(4), AMAT(3, 5)                                    
      dimension INO(3)                                                  
      common/nga/nga, mass(12)
      NK = 1000  !cambi� de 3 a 100 !Alfonsina                            
      NITER = 0                                                           
      NT = 1000  !cambi� de 3 a 100 !Alfonsina                            
      do 1 I = 1, 3                                                        
    1 if (Y(I).LT.1.D-14) NT = 1000    !Cambie de 10 iteraciones a100 !Alfonsina 0
11    NITER = NITER+1                                                     
      if (NITER.GT.NT) return                                            
      do 2 I = 1, 4                                                        
2     if (Y(I).LT.0.D0)Y(I) = 0.D0                                         
      do 3 I = 1, 2                                                        
      Y1(I) = Y(I)                                                        
3     Y2(I) = Y(I+2)                                                      
      Y1(3) = 1.D0-Y1(1)-Y1(2)                                            
      Y2(3) = 1.D0-Y2(1)-Y2(2)                                            
      if (Y1(3).LT.0.)Y1(3) = 0.                                           
      if (Y2(3).LT.0.) Y2(3) = 0.                                          
      CALL unifac(3, Y1, ACT1, DACT1, PACT)                                 
      CALL unifac(3, Y2, ACT2, DACT2, PACT)                                 
      J = 0                                                               
      do 6 I = 1, 4                                                        
      if (I.EQ.NOLD)goto 6                                               
      J = J+1                                                             
      INO(J) = I                                                          
6     continue                                                          
      do 7 I = 1, 3                                                        
      do 7 J = 1, 2                                                        
      AMAT(I, J) = DACT1(I, J)-DACT1(I, 3)                                   
7     AMAT(I, J+2) = DACT2(I, 3)-DACT2(I, J)                                 
      do 8 I = 1, 3                                                        
      AMAT(I, 5) = AMAT(I, NOLD)                                            
      do 9 J = 1, 3                                                        
9     AMAT(I, J) = AMAT(I, INO(J))                                          
8     AMAT(I, 4) = ACT1(I)-ACT2(I)                                         
      CALL GAUSL(3, 5, 3, 2, AMAT)                                          
      RES = 0.D0                                                          
      do 10 I = 1, 3                                                       
      Y(INO(I)) = Y(INO(I))-AMAT(I, 4)                                     
      DY(INO(I)) = -AMAT(I, 5)                                             
10    RES = RES+AMAT(I, 4)**2                                              
      if (RES.GT.1.D-10)goto 11                                          
      IZ = 0                                                              
      do 14 I = 1, 2                                                       
      if (Y1(I).LT.1.D-14) IZ = 1                                          
   14 if (Y2(I).LT.1.D-14) IZ = 1                                          
      if (IZ.EQ.1) goto 13                                               
      CALL GCON(3, Y1, ACT1, DACT1, ICVEX)                                  
      if (ICVEX.EQ.1) goto 15                                            
      NIC1 = NIC1+1                                                       
      IC1(NIC1) = N+1                                                     
   15 CALL GCON(3, Y2, ACT2, DACT2, ICVEX)                                  
      if (ICVEX.EQ.1) goto 13                                            
      NIC2 = NIC2+1                                                       
      IC2(NIC2) = N+1                                                     
13    DY(NOLD) = 1.D0                                                     
      NEW = NOLD                                                          
      DYMAX = 1.D0                                                        
      do 12 I = 1, 4                                                       
      if (DABS(DY(I)).LE.DYMAX)goto 12                                   
      NEW = I                                                             
      DYMAX = DABS(DY(I))                                                 
12    continue                                                          
      return                                                            
      end                                                              
      subroutine GCON(NK, X, ACT, DACT, ICVEX)                              
      IMPLICIT REAL*8(A-H, O-Z)                                          
      dimension X(3), DG(2), DDG(2, 2), ACT(3), DACT(10, 10)                  
      ICVEX = 1                                                           
      do 1 I = 1, NK                                                       
    1 if (ACT(I).LT.1.D-10) ACT(I) = 1.D-10                                
      do 5 I = 1, NK                                                       
      do 5 J = 1, NK                                                       
    5 DACT(I, J) = DACT(I, J)/ACT(I)                                        
      if (NK.EQ.3) goto 9                                                
      DDG(2, 2) = DACT(2, 2)-DACT(1, 2)-DACT(2, 1)+DACT(1, 1)                  
      goto 30                                                           
9     do 20 I = 2, NK                                                      
      II = I-1                                                            
      do 20 J = 2, NK                                                      
      JJ = J-1                                                            
   20 DDG(II, JJ) = DACT(I, J)-DACT(1, J)-DACT(I, 1)+DACT(1, 1)                
      if (X(1).LE.1.D-12.OR.X(2).LE.1.D-12) goto 30                      
      DET = DDG(1, 1)*DDG(2, 2)-DDG(2, 1)*DDG(2, 1)                           
      if (DET.LE.0.D0.OR.DDG(1, 1).LE.0.D0.OR.DDG(2, 2).LE.0.D0) ICVEX = -1  
      goto 100                                                          
   30 continue                                                          
      if (DDG(2, 2).LE.0.D0) ICVEX = -1                                     
  100 continue                                                          
      return                                                            
      end                                                              
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
      subroutine SOLBIN                                                 
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CACT/X1(10), X2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
      common/CY/Y13, Y21, STEP                                            
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
      CALL unifac(3, X1, ACT1, DACT1, PACT)                                 
      CALL unifac(3, X2, ACT2, DACT2, PACT)                                 
      do 20 I = 1, 2                                                       
      DMAT(I, 1) = DACT1(I, 1)-DACT1(I, 2)                                   
      DMAT(I, 2) = DACT2(I, 2)-DACT2(I, 1)                                   
   20 DMAT(I, 3) = ACT1(I)-ACT2(I)                                         
      CALL GAUSL(2, 3, 2, 1, DMAT)                                          
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
      CALL GCON(2, X1, ACT1, DACT1, ICVEX)                                  
      if (IOUT.NE.6.AND.ICVEX.EQ.-1) write(IOUT, 601)                     
      if (ICVEX.EQ.-1) write(6, 601)                                      
      CALL GCON(2, X2, ACT2, DACT2, ICVEX)                                  
      if (IOUT.NE.6.AND.ICVEX.EQ.-1) write(IOUT, 602)                     
      if (ICVEX.EQ.-1) write(6, 602)                                      
  601 format(' FALSE SOLUTION IN PHASE 1')                              
  602 format(' FALSE SOLUTION IN PHASE 2')                              
      return                                                            
      end                                                              
      subroutine FUNC(N, M, NDIF, X, SSQ)                                   
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10), DACT2(10, 10), PACT(2, 2)                                                     
      common/CMARQ/GRAD(2), XJTJ(2, 2)                                    
      common/CUFAC/NK, NG, P(10, 10), T                                     
      common/CA/XC(5), GE(5, 2), GC(5, 2)                                   
      dimension X(2), F(2)                                               
      JD = 0                                                              
      if (NDIF.EQ.1) JD = 4                                                
      P(1, 2) = X(1)*300.D0                                                
      P(2, 1) = X(2)*300.D0                                                
      CALL PARAM2                                                       
      SSQ = 0.                                                            
      if (NDIF.EQ.0) goto 11                                             
      do 10 I = 1, 2                                                       
      GRAD(I) = 0.                                                        
      do 10 J = 1, 2                                                       
   10 XJTJ(I, J) = 0.                                                      
   11 continue                                                          
      do 21 L = 1, 5                                                       
      Y1(1) = XC(L)                                                       
      Y1(2) = 1.D0-XC(L)                                                  
      CALL unifac(JD, Y1, ACT1, DACT1, PACT)                                
      do 17 I = 1, 2                                                       
      F(I) = ACT1(I)-GE(L, I)                                              
      GC(L, I) = ACT1(I)                                                   
   17 SSQ = SSQ+F(I)*F(I)                                                 
      if (JD.EQ.0) goto 21                                               
      do 19 I = 1, 2                                                       
      do 20 J = 1, 2                                                       
      GRAD(J) = GRAD(J)+F(I)*PACT(I, J)                                    
      do 20 K = 1, 2                                                       
   20 XJTJ(J, K) = XJTJ(J, K)+PACT(I, J)*PACT(I, K)                           
   19 continue                                                          
   21 continue                                                          
      return                                                            
      end                                                              
      subroutine MARQ(FUNC, N, M, X, XLAMB, FAC, EPSG, MAXF)                   
      IMPLICIT REAL*8(A-H, O-Z)                                          
      common/COUT/IOUT                                                  
      common/CMARQ/GRAD(2), XJTJ(2, 2)                                    
      common/CIPR/IPR                                                   
      dimension X(2), Y(2), XNY(2), A(2, 2), DX(2)                           
      IEVAL = 0                                                           
      ISTOP = 0                                                           
      IER = 0                                                             
      XMAXL = XLAMB*1.D+4                                                 
      ITER = 0                                                            
      if (IOUT.NE.6.AND.IPR.EQ.1) write(IOUT, 603)                        
      if (IPR.EQ.1) write(6, 603)                                         
      CALL FUNC(N, M, 1, X, SRES)                                           
      IEVAL = IEVAL+1                                                     
      SSQ = SRES                                                          
      if (IPR.EQ.1.AND.IOUT.NE.6) write(IOUT, 601) ITER, SSQ               
      if (IPR.EQ.1) write(6, 601) ITER, SSQ                                
   10 continue                                                          
      if (IEVAL.NE.1) CALL FUNC(N, M, 1, X, SRES)                            
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
      CALL CHOL(N, A)                                                    
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
      CALL FUNC(N, M, 0, XNY, SRES)                                         
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
      subroutine CHOL(N, A)                                              
      IMPLICIT REAL*8(A-H, O-Z)                                          
      dimension A(2, 2)                                                  
      do 50 I = 1, N                                                       
      I1 = I-1                                                            
      if (I1.EQ.0) goto 30                                               
      do 20 J = I, N                                                       
      do 20 K = 1, I1                                                      
   20 A(I, J) = A(I, J)-A(I, K)*A(J, K)                                       
   30 if (A(I, I).LT.1.D-14) A(I, I) = 1.D-14                                
      A(I, I) = DSQRT(A(I, I))                                              
      if (I.EQ.N) goto 100                                               
      J1 = I+1                                                            
      do 50 J = J1, N                                                      
   50 A(J, I) = A(I, J)/A(I, I)                                              
  100 return                                                            
      end                                                              
! subroutine GAUSL SOLVES N LINEAR ALGEBRAIC EQUATIONS BY GAUSS        
! ELIMINATION WITH ROW PIVOTING                                        
! TO SOLVE THE PROBLEM QX = U, WHERE Q IS A NXN MATRIX AND U IS NXNS, 
! ONE PLACES Q IN THE FIRST N COLUMNS OF A AND U IS PLACED IN THE      
! FOLLOWING NS COLUMNS.                                                
! THE PROGRAM RETURNS X = Q**(-1)*U AT THE PREVIOUS POSITION OF U.       
! *                                                                    
! ND IS THE ROW dimension AND NCOL IS THE COLUMN dimension OF A.       
! BOTH MUST BE TRANSFERRED TO THE subroutine.                          
! *****************                                                    
!                                                                 
      subroutine GAUSL(ND, NCOL, N, NS, A)                                  
                                                                     
      IMPLICIT REAL*8 (A-H, O-Z)                                         
      dimension A(ND, NCOL)                                              
      N1 = N+1                                                            
      NT = N+NS                                                           
      if(N .EQ. 1) GO TO 50                                            
! START ELIMINATION                                                
!                                                                 
!                                                                 
      do 10 I = 2, N                                                       
      IP = I-1                                                            
      I1 = IP                                                             
      X = DABS(A(I1, I1))                                                  
      do 11 J = I, N                                                       
      if(DABS(A(J, I1)) .LT. X) GO TO 11                                
      X = DABS(A(J, I1))                                                   
      IP = J                                                              
   11 continue                                                          
      if(IP .EQ. I1) GO TO 13                                          
!                                                                 
! ROW INTERCHANGE                                                   
!                                                                 
      do 12 J = I1, NT                                                     
      X = A(I1, J)                                                         
      A(I1, J) = A(IP, J)                                                   
   12 A(IP, J) = X                                                         
   13 do 10 J = I, N                                                       
      X = A(J, I1)/A(I1, I1)                                                
      do 10 K = I, NT                                                      
   10 A(J, K) = A(J, K) - X*A(I1, K)                                         
!                                                                 
! ELIMINATION FINISHED, NOW BACKSUBSTITUTION                       
!                                                                 
   50 do 20 IP = 1, N                                                      
      I = N1-IP                                                           
      do 20 K = N1, NT                                                     
      A(I, K) = A(I, K)/A(I, I)                                            
      if(I .EQ. 1) GO TO 20                                            
      I1 = I-1                                                            
      do 25 J = 1, I1                                                      
   25 A(J, K) = A(J, K) - A(I, K)*A(J, I)                                   
   20 continue                                                          
      return                                                            
      end                                           
      

      subroutine lubksb(a, n, np, indx, b)
      INTEGER n, np, indx(n)
      double precision a(np, np), b(n)
      INTEGER i, ii, j, ll
      double precision sum
      ii = 0
      do 12 i = 1, n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if (ii.ne.0)then
          do 11 j = ii, i-1
            sum = sum-a(i, j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii = i
        endif
        b(i) = sum
12    continue
      do 14 i = n, 1, -1
        sum = b(i)
        do 13 j = i+1, n
          sum = sum-a(i, j)*b(j)
13      continue
        b(i) = sum/a(i, i)
14    continue
      return
      END

      subroutine ludcmp(a, n, np, indx, d)
      INTEGER n, np, indx(n), NMAX
      double precision d, a(np, np), TINY
      PARAMETER (NMAX = 500, TINY = 1.0e-20)
      INTEGER i, imax, j, k
      double precision aamax, dum, sum, vv(NMAX)
      d = 1.
      do 12 i = 1, n
        aamax = 0.
        do 11 j = 1, n
          if (abs(a(i, j)).gt.aamax) aamax = abs(a(i, j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i) = 1./aamax
12    continue
      do 19 j = 1, n
        do 14 i = 1, j-1
          sum = a(i, j)
          do 13 k = 1, i-1
            sum = sum-a(i, k)*a(k, j)
13        continue
          a(i, j) = sum
14      continue
        aamax = 0.
        do 16 i = j, n

          sum = a(i, j)
          do 15 k = 1, j-1
            sum = sum-a(i, k)*a(k, j)
15        continue
          a(i, j) = sum
          dum = vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax = i
            aamax = dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k = 1, n
            dum = a(imax, k)
            a(imax, k) = a(j, k)
            a(j, k) = dum
17        continue
          d1 = -d1
          vv(imax) = vv(j)
        endif
        indx(j) = imax
        if(a(j, j).eq.0.)a(j, j) = TINY
        if(j.ne.n)then
          dum = 1./a(j, j)

          do 18 i = j+1, n
            a(i, j) = a(i, j)*dum
18        continue
        endif
19    continue
      return
      end     
                         
!******************************** F I N ***************************************
