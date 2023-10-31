! *****************************************************************************
! *   Subroutine: llecalas(Tf, Pf, Zf)                                        *
! *****************************************************************************
!
subroutine llecalas(Tf, Pf, Zf) 
!
! *****************************************************************************  
! *                                                                           *  
! *  PROGRAM   L L E C A L A S                                                *  
! *  (asociacion incorporada para el calculo de flash y                       *  
! *  curva binodal (icalc 0 y 1)                                              *
! *                                                                           *                          
! *                           ASOCIACION CRUZADA					               *		  
! *                          VERSION GENERALIZADA                             *        
! *                              Junio 2023                                   *
! *   Modificada por:                                                         *
! *                        ALFONSINA ESTER ANDREATTA                         *   
! *                          JOSE ANTONIO SCILIPOTI                           *  
! *                           JUAN PABLO ROVEZZI                              *
! *   Revisada en Octubre del 2007 en el chequeo de estabilidad               *
! *                                                                           *  
! *   Basada en las simplificaciones de los papers:                           * 	
! *      # Michelsen, et al. (Fluid Phase Equilibria, 180(2001)165-174 )      *		
! *      # Tan, et al.  (Ind. Eng. Chem. Res, 2004, 43, 203-208).			      *	   	
! *																                           *	   	
! *   Esto permitio  que todos los casos particulares de asociacion se puedan *         
! *   simplificar a un unico calculo.                                         * 
! *                                                                           *
! *   Valido para un maximo numero grupo asociativo de 12.                    *
! *   Con la implementacion en el calculo de la fraccion no asociada en el    *   
! *   componente puro por  el metodo iterativo aqui implementado se permite   *
! *   que una molecula tenga mas de un grupo asociativo (14/07/06).           *
! *                                                                           * 
! *   El calculo se limita a que el numero maximo de sitios sea dos           * 
! *   (por razones matematicas).                                              *
! *                                                                           *                                                   
! *****************************************************************************  
! *                        DATE: 24/3/1982 /TJ                                *
! *****************************************************************************
! 
!  Tf: Temperatura, 
!  Pf: presion, 
!  Zf: composicion, nro de componentes

      use InputData
      use flashout
      use CUFAC
      use CVAP
      use CGIBBS
      use CY
      use CA
      use CIPR
      IMPLICIT REAL*8(A-H, O-Z)      !C QUITAR                                    
      EXTERNAL STABIL, GMIX, FUNC                                         
      !common/CVAP/NOVAP, NDUM, IDUM(4), PRAT(10)                           
      !common/CGIBBS/NF, MAXZ, GNUL, Z(10), A(10), XVL(10, 4), SFAS(4),&
      !& GAM(10, 10), AL(10), DA(10, 10), XM(10, 4)                                       
      !common/CUFAC/NKK, NGG, Pxx(10, 10), Txx                                      
      !common/CY/Y13, Y21, STEP                                            
      !common/CA/XC(5), GE(5, 2), GC(5, 2)                                   
      !common/CIPR/IPR                                                   
      common/CQT/QT(10, 10), Q(10), R(10)                                  
      common/CACT/Y1(10), Y2(10), ACT1(10), ACT2(10), DACT1(10, 10),&
      & DACT2(10, 10), PACT(2, 2)                                                     
      common/CMODEL/MODEL                                               
      common/COUT/IOUT                                                  
      common/nga/nga, mass(12)
      common/ig/ig
      dimension DLX(10), YVAL(30), Y(10), GRAD(30), XMAT(30, 30), WORK(30, 5)  
      dimension X(2)         
      !integer::MODEL, IPR, IOUT, NOVAP, ig
      integer::MODEL, IOUT, ig            
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

      if (.not.(IOUT.EQ.6)) then 
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

      T1 = 0.                                                             
      NN = 0                                                              
      Txx = Tf
      PP = Pf
      if (Txx.EQ.0.) then
         if(IOUT.EQ.1) CLOSE (UNIT = 1)
         close (unit = 3)
         return
      endif

      if (PP/=  0..and.NOVAP/=  0) then
         do I = 1, NKK
         PRAT(I) = DLOG(PP)-ANT(1, I)+ANT(2, I)/(Txx-273.15+ANT(3, I))
         enddo                                                    
      endif

      Z(:) = Zf(:)
      ZSUM = 0.                                                           
      ZMAX = 0.                                                           
      do I = 1, NKK                                                       
         ZSUM = ZSUM+Z(I)                                                    
         if (Z(I).LT.ZMAX) cycle                                          
         ZMAX = Z(I)                                                         
         MAXZ = I
         continue
      enddo                                                            
      15 continue 

      if (.not.(Txx.EQ.T1)) then !goto 30                                               
         call PARAM2                                                       
         if (.not.(ICALC.NE.1)) then                                            
            if (NKK.NE.2.AND.NKK.NE.3) write(6, 616)                                
            if (IOUT.NE.6.AND.NKK.NE.2.AND.NKK.NE.3) write(IOUT, 616)               
            Y13 = Z(1)                                                          
            Y21 = Z(2)                                                                                                           
            if (IOUT.NE.6) write(IOUT, 633) Txx                                   
            if (.not.(NKK.EQ.3)) then                                                
               call SOLBIN                                                       
            if(IOUT.EQ.1) CLOSE (UNIT = 1)
            close (unit = 3)
            return
         endif   

         12 STEP = Z(3)/100.D0                                                  
            if (STEP.EQ.0.) STEP = .02D0                                         
            call BINOD                                                        
            if(IOUT.EQ.1) CLOSE (UNIT = 1)
            close (unit = 3)
            return
         endif 

         !16 continue

         if (.not.(ICALC.NE.2)) then                                           
            if (NKK.NE.2) write(6, 616)                                           
            if (IOUT.NE.6.AND.NKK.NE.2) write(IOUT, 616)                          
            XC(1) = 0.                                                          
            XC(2) = .2D0                                                        
            XC(3) = .5D0                                                        
            XC(4) = .8D0                                                        
            XC(5) = 1.D0                                                        
            do K = 1, 5                                                      
               Y(1) = XC(K)                                                        
               Y(2) = 1.D0-XC(K)                                                   
               call unifac(1, Y, ACT1, DACT1, PACT)                                  
               GE(K, 1) = ACT1(1)                                                   
               GE(K, 2) = ACT1(2)
            enddo                                                   
            READ(2, *) R(1), Q(1)                                               
            READ(2, *) R(2), Q(2)                                                                                    
            write(6, 627)                                                      
            do I = 1, 2                                                       
               write(6, 626) I, R(I), Q(I)   
            enddo
            if (.not.(IOUT.EQ.6)) then                                             
               write(IOUT, 627)                                                   
               do I = 1, 2
                  write(IOUT, 626) I, R(I), Q(I)
               enddo  
            endif

            13 continue                                                          
            X(1) = Z(1)/300.D0                                                  
            X(2) = Z(2)/300.D0                                                  
            do I = 1, 2                                                       
               do J = 1, 2                                                       
                  QT(I, J) = 0.                                                        
                  Pxx(I, J) = 0.                                                         
               enddo
            enddo
            QT(1, 1) = Q(1)                                                      
            QT(2, 2) = Q(2)                                                      
            NK = 2                                                              
            NGG = 2                                                              
            XLAMB = 1.                                                          
            call MARQ(FUNC, 2, 10, X, XLAMB, 3.D0, 1.D-7, 99)                        
            write(6, 633) Txx                                                    
            if (IOUT.NE.6) write(IOUT, 633) Txx                                   
            write(6, 617) Pxx(1, 2), Pxx(2, 1)       
            !(///, ' ** UNIQUAC PARAMETERS FROM UNIFAC **',
            !//, 5X, 'A12/R = ', F12.3, ' K , A21/R = ', F12.3, ' K', ///)                                  
            if (IPR.EQ.1) write(6, 618)                                         
            do L = 1, 5                                                       
               do I = 1, 2                                                       
                  GE(L, I) = DEXP(GE(L, I))                                             
                  GC(L, I) = DEXP(GC(L, I))
               enddo
            enddo 

            if (IPR.EQ.1) write(6, 619) ((GE(L, I), L = 1, 5), I = 1, 2)                 
            if (IPR.EQ.1) write(6, 619) ((GC(L, I), L = 1, 5), I = 1, 2)                 
      
            if (.not.(IOUT.EQ.6)) then                                             
               write(IOUT, 617) Pxx(1, 2), Pxx(2, 1)                                     
               if (IPR.EQ.1) write(IOUT, 618)                                      
               if (IPR.EQ.1) write(IOUT, 619) ((GE(L, I), L = 1, 5), I = 1, 2)              
               if (IPR.EQ.1) write(IOUT, 619) ((GC(L, I), L = 1, 5), I = 1, 2)
            endif              
   
            22 continue                                                          
            if(IOUT.EQ.1) CLOSE (UNIT = 1)
            close (unit = 3)
            return
         endif

         !19 continue                                                          
         do I = 1, NKK                                                       
            do J = 1, NKK                                                       
               GAM(I, J) = 0.D0                                                     
               if (.not.(J.EQ.I)) then                                               
                  call GAMINF(I, J, G)                                                
                  GAM(I, J) = G
               endif
            enddo                                                        
         enddo   

         !20 continue                                                          
      endif
      30 T1 = Txx                                                              
      NN = NN+1                                                                                                       
      do I = 1, NKK                                                       
         Z(I) = Z(I)/ZSUM  
      enddo

      if (IOUT.NE.6) write(IOUT, 602) NN                                  
      if (IOUT.NE.6) write(IOUT, 605) Txx, PP, ZSUM, (Z(I), I = 1, NKK)              
      call unifac(1, Z, AL, DA, PACT)                                       
      SFAS(1) = 1.                                                        
      GNUL = 0.                                                           
      do I = 1, NKK                                                       
         XVL(I, 1) = 1.                                                       
         Z(I) = Z(I)+1.D-20                                                  
         DLX(I) = DLOG(Z(I))                                                 
         A(I) = AL(I)+DLX(I)                                                 
         GNUL = GNUL+Z(I)*AL(I)
      enddo                                              
      NF = 1     

   50 call STIG(Y, S)                                                    
      !if (S.GT.-1.D-7) goto 70
      if (.not.(S.GT.-1.D-7)) then
         if (IOUT.NE.6) write(IOUT, 603)                                     
         do I = 1, NKK                                                       
            YVAL(I) = 1.D-5*Y(I)/Z(I)
         enddo
         60 continue                                                          
         !goto 100
      else
         70 continue
         do I = 1, NKK                                                       
            YVAL(I) = DLOG(Y(I))  
         enddo

         XLAM = 1.                                                                              
         if (IOUT.NE.6.AND.NF.EQ.1.AND.IPR.GT.0) write(IOUT, 606)            
         if (IOUT.NE.6.AND.NF.GT.1.AND.IPR.GT.0) write(IOUT, 609) NF         
         call TMSSJ(30, NKK, IPR, 15, XLAM, 1.D-12, FUN, YVAL, GRAD, XMAT, WORK, 1)     
  
         if (.not.(FUN.LT.-1.D-7)) then                                        
            write(output, *) 1
            write(output, 2613) (Z(j), J = 1, NKK)
            write(output, 2613) (AL(j), j = 1, NKK)        
            write(output, *) "SYSTEM IS STABLE"                                                   
	         write(7, 46) Txx, (xM(l, 1), l = 1, NKK)   
	         write(7, 46) Txx, (xM(l, 2), l = 1, NKK)
	         write(7, *)                    
            if (IOUT.NE.6) write(IOUT, 604)
         endif

         80 if (IOUT.NE.6) write(IOUT, 603)                                     
         do I = 1, NKK                                                       
            YVAL(I) = 1.D-5*DEXP(YVAL(I))/Z(I)   
         enddo
      endif
                       
  100 NF = NF+1

  !104 do 105 I = 1, NKK                                                      
  !      if (YVAL(I).GT.1.D0) goto 106                                      
  !105 continue                                                          
  !    goto 109                                                          
  !106 do 107 I = 1, NKK                                                      
  !107   YVAL(I) = YVAL(I)/10.                                               
  !    goto 104                                                          
  !109 continue
  
  outer:  do while (YVAL(I).GT.1.D0)
            do I = 1, NKK
                if (.not.(YVAL(I).GT.1.D0)) then
                  exit outer
               endif
            enddo
            do I = 1, NKK
               YVAL(I) = YVAL(I)/10.
            enddo
         enddo outer

      SFAS(NF) = 1.                                                       
      XLAM = .2                                                           
      if (NF.EQ.2) XLAM = .5                                               
      M = (NF-1)*NKK                                                                                          
      if (IOUT.NE.6.AND.IPR.GT.0) write(IOUT, 607) NF                     
      call TMSSJ(30, M, IPR, 60, XLAM, 1.D-16, FUN, YVAL, GRAD, XMAT, WORK, 2)     
      NT = NF*NKK                                                           
      NB = NT-NKK                                                           
      do I = 1, NB                                                     
         YVAL(NT+1-I) = YVAL(NB+1-I)
      enddo
                                                                                          
      NVAP = 0

      do J = 1, NF                                                     
         if (IDUM(J).EQ.1) NVAP = J
      enddo
                                                
  111 continue                                                                                         
      if (IOUT.NE.6.AND.NVAP.EQ.0) write(IOUT, 630)                       
      if (IOUT.NE.6.AND.NVAP.NE.0) write(IOUT, 631) NVAP                                                   
      if (IOUT.NE.6) write(IOUT, 614) NF                                  
      if (IOUT.NE.6) write(IOUT, 611)(J, SFAS(J), J = 1, NF)                   
      if (IOUT.NE.6) write(IOUT, 612) (J, J = 1, NF)                          
      SUM = 0.                                                            
      do I = 1, NKK                                                      
         DLX(I) = XVL(I, NF)*Z(I)/SFAS(NF)                                    
         SUM = SUM+DLX(I)
      enddo

      SUM = DLOG(SUM)                                                     
      call unifac(1, DLX, A, DA, PACT)                                      
      do I = 1, NKK                                                      
         DLX(I) = DLOG(DLX(I))                                               
         A(I) = A(I)+DLX(I)-SUM
      enddo

!-------------------------------------------------------------------------------
      do j = 1, nf
         do i = 1, NKK
            xmj(i) = xm(i, j)
         enddo
         call unifac(1, xmj, actgam, de, pe)
         do i = 1, NKK
            agam(i, j) = actgam(i)
         enddo
      enddo
      
      write (output, *) NF
      compfases(:, :) = xm(:, :)
      do i = 1, NF !escribe resultados para el output que sera leido por excel
        write(output, 2613) (XM(j, i), J = 1, NKK)
        write(output, 2613) (agam(j, i), j = 1, NKK)  
      enddo
    
      if (.not.(IOUT.EQ.6)) then                                            
         do I = 1, NKK                                                      
            write(IOUT, 613) I, (XM(I, J), J = 1, NF)    
            write(iout, 1613) i, (agam(i, j), j = 1, nf)
         enddo
         46	format (2X, F12.2, 8X, F12.6, 8X, F12.6 , 8X, F12.6, 8X, F12.6, &                                                 
         8X, F12.6, 8X, F12.6, 8X, F12.6, 8X, F12.6, 8X, F12.6, 8X,&
         F12.6, 8X, F12.6)
!-------------------------------------------------------------------------------
      endif
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
      

!###############################################################################                                                                 

!******************************** F I NKK ***************************************
