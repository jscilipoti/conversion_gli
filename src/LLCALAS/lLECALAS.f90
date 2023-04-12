    subroutine llecalas(Tf,Pf,Zf) !Temperatura, presion, composici�n, nro de componentes
!C  *******************************************************************  
!C  *                                                                 *  
!C  *           PROGRAM   L L E C A L A S  (asociacion incorporada    *  
!C  *                                      para el calculo de flash y *  
!C  *                                      curva binodal (icalc 0 y 1)*
!C  *                                                                 *
!C  *                          ASOCIACI�N CRUZADA					   *		  
!C  *                       VERSI�N GENERALIZADA                      *        
!C  *              FEBRERO 2006 MODIFICADA POR                        *
!C  *                   ALFONSINA  ESTER ANDREATTA                    *        
!C  *        BASADA EN LAS SIMPLLIFICACIONES DE LOS PAPERS:           *
!c  *       Revisada en Octubre del 2007 en el chequeo de estabilidad *
!C  *																   *   	
!c  * Michelsen, et al. (Fluid Phase Equilibria, 180(2001)165-174 )   *		
!C  * Tan, et al.  (Ind. Eng. Chem. Res, 2004,43,203-208).			   *	   	
!C  *																   *	   	
!C  *        Esto permiti�  que todos los casos particulares          *         
!c  *       de asociaci�n se puedan simplificar a un �nico c�lculo. 
!c
!c   V�lido para un m�ximo n�mero grupo asociativo de 12
!c   Con la implementaci�n en el c�lculo de la fracci�n no asociada en el componente puro 
!c   por  el metodo iterativo aqu� implementado se permite que una mol�cula
!c   tenga m�s de un grupo asociativo 14/07/06
!C  El c�lculo se limita a que el n�mero m�ximo de sitios sea 2(por razones matem�ticas)
!c                                                       
!C  *******************************************************************  
!C  *                                           DATE: 24/3 - 1982 /TJ *  
      use InputData
      use flashout
      IMPLICIT REAL*8(A-H,O-Z)                                          
      EXTERNAL STABIL,GMIX,FUNC                                         
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
      COMMON/CUFAC/N,NG,P(10,10),T                                      
      COMMON/CY/Y13,Y21,STEP                                            
      COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
      COMMON/CIPR/IPR                                                   
      COMMON/CQT/QT(10,10),Q(10),R(10)                                  
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CMODEL/MODEL                                               
      COMMON/COUT/IOUT                                                  
      common/nga/nga,mass(12)
      common/ig/ig
      DIMENSION DLX(10),YVAL(30),Y(10),GRAD(30),XMAT(30,30),WORK(30,5)  
      DIMENSION X(2)         
      integer::MODEL,IPR,IOUT,NOVAP,ig            
      character(len=36):: name1 
      integer:: parameters 
      real*8::Zf(10)
      dimension xmj(10),actgam(10),de(10,10),pe(2,2)
   
    !Variables le�das desde main se pasan a los commons de llcalas
    model = modelo
    IPR = iprm
    IOUT = ioutm
    NOVAP = novapm
    IG = igm

!   icalc:  0-' **** FLASH CALCULATION ****'                            
!           1-' **** BINODAL CURVE CALCULATION ****'
!           2-' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** '
!   model:  0-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'     
!           1-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'
!   ipr:    1-' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'
!   iout:   1-'open 'lleasoccuzada.OUT''
!   novap:  0-'VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS'
!   ig:     0-'write compositions'
!           1-'write compositions and activities'
!   ipareq: 1-'liquid-liquid parameters table (UNIFAC)'
!           2-'liquid-vapor parameters table (UNIFAC)'
!           3-'infinite dilution parameters table (UNIFAC)'
!           4-'GC-EOS parameters'
  
   
    !  call ab_ban1(model) ahora se llama desde el main
	  IF (IOUT.EQ.1) OPEN (UNIT=1,FILE='lleasoccuzada.OUT',FORM='FORMATTED')
	  OPEN (UNIT=3,FILE='output.OUT',FORM='FORMATTED')
	  output=3
      !WRITE(6,608)                                                      
      !WRITE(6,628) IOUT                                                 
      !WRITE(6,610)                                                      
      IF(IOUT.EQ.0) IOUT=6                                              
      !IF(ICALC.EQ.0) WRITE(6,620)                                       
      !IF(ICALC.EQ.1) WRITE(6,621)                                       
      !IF(ICALC.EQ.2) WRITE(6,622)                                       
      !IF(NOVAP.NE.0) WRITE(6,629)                                       
      !IF(MODEL.EQ.0) WRITE(6,624)                                       
      !IF(MODEL.EQ.1) WRITE(6,625)                                       
      !WRITE(6,623) NTEXT                                                
      IF(IOUT.EQ.6) GOTO 5                                              
      WRITE(IOUT,608)                                                   
      WRITE(IOUT,610)                                                   
      IF(ICALC.EQ.0) WRITE(IOUT,620)                                    
      IF(ICALC.EQ.1) WRITE(IOUT,621)                                    
      IF(ICALC.EQ.2) WRITE(IOUT,622)                                    
      IF(NOVAP.NE.0) WRITE(IOUT,629)                                    
      IF(MODEL.EQ.0) WRITE(IOUT,624)                                    
      IF(MODEL.EQ.1) WRITE(IOUT,625)                                    
      WRITE(IOUT,623) NTEXT                                             
    5 CONTINUE                                                          
    !  CALL PARIN2     hago que se llame desde el main                                                  
!      IF(NOVAP/=0) then                                             
!        DO 6 J=1,N                                                        
!!C   6 READ(2,502) (ANT(K,J),K=1,3)                                      
!    6       READ(2,*) (ANT(K,J),K=1,3)                                        
!        DO 7 J=1,N                                                        
!            ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
!    7       ANT(2,J)=2.302585*ANT(2,J)                                        
!      endif                                                          
      T1=0.                                                             
      NN=0                                                              
!C  10 READ(2,502) T,PP   se lee en el main                                               
!   10 IF (ICALC.EQ.0) READ(2,*) T,PP                                   
 !     IF (ICALC.GE.1) READ(2,*) T
      T = Tf
      PP = Pf
      IF(T.EQ.0.) GOTO 10000 
      !IF(PP.EQ.0..OR.NOVAP.EQ.0) then !reemplazado por la siguiente estructura
      IF(PP/=0..and.NOVAP/=0) then                                 
        DO 1 I=1,N                                                        
    1       PRAT(I)=DLOG(PP)-ANT(1,I)+ANT(2,I)/(T-273.15+ANT(3,I))            
!C   4 READ(2,502) (Z(I),I=1,N)                                          
      endif
  !  4 READ(2,*) (Z(I),I=1,N)         
      Z(:) = Zf(:)
      ZSUM=0.                                                           
      ZMAX=0.                                                           
      DO 15 I=1,N                                                       
        ZSUM=ZSUM+Z(I)                                                    
        IF(Z(I).LT.ZMAX) cycle !GOTO 15                                          
        ZMAX=Z(I)                                                         
        MAXZ=I                                                            
   15 CONTINUE                                                          
      IF(T.EQ.T1) GOTO 30                                               
      CALL PARAM2                                                       
      IF(ICALC.NE.1) GOTO 16                                            
      IF(N.NE.2.AND.N.NE.3) WRITE(6,616)                                
      IF(IOUT.NE.6.AND.N.NE.2.AND.N.NE.3) WRITE(IOUT,616)               
      Y13=Z(1)                                                          
      Y21=Z(2)                                                          
      !WRITE(6,633) T                                                    
      IF(IOUT.NE.6) WRITE(IOUT,633) T                                   
      IF(N.EQ.3) GOTO 12                                                
      CALL SOLBIN                                                       
      GOTO 10000                                                        
   12 STEP=Z(3)/100.D0                                                  
      IF(STEP.EQ.0.) STEP=.02D0                                         
      CALL BINOD                                                        
      GOTO 10000                                                        
   16 CONTINUE                                                          
      IF(ICALC.NE.2) GOTO 19                                            
      IF(N.NE.2) WRITE(6,616)                                           
      IF(IOUT.NE.6.AND.N.NE.2) WRITE(IOUT,616)                          
      XC(1)=0.                                                          
      XC(2)=.2D0                                                        
      XC(3)=.5D0                                                        
      XC(4)=.8D0                                                        
      XC(5)=1.D0                                                        
      DO 17 K=1,5                                                       
        Y(1)=XC(K)                                                        
        Y(2)=1.D0-XC(K)                                                   
        CALL unifac(1,Y,ACT1,DACT1,PACT)                                  
        GE(K,1)=ACT1(1)                                                   
   17   GE(K,2)=ACT1(2)                                                   
      READ(2,*) R(1),Q(1)                                               
      READ(2,*) R(2),Q(2)                                               
!C     READ(2,502) R(1),Q(1)                                             
!C     READ(2,502) R(2),Q(2)                                             
      WRITE(6,627)                                                      
      DO 14 I=1,2                                                       
   14   WRITE(6,626) I,R(I),Q(I)                                          
      IF(IOUT.EQ.6) GOTO 13                                             
      WRITE(IOUT,627)                                                   
      DO 11 I=1,2                                                       
   11   WRITE(IOUT,626) I,R(I),Q(I)                                       
   13 CONTINUE                                                          
      X(1)=Z(1)/300.D0                                                  
      X(2)=Z(2)/300.D0                                                  
      DO 18 I=1,2                                                       
        DO 18 J=1,2                                                       
            QT(I,J)=0.                                                        
   18       P(I,J)=0.                                                         
      QT(1,1)=Q(1)                                                      
      QT(2,2)=Q(2)                                                      
      NK=2                                                              
      NG=2                                                              
      XLAMB=1.                                                          
      CALL MARQ(FUNC,2,10,X,XLAMB,3.D0,1.D-7,99)                        
      WRITE(6,633) T                                                    
      IF(IOUT.NE.6) WRITE(IOUT,633) T                                   
      WRITE(6,617) P(1,2),P(2,1)       !(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                                  
      IF(IPR.EQ.1) WRITE(6,618)                                         
      DO 21 L=1,5                                                       
        DO 21 I=1,2                                                       
            GE(L,I)=DEXP(GE(L,I))                                             
   21       GC(L,I)=DEXP(GC(L,I))                                             
      IF(IPR.EQ.1) WRITE(6,619) ((GE(L,I),L=1,5),I=1,2)                 
      IF(IPR.EQ.1) WRITE(6,619) ((GC(L,I),L=1,5),I=1,2)                 
      IF(IOUT.EQ.6) GOTO 22                                             
      WRITE(IOUT,617) P(1,2),P(2,1)                                     
      IF(IPR.EQ.1) WRITE(IOUT,618)                                      
      IF(IPR.EQ.1) WRITE(IOUT,619) ((GE(L,I),L=1,5),I=1,2)              
      IF(IPR.EQ.1) WRITE(IOUT,619) ((GC(L,I),L=1,5),I=1,2)              
   22 CONTINUE                                                          
      GOTO 10000                                                        
   19 CONTINUE                                                          
      DO 20 I=1,N                                                       
        DO 20 J=1,N                                                       
            GAM(I,J)=0.D0                                                     
            IF(J.EQ.I) GOTO 20                                                
            CALL GAMINF(I,J,G)                                                
            GAM(I,J)=G                                                        
   20 CONTINUE                                                          
   30 T1=T                                                              
      NN=NN+1                                                           
      !WRITE(6,602) NN                                                   
      DO 35 I=1,N                                                       
   35 Z(I)=Z(I)/ZSUM                                                    
      !WRITE(6,605) T,PP,ZSUM,(Z(I),I=1,N)                               
      IF(IOUT.NE.6) WRITE(IOUT,602) NN                                  
      IF(IOUT.NE.6) WRITE(IOUT,605) T,PP,ZSUM,(Z(I),I=1,N)              
      CALL unifac(1,Z,AL,DA,PACT)                                       
      SFAS(1)=1.                                                        
      GNUL=0.                                                           
      DO 40 I=1,N                                                       
        XVL(I,1)=1.                                                       
        Z(I)=Z(I)+1.D-20                                                  
        DLX(I)=DLOG(Z(I))                                                 
        A(I)=AL(I)+DLX(I)                                                 
   40   GNUL=GNUL+Z(I)*AL(I)                                              
      NF=1                                                              
   50 CALL STIG(Y,S)                                                    
      IF(S.GT.-1.D-7) GOTO 70                                           
      !WRITE(6,603)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,603)                                     
      DO 60 I=1,N                                                       
        YVAL(I)=1.D-5*Y(I)/Z(I)                                           
   60 CONTINUE                                                          
      GOTO 100                                                          
   70 DO 75 I=1,N                                                       
   75   YVAL(I)=DLOG(Y(I))                                                
      XLAM=1.                                                           
      !IF(NF.EQ.1.AND.IPR.GT.0) WRITE(6,606)                             
      !IF(NF.GT.1.AND.IPR.GT.0) WRITE(6,609) NF                          
      IF(IOUT.NE.6.AND.NF.EQ.1.AND.IPR.GT.0) WRITE(IOUT,606)            
      IF(IOUT.NE.6.AND.NF.GT.1.AND.IPR.GT.0) WRITE(IOUT,609) NF         
      CALL TMSSJ(30,N,IPR,15,XLAM,1.D-12,FUN,YVAL,GRAD,XMAT,WORK,1)     
      IF(FUN.LT.-1.D-7) GOTO 80                                         
      !WRITE(6,604)         
      write(output,*) 1
        write(output,2613) (Z(j),J=1,N)
        write(output,2613) (AL(j),j=1,N)        
      write(output,*) "SYSTEM IS STABLE"                                                   

	 write(7,46) T,  (xM(l,1),l=1,N)    !Alfonsina
	 write(7,46) T,  (xM(l,2),l=1,N)!Alfonsina
	 write(7,*)                    !Alfonsina



      IF(IOUT.NE.6) WRITE(IOUT,604)                                     
    !  GOTO 10                                                           
   !80 WRITE(6,603)                                                      
   80 IF(IOUT.NE.6) WRITE(IOUT,603)                                     
      DO 90 I=1,N                                                       
   90   YVAL(I)=1.D-5*DEXP(YVAL(I))/Z(I)                                  
  100 NF=NF+1                                                           
  104 DO 105 I=1,N                                                      
        IF(YVAL(I).GT.1.D0) GOTO 106                                      
  105 CONTINUE                                                          
      GOTO 109                                                          
  106 DO 107 I=1,N                                                      
  107   YVAL(I)=YVAL(I)/10.                                               
      GOTO 104                                                          
  109 CONTINUE                                                          
      SFAS(NF)=1.                                                       
      XLAM=.2                                                           
      IF(NF.EQ.2) XLAM=.5                                               
      M=(NF-1)*N                                                        
      !IF(IPR.GT.0) WRITE(6,607) NF                                      
      IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,607) NF                     
      CALL TMSSJ(30,M,IPR,60,XLAM,1.D-16,FUN,YVAL,GRAD,XMAT,WORK,2)     
      NT=NF*N                                                           
      NB=NT-N                                                           
      DO 110 I=1,NB                                                     
  110   YVAL(NT+1-I)=YVAL(NB+1-I)                                         
      !WRITE(6,614) NF                                                   
      NVAP=0                                                            
      DO 111 J=1,NF                                                     
        IF(IDUM(J).EQ.1) NVAP=J                                           
  111 CONTINUE                                                          
      !IF(NVAP.EQ.0) WRITE(6,630)                                        
      !IF(NVAP.NE.0) WRITE(6,631) NVAP                                   
      IF(IOUT.NE.6.AND.NVAP.EQ.0) WRITE(IOUT,630)                       
      IF(IOUT.NE.6.AND.NVAP.NE.0) WRITE(IOUT,631) NVAP                  
      !WRITE(6,611) (J,SFAS(J),J=1,NF)                                   
      !WRITE(6,612) (J,J=1,NF)                                           
      IF(IOUT.NE.6) WRITE(IOUT,614) NF                                  
      IF(IOUT.NE.6) WRITE(IOUT,611)(J,SFAS(J),J=1,NF)                   
      IF(IOUT.NE.6) WRITE(IOUT,612) (J,J=1,NF)                          
      SUM=0.                                                            
      DO 115 I=1,N                                                      
        DLX(I)=XVL(I,NF)*Z(I)/SFAS(NF)                                    
  115   SUM=SUM+DLX(I)                                                    
      SUM=DLOG(SUM)                                                     
      CALL unifac(1,DLX,A,DA,PACT)                                      
      DO 120 I=1,N                                                      
        DLX(I)=DLOG(DLX(I))                                               
  120   A(I)=A(I)+DLX(I)-SUM                                              
!c-----
      do 1130 j=1,nf
        do 1131 i=1,n
 1131       xmj(i)=xm(i,j)
        call unifac(1,xmj,actgam,de,pe)
        do 1132 i=1,n
 1132       agam(i,j)=actgam(i)
 1130 continue
      write (output,*) NF
      compfases(:,:) = xm(:,:)
      do i=1,NF !escribe resultados para el output que ser� le�do por excel
        write(output,2613) (XM(j,i),J=1,N)
        write(output,2613) (agam(j,i),j=1,N)  
      enddo
!      do i=1,NF !escribe resultados para el output que ser� le�do por excel
!        write(output,2613) (agam(j,i),j=1,N)      
!      enddo
      
      !DO 130 I=1,N                                                      
        !WRITE(6,613) I,(XM(I,J),J=1,NF)     !composition        
  !130   write(6,1613) i,(agam(i,j),j=1,nf) !Ln(gamma)
     
      IF(IOUT.EQ.6) GOTO 132                                            
      DO 131 I=1,N                                                      
      WRITE(IOUT,613) I,(XM(I,J),J=1,NF)    
  131 write(iout,1613) i,(agam(i,j),j=1,nf)
   46	FORMAT (2X,F12.2, 8X,F12.6, 8X,F12.6 , 8X,F12.6, 8X,F12.6, & !Alfonsina                                                 
     8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6)
!c-----
  132 CONTINUE                                                          
  !   GOTO 50                                                           


  501 FORMAT(36A2)                                                      
  502 FORMAT(8F10.2)                                                    
  503 FORMAT(20I3)                                                      
  602 FORMAT(///,' * * * FLASH NUMBER',I3,' * * *',//)                  
  603 FORMAT(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')            
  604 FORMAT(/,' * SYSTEM IS STABLE *',/)                               
  605 FORMAT(' TEMPERATURE =',F10.4,' K, PRESSURE =',F7.3,' ATM, FEED ='&
     ,F10.2,' MOLES',/,' FEED COMPOSITION (MOLE PERCENT):',/,1X,15(2PF7&
     .3))                                                              
  606 FORMAT(//,' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')       
  607 FORMAT(/,' PHASE SPLIT CALCULATION,',I2,' PHASES:')               
  608 FORMAT(1H1)                                                       
  609 FORMAT(//,' DIFFERENTIAL STABILITY TEST FOR',I2,'-PHASE SYSTEM')  
  610 FORMAT(///)                                                       
  611 FORMAT(/,'  PHASE FRACTIONS (PERCENT):',4(5X,I3,2PF7.3,5X))       
  612 FORMAT(/,'  COMPOSITION  ',10X,4(8X,I3,9X))                       
  613 FORMAT('   X(',I2,')            ',5(8X,F12.8))                    
!c-----
 1613 format('  ln(G',i2,')            ',5(8x,f12.8))
 2613 format(5(2x,f12.8))

!c-----
  614 FORMAT(//,' RESULT OF',I2,'-PHASE CALCULATION:')                  
  616 FORMAT(//,' * WRONG INPUT SPECIFICATION *',//)                    
  617 FORMAT(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                         
  618 FORMAT(//,' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'//)                                       
  619 FORMAT(10F12.5)                                                   
  620 FORMAT(' **** FLASH CALCULATION ****')                            
  621 FORMAT(' **** BINODAL CURVE CALCULATION ****',//)                 
  622 FORMAT(' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** ',//)                                                              
  623 FORMAT(1X,'COMPONENTS : ',40A2,//)                                
  624 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'//)     
  625 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'//)    
  626 FORMAT(I5,2F15.4)                                                 
  627 FORMAT(//,' SPECIFIED UNIQUAC R AND Q',/)                         
  628 FORMAT(/,' IOUT = ',I2,/' IF IOUT = 0: OUTPUT ONLY ON UNIT 6',/,  ' IF IOUT = 1: OUTPUT ON BOTH UNIT 6 AND 1')                      
  629 FORMAT(/,' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS',//)        
  630 FORMAT(' NO VAPOR PHASE')                                         
  631 FORMAT(' PHASE',I2,' IS A VAPOR PHASE')                           
  633 FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')                     
10000 continue
      IF (IOUT.EQ.1) CLOSE (UNIT=1)
      close (unit=3)
!c      call salida(name)
                                                                    
      endsubroutine llecalas   
                                                                  
      SUBROUTINE unifac(NDIF,X,ACT,DACT,PACT)                           
      IMPLICIT REAL*8(A-H,O-Z)                                          
!c------
      common/asoc/nktt,igamt(20,12),nytt(20,12)   
      common/nga/nga,mass(12)
      common/grupas1/rkass(6,12,6,12),enass(6,12,6,12),deloh(6,12,6,12)!Alfonsin
!c------
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
      COMMON/CUFAC/NK,NG,P(10,10),T                                     
      COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
      COMMON/CQT/QT(10,10),Q(10),R(10)                                  
      DIMENSION X(10),GAM(10),ACT(10),DACT(10,10),THETA(10),PHI(10),RI(10),QI(10),ETA(10),QIL(10),RIL(10),QID(10),ETAL(10),TETAR(10) 
      DIMENSION U(10,10),V(10,10),PACT(2,2),DTAU(2,2,2)                 
!c------
      dimension goh(10),xgamk(20),dxohdx(10),dxxdx(10,10),dasdx1(10,10),dasdx2(10,10),dasdx(10,10)
	common/ioh2/rngoh(12,12)

      dimension dif(12,12), dif1(10,12,12) !Alfonsina
      common/zzzas/xoh(6,12),xohi0(12,6,12),xoh_old(6,12),xohi(6,12),xohi_old(6,12), xohi0_old(12,6,12)  !Alfonsina
	dimension m_lambda(nga*2,nga*2),m_lambda1(nga*2,nga*2) !Alfonsina
	dimension psin(12) !Alfonsina
	dimension indx(20)
      double precision  m_lambda,m_lambda1,xoh,xohi0,xoh_old,xohi0_old  !Alfon
	double precision del, del1, dif, dif1, d1,psin, xgam, xnoh1 !Alfonsina
	integer order !Alfonsina
      double precision sum1, sum2, sum3, sum4, SUMA1J, sumaj !Alfonsina
      dimension xnohi0(12,12),tgt(12),dnohdx(12,12),actas(12) !Alfonsina
	dimension xnoh1(12), xnoh(12),das1(3),das3(3),dxkdni(12,6,12), dxkdnic(12,6,12), dgasdx(12)  !Alfonsina
      dimension dgasdxij (12,12), drhodx(12), drhodni(12,6,12)
	


      dk=1.381e-23
      deloh=0.0
      xnoh=0.0
      xnoh1=0.0
	xoh=0.0
      xgam=0.0
      do 7777 i=1,10
	xohi0=0
      xnohi0=0.0
      tgt(i)=0.0
 7777 continue

      THETS=0.                                                          
      PHS=0.                                                            
      DO 10 I=1,NK                                                      
      THETA(I)=X(I)*Q(I)                                                
      PHI(I)=R(I)*X(I)                                                  
      THETS=THETS+THETA(I)                                              
   10 PHS=PHS+PHI(I)                                                    
      DO 20 I=1,NK                                                      
      RI(I)=R(I)/PHS                                                    
      RIL(I)=DLOG(RI(I))                                                
      QI(I)=Q(I)/THETS                                                  
   20 QIL(I)=DLOG(QI(I))                                                

      do 33 i=1,nk
      goh(i)=0.
      tgt(i)=0.0
      xnohi0=0.0
      xgam=0.0
     
!CCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do j=1,nk
		tgt(j)=0.0
		xgam=0.0
	end do
   33 continue
      do i=1,nk
	if(nga.gt.0) then
	
		do k=1,nktt
			tgt(i)=tgt(i)+nytt(k,i)
		end do

		do j=1,nga  
			xnohi0(i,j)=rngoh(i,j)/R(i)  
		end do
		xgam=xgam+R(i)*x(i)
		end if  
	end do
  
      xnoh1=0d0
	do ja=1,nga
	    do i=1,nk
		xnoh1(ja)=xnoh1(ja)+rngoh(i,ja)*x(i)
	  end do
      end do
	
      do ja=1,nga
		xnoh(ja)=xnoh1(ja)/xgam
	end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c------
      DO 40 I=1,NG                                                      
      ETA(I)=0.                                                         
      DO 45 J=1,NK                                                      
   45 ETA(I)=ETA(I)+S(I,J)*X(J)                                         
   40 ETAL(I)=DLOG(ETA(I))                                              
      DO 55 I=1,NG                                                      
      TETAR(I)=0.                                                       
      DO 55 J=1,NK                                                      
   55 TETAR(I)=TETAR(I)+QT(I,J)*X(J)                                    
      DO 60 I=1,NK                                                      
      QID(I)=1.-RI(I)/QI(I)                                             
      XX=F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                             
      XX=XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                              
      ACT(I)=XX                                                         
      DO 661 J=1,NG                                                     
      U(J,I)=S(J,I)/ETA(J)                                              
      V(J,I)=U(J,I)*TETAR(J)                                            
  661 ACT(I)=ACT(I)-V(J,I)-QT(J,I)*ETAL(J)                              
!c------
  
  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	!**************************calculo de las fuerzas de asociacion*******************
      if(nga.ne.0) then
	    DO J=1,NGA 
              IF(MASS(J).EQ.0) GO TO 201
                DO m=1,NGA
                  IF(MASS(m).EQ.0) GO TO 101
                        DO L=1,MASS(J)
                                DO K=1,MASS(m)
                                  IF(ENASS(K,m,L,J).EQ.0) THEN
                                         CONTINUE
                                    ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      deloh(k,m,l,j)=(DEXP(ENASS(K,m,L,J)/T) - 1 )*RKASS(K,m,L,J)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                                  END IF
                                END DO
                        END DO

 101            CONTINUE
                END DO
 201          CONTINUE
        END DO
  

	end if  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	!***********************************calculo de xoh**************************
!c C�lculo de la  fracci�n no asociada Paper:Ind. Eng. Chem. Res, 2004,43,203-208   
!c !Inicializaci�
      if(nga.ne.0) then
      xoh=1.0d0 
      del=1.d0
!c Iteraciones con tolerancia de 10-9
      do while (del>1.0d-10)
      xoh_old=xoh
	do m=1, nga
		do j=1,2
      sum1=0.D0
	do k=1, nga
	sum2=0.D0
	do l=1,2

	sum2=sum2+xoh_old(l,k)*deloh(l,k,j,m)
	end do
	sum1=sum1+sum2*xnoh1(k)
      end do
      xoh(j,m)=1.D0/(1.D0+sum1/xgam)          
	dif(j,m)=dabs((xoh(j,m)-xoh_old(j,m))/xoh(j,m))
	end do
	end do
	del=maxval(dif)
      end do
      end if

!	write(4,*)"T=", t           
!	write(4,*)"xoh (1,1)=", xoh (1,1)  
!	write(4,*)"xoh (1,2)=", xoh(2,1) 
 
	 

!cc Fin del C�lculo de la  fracci�n no asociada 
!c	!*****************************calculo de xohi0**************************************
!C	xohi0=1d0
!c	do i=1, nc
!C		do j=1, nga
!C	       do l=1,mass(j)
!C	do k=1,nga
!C	 do m=1,mass(k)
!C			If (rngoh(i,j).eq.0d0) then
!C					xohi0(i,l,j)=1.d0
!C			elseif (deloh(l,j,m,k).gt.0d0.and.rngoh(i,j).ne.0d0.and.
!C     @		mass(j).eq.2)then
!C		xohi0(i,l,j)=(-1d0+dsqrt(1d0+4d0*xnohi0(i,j)*deloh(l,j,m,k)))/
!C     @			(2d0*xnohi0(i,j)*deloh(l,j,m,k))
!C	
!C			end if
!C		 end do
!C	end do
!C	end do
!C		end do
!c	end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c****************************calculo de xohi0(esta implementaci�n permite que una mol�cula
!c tenga m�s de un grupo asociativo 14/07/06**************************************
!c !Inicializaci�
     
      if(nga.ne.0) then
      xohi0=1.D0
      del1=1.D0
!c Iteraciones con tolerancia de 10-12
      do while (del1>1.0d-10)
      xohi0_old=xohi0
	do m=1, nga
	if	(rngoh(i,m).gt.0d0) then
		do j=1,2
      sum3=0.D0
	do k=1, nga
	sum4=0.D0
	do l=1,2 
	sum4=sum4+ xohi0_old(i,l,k)*deloh(l,k,j,m)*xnohi0(i,k)
	end do
	sum3=sum3+sum4
      end do
      xohi0(i,j,m)=1.D0/(1.D0+sum3)    
 	dif1(i,j,m)=dabs((xohi0(i,j,m)-xohi0_old(i,j,m))/xohi0(i,j,m))
	end do
	else
	end if
	end do
	del1=maxval(dif1)
      end do
      end if

!c*****************************fin del calculo de xohi0**************************************
!C�lculo del gama de asociaci�n ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C actas(M) = LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I  

       if(nga.ne.0) then	

      SUMAJ = 0.D0
      DO J=1,NGA 
      IF(MASS(J).NE.0) THEN      
      DO K=1,MASS(J)
	If(XOH(K,J).gt.1d-13)then
      SUMAJ = SUMAJ + RNGOH(i,j)*(dlog(XOH(K,J)/XOHI0(I,K,J))+0.5D0*(XOHi0(i,K,J)-1))+0.5D0*R(i)*xnoh(j)*(1-xoh(k,j))
	end if
      END DO
      ELSE
      CONTINUE
      END IF
      END DO
      actas(I) = SUMAJ
	end if
      act(i)=act(i)+actas(i)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   60 continue
      NDUM=0                                                            
      IF(NOVAP.EQ.0) GOTO 69                                            
      SS=0                                                              
      DO 61 I=1,NK                                                      
   61 SS=SS+X(I)*(PRAT(I)-ACT(I))                                       
      IF(SS.GT.0.) GOTO 69                                              
      NDUM=1                                                            
      DO 62 I=1,NK                                                      
      ACT(I)=PRAT(I)                                                    
      DO 62 J=1,NK                                                      
   62 DACT(I,J)=0.                                                      
      GOTO 100                                                          
   69 CONTINUE                                                          
      IF(NDIF.EQ.4) GOTO 90                                             
      IF(NDIF.LT.2) GOTO 100                                            
      DO 70 I=1,NK                                                      
      DO 70 J=I,NK                                                      
      XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))         
      DO 75 K=1,NG                                                      
   75 XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)                      

!********************************calculo de dxkdni Alfonsina**************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!cCalcula los elementos de la matriz deltapq para el c�lculo de la derivada de la fracci�n 
!c no asociada respecto a la fracci�n molar del componente
      psin=0.0d0
	if(nga.ne.0) then
	m_lambda1=0.0d0
	m_lambda=0.0d0
      z=0; y=0
	do n=1,2
	do m=1,nga
     	z=z+1
	do l=1, 2
	do k=1, nga
		y=y+1
      m_lambda(z,y)=xnoh(k)*deloh(l,k,n,m)*xoh(n,m)**2 
	   if (z.eq.y)  then
      m_lambda(z,y)=m_lambda(z,y)+ 1.0d0
      end if
	end do
	end do
	y=0
	end do
	end do
      order=nga*2
	end if 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!c Calculo de  los elementos de la matriz [Yp] para el c�lculo de la derivada de la fracci�n 
!c no asociada respecto a la fracci�n molar del componente
      if(nga.ne.0) then
      do k=1,nga

	do ll=1,2
	do m=1,nga
	drhodni(j,ll,m)=(((rngoh(j,m)*xgam-xnoh1(m)*R(j)))/xgam**2)	  
	end do
	end do

      z=0     
	do ll=1,2
	do m=1,nga
	sum3=0.0d0
      do l=1,nga
	sum4=0.0d0


	do kk=1,2
	sum4=sum4+ (xoh(kk,l)*deloh(kk,l,ll,m))*drhodni(j,kk,l)
	end do
	sum3=sum3+sum4
	end do
	z=z+1
	psin(z)=-(xoh(ll,m)**2)*sum3
	end do
	end do


      N=order
	NP=order
      m_lambda1=m_lambda
      call  ludcmp(m_lambda1,N,NP,indx,d1)
       call lubksb(m_lambda1,N,NP,indx,psin)
!c colectando las derivadas en su correspondiente sub�ndice
      z=0
	do m=1,2
	do l=1, nga
	z= z+1
	dxkdni(k,m,l)=psin(z)
	end do
	end do 
	end do


	do l=1,nga
	do m=1,2
	do kk=1,nga
	if (rngoh(i,kk).ne.0) then
      dxkdnic(j,m,l)=dxkdni(kk,m,l)   
	end if
	end do
      end do
	end do

!c fin del c�lculo de la derivada de la fracci�n no asociada respecto a la 
!c fracci�n molar del componente
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C dgasdx(M) = derivada LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I
      if(nga.ne.0) then
   
      DO l=1,NGA 
	   SUMA1J = 0.D0

      IF(MASS(l).NE.0) THEN      
      DO K=1,MASS(l)
	If(XOH(K,l).gt.1d-13)then

      SUMA1J = SUMA1J + RNGOH(i,l)*1.D0/XOH(K,l)*dxkdnic(j,k,l)+0.5D0*&
     r(i)*(drhodni(j,k,l)-xnoh(l)*dxkdnic(j,k,l)-drhodni(j,k,l)*&
     XOH(K,l))


	end if
      END DO
      ELSE
      CONTINUE
      END IF
	xx=xx+ SUMA1J
 	end do   
	end if
		end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DACT(I,J)=XX                                                      
   70 DACT(J,I)=XX    
       IF(NDIF.LT.3) GOTO 100                                           
      DO 80 I=1,NK                                                      
      GAM(I)=DEXP(ACT(I))                                               
   80 ACT(I)=GAM(I)*X(I)                                                
      DO 85 I=1,NK                                                      
      DO 85 J=1,NK                                                      
      DACT(I,J)=ACT(I)*(DACT(I,J)-1.D0)                                 
      IF(J.EQ.I)DACT(I,J)=DACT(I,J)+GAM(I)                              
   85 CONTINUE                                                          
      GOTO 100                                                          
   90 CONTINUE                                                          
      DO 91 I=1,2                                                       
      DO 91 K=1,2                                                       
      DTAU(I,K,K)=0.                                                    
      DO 91 L=1,2                                                       
      IF(L.EQ.K) GOTO 91                                                
      H1=TETAR(L)-QT(L,I)*ETA(L)/S(L,I)                                 
      H2=QT(K,I)-S(L,I)*TETAR(K)/ETA(L)                                 
      DTAU(I,K,L)=-H1*H2/ETA(L)                                         
   91 CONTINUE                                                          
      DO 92 I=1,NK                                                      
      PACT(I,1)=-DTAU(I,1,2)*TAU(1,2)/T*300.D0                          
   92 PACT(I,2)=-DTAU(I,2,1)*TAU(2,1)/T*300.D0                          
  100 RETURN                                                            
      END  
	
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      SUBROUTINE PARIN2           
      use InputData                                      
      IMPLICIT REAL*8(A-H,O-Z)    
    !PARAMETER(NMS=2) !,NCOM=3,NSCM=10)                                      
!c------
      common/asoc/nktt,igamt(20,12),nytt(20,12)  !,p4,p5
      common/grupas1/rkass(6,12,6,12),enass(6,12,6,12), deloh(6,12,6,12)!Alfon
      common/nga/nga,mass(12) 
		common/ioh2/rngoh(12,12)

	common/ioh2sis/rngoht(3,2)
!c------
      COMMON/CUFAC/NK,NG,P(10,10),T                                     
      COMMON/CQT/QT(10,10),Q(10),R(10)                                  
      COMMON/CMODEL/MODEL                                               
      COMMON/COUT/IOUT                                                  
      DIMENSION RT(10,10),A(100,100),NGM(10) !,MAINSG(57)                   
      DIMENSION MS(10,10,2),NY(10,20),JH(150),IH(20)  

    DIMENSION NPUNTA(NMG),NGRUPA(NMG),NPUNT(NMG),NGRUP(NMG) !xnoh1(12),&
    !RR(NMG),
!    ,ENASS(6,12,6,12),&
!    RKASS(6,12,6,12),deloh(6,12,6,12),R(NC),xnohi0(12,12),&
!    X(NC),xnoh(12),dif(12,12),xoh_old(6,12),xoh(6,12),&
!    xohi0(12,6,12),xohi0_old(12,6,12),xxohi0(12,6,12), dif1(10,12,12),&
!    act(NC),actas(NC),NGRUPASD(NMG),NPUNTASD(NMG)
    INTEGER CS(NMG),TS(NMG,NMS),NPUNTMG(NMG),NMAINGR(NMG),J,K,CANIL(10) !,A,B,rngoh(NC,NSCM),GA,&
             !
        real*8::par             
      logical VL
      external mainsgfunc
                         
      REAL*4 RR(150),QQ(150)                                              
      REAL*4    A1(32),A2(32),A3(32),A4(32),A5(32),A6(32),A7(32),A8(32),&
     A9(32),A10(32),A11(32),A12(32),A13(32),A14(32),A15(32),A16(32),A17(32)&
     ,A18(32),A19(32),A20(32),A21(32),A22(32),A23(32),A24(32),A25(32)&
     ,A26(32),A27(32),A28(32),A29(32),A30(32),A31(32),A32(32)        
   !   DATA MAINSG/4*1,4*2,2*3,3*4,5,6,7,8,9,2*10,11,12,2*13,2*14,4*15,3*16,3*17,2*18,19,20,2*21,22,3*23,24,25,26,3*27,28,29,30,31,32/     
      DATA RR/.9011,.6744,.4469,.2195,1.3454,1.1167,.8886,1.1173,.5313,.3652,&
      1.2663,1.0396,.8121,1.,3.2499,3.2491,.92,.8952,1.6724,1.4457,&
     .998,3.168,1.3013,1.528,1.9031,1.6764,1.145,.9183,.6908,.9183,1.4654,&
     1.238,1.006,2.2564,2.0606,1.8016,2.87,2.6401,3.39,1.1562,1.870132&
     ,1.6434,1.06,2.0086,1.7818,1.5544,1.4199,2.4088,4.0013,2.9993,2.83&
     ,2.667,3.3092,2.4317,3.0856,4.0358,2.8266,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0/                    
      DATA QQ/.848,.54,.228,0.,1.176,.867,.676,.988,.4,.12,.968,.66,.348&
     ,1.2,3.128,3.124,1.4,.68,1.488,1.18,.948,2.484,1.224,1.532,1.728,1.42&
     ,1.088,.78,.468,1.1,1.264,.952,.724,1.988,1.684,1.448,2.41,2.184&
     ,2.91,.844,1.724,1.416,.816,1.868,1.56,1.248,1.104,2.248,3.568,2.113&
     ,1.833,1.553,2.86,2.192,2.736,3.2,2.472,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
     0.0,0.0,0.0/                        
      DATA A1/0.,292.3,156.5,104.4,328.2,-136.7,-131.9,342.4,-159.8,66.56&
     ,146.1,14.78,1744.,-320.1,1571.,73.8,27.9,21.23,89.97,-59.06,29.08&
     ,175.8,94.34,193.6,108.5,81.49,-128.8,147.3,-11.91,14.91,67.84,36.42/                                                              
      DATA A2/74.54,0.,-94.78,-269.7,470.7,-135.7,-135.7,220.6,1.,306.1,&
     517.,1.,-48.52,485.6,76.44,-24.36,-52.71,-185.1,-293.7,1.,34.78,1.&
     ,375.4,5*1.,176.7,132.1,42.73,60.82/                              
      DATA A3/-114.8,340.7,0.,-146.8,-9.21,-223.,-252.,372.8,-473.2,-78.31&
     ,-75.3,-10.44,75.49,114.8,52.13,4.68,1.,288.5,-4.7,777.8,56.41,-218.9&
     ,113.6,7.18,247.3,-50.71,-255.3,1.,-80.48,-17.78,59.16,29.77/
      DATA A4/-115.7,4102.,167.,0.,1.27,-162.6,-273.6,203.7,-470.4,-73.87&
     ,223.2,-184.9,147.3,-170.,65.69,122.9,1.,33.61,134.7,-47.13,-53.29&
     ,-15.41,-97.05,-127.1,453.4,-30.28,-124.6,3*1.,26.59,55.97      /
      DATA A5/644.6,724.4,703.9,4000.,0.,-281.1,-268.8,-122.4,-63.15,216.&
     ,-431.3,444.7,118.4,180.6,137.1,455.1,669.2,418.4,713.5,1989.,2011.&
     ,529.,483.8,332.6,-289.3,-99.56,-319.2,837.9,4*1.              /
      DATA A6/329.6,1731.,511.5,136.6,937.3,2*0.,247.,-547.,401.7,643.4,&
     -94.64,728.7,-76.64,-218.1,351.5,-186.1,-465.7,-260.3,3*1.,264.7,9*1./                                                              
      DATA A7/310.7,1731.,577.3,906.8,991.3,2*0.,104.9,-547.2,-127.6,231.4&
     ,732.3,349.1,-152.8,-218.1,351.5,-401.6,-465.7,512.2,3*1.,264.7,9*1./                                                             
      DATA A8/1300.,896.,859.4,5695.,28.73,-61.29,5.89,0.,-595.9,634.8,623.7&
     ,211.6,652.3,385.9,212.8,770.,740.4,793.2,1205.,390.7,63.48,-239.8&
     ,13.32,439.9,-424.3,1.,203.,1153.,-311.,-262.6,1.11,1.       /
      DATA A9/2255.,1.,1649.,292.6,-195.5,-153.2,-153.2,344.5,0.,-568.,3*1.&
     ,-337.3,4*1.,1616.,2*1.,-860.3,1.,-230.4,523.,1.,-222.7,5*1.  /
      DATA A10/472.6,343.7,593.7,916.7,67.07,-47.41,353.8,-171.8,-825.7,&
     0.,128.,48.93,-101.3,58.84,52.38,483.9,550.6,342.2,550.,190.5,-349.2&
     ,857.7,377.,211.6,82.77,2*1.,417.4,4*1.          /              
      DATA A11/158.1,-214.7,362.3,1218.,1409.,-344.1,-338.6,-349.9,1.,-37.36,&
      0.,-311.6,1051.,1090.,1.,-47.51,16*1./                       
      DATA A12/383.,1.,31.14,715.6,-140.3,299.3,-241.8,66.95,1.,120.3,1724.&
     ,0.,-115.7,-46.13,2*1.,808.8,203.1,70.14,5*1.,-75.23,1.,-201.9,123.2,1.,&
     -281.9,2*1./                                             
      DATA A13/139.4,1647.,461.8,339.1,-104.,244.4,-57.98,-465.7,1.,1247.&
     ,.75,1919.,0.,1417.,1402.,337.1,437.7,370.4,438.1,1349.,1.,681.4,&
     152.4,1.,-1707.,2*1.,639.7,4*1.          /                        
      DATA A14/972.4,-577.5,6.,5688.,195.6,19.57,487.1,-6.32,-898.3,258.70&
     ,-245.8,57.7,-117.6,0.,461.3,1.,-132.9,176.5,129.5,-246.3,2.41,3*1.,&
     29.86,7*1./                      
      DATA A15/662.1,289.3,32.14,213.1,262.5,1970.,1970.,64.42,1.,5.202,2*1.&
     ,-96.62,-235.7,0.,225.4,-197.7,-20.93,113.9,3*1.,-94.49,9*1. /
      DATA A16/42.14,99.61,-18.81,-114.1,62.05,-166.4,-166.4,315.9,1.,1000.&
     ,751.8,1.,19.77,1.,301.1,0.,-21.35,-157.1,11.8,13*1.       /   
      DATA A17/-243.9,337.1,2*1.,272.2,128.6,507.8,370.7,1.,-301.,1.,-347.9,&
      1670.,108.9,137.8,110.5,0.,1.,17.97,13*1.         /           
      DATA A18/7.5,4583.,-231.9,-12.14,-61.57,2*1544.,356.8,1.,12.01,1.,&
     -249.3,48.15,-209.7,-154.3,249.2,1.,0.,51.9,1.,-15.62,-216.3,4*1.,&
     -114.7,5*1. /                                                     
      DATA A19/-5.55,5831.,3000.,-141.3,-41.75,224.6,-207.,502.9,4894.,&
     -10.88,1.,61.59,43.83,54.57,47.67,62.42,56.33,-30.1,0.,-255.4,-54.86,&
     8455.,-34.68,514.6,8*1. /                                       
      DATA A20/924.8,1.,-878.1,-107.3,-597.1,2*1.,-97.27,1.,902.6,2*1.,&
     874.3,629.,4*1.,475.8,0.,-465.2,1.,794.4,1.,-241.7,1.,-906.5,5*1. /
      DATA A21/696.8,405.9,29.13,1208.,-189.3,2*1.,198.3,1.,430.6,3*1.,&
      -149.2,3*1.,70.04,492.,346.2,0.,5*1.,-169.7,5*1.        /          
      DATA A22/902.2,1.,1.64,689.6,-348.2,1.,1.,-109.8,-851.6,1010.,2*1.,&
      942.2,4*1.,-75.5,1302.,2*1.,0.,1.,175.8,164.4,1.,-944.9,5*1./    
      DATA A23/556.7,425.7,-1.77,3629.,-30.7,150.8,150.8,1538.6,1.,400.,&
      2*1.,446.3,1.,95.18,3*1.,490.9,-154.5,2*1.,0.,1.,481.3,7*1.      /
      DATA A24/575.7,1.,-11.19,-175.6,-159.,2*1.,32.92,-16.13,-328.6,8*1.,&
      534.7,2*1.,179.9,1.,0.,-246.,7*1.           /                   
      DATA A25/527.5,1.,358.9,337.7,536.6,2*1.,-269.2,-538.6,211.6,1.,&
      -278.2,572.7,343.1,5*1.,124.8,1.,125.3,139.8,963.,0.,7*1.    /      
      DATA A26/269.2,1.,363.5,1023.,53.37,20*1.,0.,6*1.          /      
      DATA A27/-300.,1.,-578.2,-390.7,183.3,2*1.,-873.6,-637.3,2*1.,-208.4,&
      5*1.,18.98,1.,-387.7,134.3,924.5,4*1.,0.,5*1.     /            
      DATA A28/-63.6,3*1.,-44.44,2*1.,1429.,1.,148.,1.,-13.91,-2.16,14*1.,0.,&
      4*1./                                                        
      DATA A29/928.3,500.7,364.2,4*1.,-364.2,20*1.,0.,3*1.    /         
      DATA A30/331.,115.4,-58.1,4*1.,-117.4,3*1.,173.8,17*1.,0.,2*1.   /
      DATA A31/561.4,784.4,21.97,238.,3*1.,18.41,22*1.,0.,1.       /    
      DATA A32/956.5,265.4,84.16,132.2,27*1.,0./                        
      DO 5 I=1,32                                                       
      A(I,1)=A1(I)                                                      
      A(I,2)=A2(I)                                                      
      A(I,3)=A3(I)                                                      
      A(I,4)=A4(I)                                                      
      A(I,5)=A5(I)                                                      
      A(I,6)=A6(I)                                                      
      A(I,7)=A7(I)                                                      
      A(I,8)=A8(I)                                                      
      A(I,9)=A9(I)                                                      
      A(I,10)=A10(I)                                                    
      A(I,11)=A11(I)                                                    
      A(I,12)=A12(I)                                                    
      A(I,13)=A13(I)                                                    
      A(I,14)=A14(I)                                                    
      A(I,15)=A15(I)                                                    
      A(I,16)=A16(I)                                                    
      A(I,17)=A17(I)                                                    
      A(I,18)=A18(I)                                                    
      A(I,19)=A19(I)                                                    
      A(I,20)=A20(I)                                                    
      A(I,21)=A21(I)                                                    
      A(I,22)=A22(I)                                                    
      A(I,23)=A23(I)                                                    
      A(I,24)=A24(I)                                                    
      A(I,25)=A25(I)                                                    
      A(I,26)=A26(I)                                                    
      A(I,27)=A27(I)                                                    
      A(I,28)=A28(I)                                                    
      A(I,29)=A29(I)                                                    
      A(I,30)=A30(I)                                                    
      A(I,31)=A31(I)                                                    
    5 A(I,32)=A32(I)                                                    
      IF(IOUT.EQ.0) IOUT=6                                              
!C     READ(2,501) NK                                                    
      READ(2,*) NK                                                      
      NC = NK                                                           
      DO 15 I=1,10                                                      
        DO 15 J=1,NK                                                      
            QT(I,J)=0.D0                                                      
   15 RT(I,J)=0.D0                                                      
      IF(MODEL.NE.1) GOTO 19                                            
      NG=NK                                                             
      DO 16 I=1,NK                                                      
!C  16 READ(2,502) RT(I,I),QT(I,I),(P(I,J),J=1,NK)                       
   16 READ(2,*) RT(I,I),QT(I,I),(P(I,J),J=1,NK)                         
   19 CONTINUE                                                          
      IF(MODEL.EQ.1) GOTO 21                                            
!C     READ(2,501) IOWNRQ,IOWNP                                          

!Lectura de par�metros R, Q, int
    read(2,*) IOWNRQ,IOWNP                                            

    if(IOWNRQ/=0)then !Lectura de R y Q
        do I=1,IOWNRQ                                                   
            read(2,*) K,RR(K),QQ(K)  !K= n�m de grupo                                     
        enddo
    endif            

10  continue
  
    if(IOWNP/=0)then !Lectura de par�metros de int
        do I=1,IOWNP                                                   
            read(2,*) J,K,par !A(J,K) !j y k= n�m de grupos
            j=mainsgfunc(j,ipareq)
            k=mainsgfunc(k,ipareq)
            a(j,k)=par
        enddo
    endif
                                                  
14  CONTINUE                                                          
!c------
!ccccccccccccccccccccccccccAlfonsinaccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!rngoh=0.
	!LECTURA DE LA COMPOSICION GRUPAL ASOCIATIVA
           do ja=1, nga
		read(2,*) (rngoh(i,ja),i=1,nc) 
	end do     
	


      	!LECTURA DEL NUMERO DE SITIOS Y PARAMETROS ASOCIATIVOS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC
  


	IF(NGA.GT.0) READ(2,*)(MASS(I),I=1,NGA)
    MS(:,:,:)=0                                                       
    JH(:)=0                                                           
    DO 50 I=1,NK !Lectura de compuestos
      write(*,*) "----"
       ! do j=1,10
            READ(2,*) (MS(I,J,1),MS(I,J,2),j=1,size(ms(1,:,1)))    
       ! enddo
      !  READ(2,*) NGI,(MS(I,J,1),MS(I,J,2),J=1,NGI)    
        do 50 j=1,size(ms(1,:,1))
            igamt(j,i)=ms(i,j,1)
            nytt(j,i)=ms(i,j,2)
            write(*,*)ms(i,j,1)
            write(*,*)ms(i,j,2)
50  continue	

!****Jose S****

    NUM = 0
    NMGR = 0
    NGA = 0
	DO J=1,NC
	  DO I=1,10 
	      IF(MS(J,I,1).EQ.0)CYCLE
	      IREC2=MS(J,I,1)+(ipareq-1)*150      
      	    READ(14,500,REC=IREC2)MGR,RRT,ICS,ITS1,ITS2 !500  FORMAT(2x,I2,37X,D15.8,120X,3I2)
      	    if(model==2) then
              IF (ICS.NE.0) THEN !GRUPOS ASOCIATIVOS
                NG = MS(J,I,1)
                CALL BUSCARAS (NG,NGRUPA,NMG,VL)
                IF(VL)THEN
                    IF(MS(J,I,1).EQ.10) THEN
                        RNGOH(J,NPUNTA(NG))=CANIL(J)
                    ELSE
                        RNGOH(J,NPUNTA(NG))=MS(J,I,2)
                    ENDIF                  
                ELSE
                    NGA=NGA+1
                    NPUNTA(NG) = NGA
                    NGRUPA(NGA) = NG
                    CS(NGA) = ICS
                    MASS(nga) = ICS
                    TS(NGA,1) = ITS1
                    TS(NGA,2) = ITS2
                    IF(MS(J,I,1).EQ.10) THEN
                        RNGOH(J,NGA)=CANIL(J)
                    ELSE
                        RNGOH(J,NGA)=MS(J,I,2)
                    ENDIF
                ENDIF
              ENDIF
            endif        	    
      	    IF (NPUNT(MS(J,I,1)).EQ.0) THEN !Ordena todos los grupos
                NUM = NUM + 1
                NPUNT(MS(J,I,1)) = NUM
                NGRUP(NUM) = MS(J,I,1)
               ! RR(NUM)=RRT
                CALL BUSCARAS (MGR,NMAINGR,NMG,VL)
                IF(.NOT.VL) THEN !Crea vectores NPUNTMG y NMAINGR
                    NMGR = NMGR+1
                    NPUNTMG(MGR) = NMGR
                    NMAINGR(NMGR) = MGR
                ENDIF
            END IF
        ENDDO
    ENDDO


!C....................................................................................
!C.....Genera las matrices ENASS Y RKASS con los par�metros de energ�a y volumen 
!c.....de asociaci�n, respectivamente, seg�n los grupos asociativos de los componentes
!c.....del sistema que se est� corriendo. (si el modelo elegido es A-UNIFAC)
!C....................................................................................     
    enass(:,:,:,:)=0.0
    rkass(:,:,:,:)=0.0
    if(model==2)then
      DO J=1,nga
        DO K=1,nga
            DO BB=1,CS(J)
                DO AA=1,CS(K)
                    IREC1 = MAINSGfunc(NGRUPA(K),IPAREQ)+(ipareq-1)*70
                    IF(TS(J,BB).EQ.1.OR.TS(K,AA).EQ.1)THEN !Si alguno de ambos sitios es del tipo 1
                       CALL LEEPAR (J,IREC1,IPAREQ,NGRUPA,ENASST,RKASST)
                       ENASS(AA,K,BB,J)=ENASST
                       RKASS(AA,K,BB,J)=RKASST                        
                    ELSEIF (TS(J,BB).NE.TS(K,AA)) THEN
                       CALL LEEPAR (J,IREC1,IPAREQ,NGRUPA,ENASST,RKASST)
                       ENASS(AA,K,BB,J)=ENASST
                       RKASS(AA,K,BB,J)=RKASST
                    ELSE !ELSEIF (TS(J,B).EQ.TS(K,A)) THEN
                       ENASS(AA,K,BB,J)=0.0
                       RKASS(AA,K,BB,J)=0.0     
                    ENDIF   
                ENDDO
            ENDDO
        ENDDO
      ENDDO
    endif

!****Jose S****


!!C
!!c     lectura parametros ENERG�TICOS DE asociacion
!!c
!	DO J=1,NGA
!			IF(MASS(J).EQ.0) GO TO 201
!		DO L=1,NGA
!				IF(MASS(L).EQ.0) GO TO 103
!			DO I=1,MASS(J)
!				DO K=1,MASS(L)
!						READ(2,*)ENASS(I,J,K,L) 
!				END DO
!			END DO
!  103			CONTINUE
!			END DO
!  201		   CONTINUE
!	END DO
!!C
!!c     lectura parametros VOLUM�TRICOS DE asociacion
!!C
!      DO J=1,NGA
!			IF(MASS(J).EQ.0) GO TO 3330
!		DO L=1,NGA
!				IF(MASS(L).EQ.0) GO TO 5550
!			DO I=1,MASS(J)
!				DO K=1,MASS(L)
!						READ(2,*) RKASS(I,J,K,L)
! 179						FORMAT (F16.10)
!				END DO
!			END DO
! 5550			CONTINUE 
!		END DO
! 3330			CONTINUE
!	END DO


!ccccccccccccccccccccccccccAlfonsinaccccccccccccccccccccccccccccccccccccccccccccccccccccc


      IC=1                                                              
      DO 71 I=1,NK                                                      
        DO 70 J=1,10                                                      
            IF(MS(I,J,1).EQ.0) GOTO 71                                        
            IH(IC)=MS(I,J,1)                                                  
            IF(IC.EQ.1) GOTO 69                                               
            IF(IH(IC).EQ.IH(IC-1)) GOTO 70                                    
            IF(IH(IC).GT.IH(IC-1)) GOTO 69                                    
            IF(IC.GT.2) GOTO 55                                               
            IHH=IH(1)                                                         
            IH(1)=IH(2)                                                       
            IH(2)=IHH                                                         
            GOTO 69                                                           
   55       I1=IC-1                                                           
            DO 65 I2=1,I1                                                     
                IF(IH(IC).GT.IH(I2)) GOTO 65                                      
                IF(IH(IC).EQ.IH(I2)) GOTO 70                                      
                I4=IC-I2                                                          
                DO 61 I3=1,I4                                                     
   61               IH(IC+1-I3)=IH(IC-I3)                                             
                IH(I2)=MS(I,J,1)                                                  
   65       CONTINUE                                                          
   69       IC=IC+1                                                           
            IF(IC.GT.20) WRITE(6,607)                                         
            IF(IOUT.NE.6.AND.IC.GT.20) WRITE(IOUT,607)                        
   70   CONTINUE                                                          
   71 CONTINUE                                                          
      IC=IC-1                                                           
!c------
      nktt=ic
!c------
      DO 73 I=1,IC                                                      
   73   JH(IH(I))=I                                                       
      DO 72 I=1,10                                                      
        DO 72 J=1,20                                                      
   72       NY(I,J)=0                                                         
      DO 75 I=1,NK                                                      
        DO 74 J=1,10                                                      
            IF(MS(I,J,1).EQ.0) GOTO 75                                        
            N1=MS(I,J,1)                                                      
            N2=MS(I,J,2)                                                      
            IF(N1.EQ.0) GOTO 75                                               
            N3=JH(N1)                                                         
   74   NY(I,N3)=N2                                                       
   75 CONTINUE                                                          
      I=0                                                               
      NGMGL=0                                                           
      DO 80 K=1,IC                                                      
        NSG=IH(K)                                                         
        NGMNY=MAINSGfunc(NSG,ipareq)                                                 
        IF(NGMNY.NE.NGMGL) I=I+1                                          
        NGM(I)=NGMNY                                                      
        NGMGL=NGMNY                                                       
      DO 80 J=1,NK                                                      
      RT(I,J)=RT(I,J)+NY(J,K)*RR(NSG)                                   
   80 QT(I,J)=QT(I,J)+NY(J,K)*QQ(NSG)                                   
      NG=I                                                              
      !WRITE(6,608) (IH(K),K=1,IC)                                       
      !WRITE(6,609) (MAINSGfunc(IH(K),ipareq),K=1,IC)                               
      !WRITE(6,610)                                                      
      !DO 90 I=1,NK                                                      
   !90 WRITE(6,611) I,(NY(I,K),K=1,IC)                                   
      !WRITE(6,699)                                                      
      IF(IOUT.EQ.6) GOTO 85                                             
      !WRITE(IOUT,608) (IH(K),K=1,IC)                                    
      !WRITE(IOUT,609) (MAINSGfunc(IH(K),ipareq),K=1,IC)                            
      !WRITE(IOUT,610)                                                   
      !DO 91 I=1,NK                                                      
   !91 WRITE(IOUT,611) I,(NY(I,K),K=1,IC)                                
      !WRITE(IOUT,699)                                                   
   85 CONTINUE                                                          
      DO 20 I=1,NG                                                      
      DO 20 J=1,NG                                                      
      NI=NGM(I)                                                         
      NJ=NGM(J)                                                         
   20 P(I,J)=A(NI,NJ)                                                   
      !WRITE(6,612)                                                      
      DO 95 K=1,IC                                                      
   95 NN=IH(K)                                                          
   !95 WRITE(6,613) NN,RR(NN),QQ(NN)                                     
      !WRITE(6,699)                                                      
      IF(IOUT.EQ.6) GOTO 99                                             
      WRITE(IOUT,612)                                                   
      DO 96 K=1,IC                                                      
      NN=IH(K)                                                          
   96 WRITE(IOUT,613) NN,RR(NN),QQ(NN)                                  
      WRITE(IOUT,699)                                                   
   99 CONTINUE                                                          
   21 CONTINUE                                                          
      !WRITE(6,604)   
                                                    
      !DO 25 I=1,NG                                                      
   !25 WRITE(6,603) (P(I,J),J=1,NG)                                      
      !WRITE(6,699)                                                      
      !IF(MODEL.EQ.0) WRITE(6,605)                                       
      !IF(MODEL.EQ.1) WRITE(6,627)                                       
      IF(IOUT.EQ.6) GOTO 26                                             
      WRITE(IOUT,604)                                                   
      DO 27 I=1,NG                                                      
   27 WRITE(IOUT,603) (P(I,J),J=1,NG)                                   


!ccccccccccccccccc Escritura de los par�metros de asociaci�n ALFONSINAccccccccccccccccc
       IF (NGA.GT.0) THEN

	WRITE(1,218) (I, I=1,NC)
 218	FORMAT(/,X,'"ASSOC GROUP COMPOSITION" ',/,23X,'COMPONENTES',/,' GRUPO 	  #  SIT ASOC  ',2I5,/) 

	DO ja=1,NGA
	WRITE(1,219) ja,MASS(ja)  , (rngoh(i,ja),i=1,nc) 
 219	FORMAT(3X,I3,9X,I3,6X,2f5.1)
	END DO

	WRITE(1,220)
 220	FORMAT(/,X,'PARAMETROS DE ENERGIA DE ASOCIACION (Kelvin)  ',/)  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
	DO J=1,NGA
			IF(MASS(J).EQ.0) GO TO 202
		DO L=1,NGA
				IF(MASS(L).EQ.0) GO TO 104
			DO I=1,MASS(J)
				DO K=1,MASS(L)
						WRITE(1,221) I,J,K,L,ENASS(I,J,K,L)
  221	FORMAT(X,' ENASS( ',I3,I3,I3,I3,' ) = ',F10.4)
				END DO
			END DO
  104			CONTINUE
			END DO
  202			CONTINUE
	END DO


	WRITE(1,222)
  222	FORMAT (/,X,'PARAMETROS DE VOLUMEN DE ASOCIACI�N (cm3/mol) ',/)

	DO J=1,NGA
			IF(MASS(J).EQ.0) GO TO 301
		DO L=1,NGA
				IF(MASS(L).EQ.0) GO TO 5011
			DO I=1,MASS(J)
				DO K=1,MASS(L)
					WRITE(1,223) I,J,K,L,RKASS(I,J,K,L)
  223	FORMAT(X,' RKASS( ',I3,I3,I3,I3,' ) = ',F10.4)
				END DO
			END DO
 5011			CONTINUE
			END DO
  301		   CONTINUE
	END DO
	END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC
!c------
      WRITE(IOUT,699)                                                   
      IF(MODEL.EQ.0) WRITE(IOUT,605)                                    
      IF(MODEL.EQ.1) WRITE(IOUT,627)                                    
   26 CONTINUE                                                          
      DO 30 I=1,NK                                                      
      Q(I)=0.D0                                                         
      R(I)=0.D0                                                         
      DO 30 K=1,NG                                                      
      Q(I)=Q(I)+QT(K,I)                                                 
   30 R(I)=R(I)+RT(K,I)                                                 
      !DO 40 I=1,NK                                                      
   !40 WRITE(6,606) I,R(I),Q(I)                                          
      IF(IOUT.EQ.6) GOTO 42                                             
      DO 41 I=1,NK                                                      
   41 WRITE(IOUT,606) I,R(I),Q(I)                                       
   42 CONTINUE      
  500 FORMAT(2x,I2,37X,D15.8,120X,3I2)                                                     
  501 FORMAT(20I3)                                                      
  502 FORMAT(8F10.2)                                                    
  503 FORMAT(I3,2F10.2)                                                 
  504 FORMAT(2I3,F10.2)                                                 
  603 FORMAT(1X,10F12.3)                                                
  604 FORMAT('  INTERACTION PARAMETERS',/)                              
  605 FORMAT(' UNIFAC MOLECULAR R AND Q',/)                             
  606 FORMAT(I5,2F15.4)                                                 
  607 FORMAT(' ** WARNING: NUMBER OF SUB GROUPS MUST NOT EXCEED 20 **') 
  608 FORMAT(//,' SUB GROUPS :',20I3)                                   
  609 FORMAT(' MAIN GROUPS:',20I3)                                      
  610 FORMAT(' COMPONENT')                                              
  611 FORMAT(6X,I2,5X,20I3)                                             
  612 FORMAT(' GROUP R- AND Q-VALUES',/)                                
  613 FORMAT(1X,I3,2F10.4)                                              
  627 FORMAT(' SPECIFIED UNIQUAC R AND Q',/)                            
  699 FORMAT(//)                                                        
!c------
 1603 format('  ASSOCIATION PARAMETERS',//,10X,'K(OH)   :',F7.3,/,10X,'E(OH)/k :',F7.1,' K-1')
!c------
      RETURN                                                            
    endsubroutine PARIN2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
      SUBROUTINE PARAM2                                                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CUFAC/NK,NG,P(10,10),T                                     
      COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
      COMMON/CQT/QT(10,10),Q(10),R(10)                                  
      DO 30 I=1,NG                                                      
        DO 30 J=1,NG                                                      
   30       TAU(I,J)=DEXP(-P(I,J)/T)                                          
   40 CONTINUE                                                          
      DO 50 I=1,NK                                                      
        DO 50 K=1,NG                                                      
            S(K,I)=0.D0                                                       
      DO 50 M=1,NG                                                      
   50   S(K,I)=S(K,I)+QT(M,I)*TAU(M,K)                                    
      DO 60 I=1,NK                                                      
        F(I)=1.D0                                                         
      DO 60 J=1,NG                                                      
   60   F(I)=F(I)+QT(J,I)*DLOG(S(J,I))                                    
      RETURN                                                            
      END                                                               
      SUBROUTINE GAMINF(N1,N2,GAM)                                      
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CUFAC/NK,NG,P(10,10),T                                     
      COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
      COMMON/CQT/QT(10,10),Q(10),R(10)                                  
!c------
      common/asoc/nktt,igamt(20,12),nytt(20,12)  
      common/nga/nga,mass(12)
      common/grupas1/rkass(6,12,6,12),enass(6,12,6,12), deloh(6,12,6,12) !Alfonsina
!c------
      dimension xohi0(10),xnohi0(10),tgt(10),xgamk(20),x(2)
	common/ioh2/rngoh(12,12)
	common/ioh2sis/rngoht(3,2)

!c------
      dk=1.381e-23
      deloh=0.0
      xnoh=0.0
      xnoh1=0.0
	xoh=0.0
      xgam=0.0
      xgamt=0.0
      do 7777 i=1,10
      xnohi0(i)=0.0
      tgt(i)=0.0
 7777 continue
      do 8888 k=1,20
      xgamk(k)=0.0
 8888 continue
!c------
!c      x(n1)=1.0
!c      x(n2)=0.0
!c------


!c      do 33 jc=1,2
!c      if(jc.eq.1) then
!c      i=n1
!c      else if(jc.eq.2) then
!c      i=n2
!c     end if
!c      goh(i)=0
!c      tgt(i)=0.0
!c      xnohi0(i)=0.0
!c      xgam=0.0
!c      do 33 j=1,nktt
!c      if(((igamt(j,i).eq.14).or.(igamt(j,i).eq.17)).and.(nytt(j,i). 
!c     *ne.0)) goh(i)=nytt(j,i)
!c   33 continue
!c      do 32 jc=1,2
!c      if(jc.eq.1) then
!c      i=n1
!c      else if(jc.eq.2) then
!c      i=n2
!c      end if
!c      if(nga.eq.1) then
!c      do 5 k=1,nktt
!c      tgt(i)=tgt(i)+nytt(k,i)
!c   5 continue
      !alterar
!c	xnoh1=xnoh1+goh(i)*x(i)
	!write (1,*) "xnoh1parte2=", xnoh1 
!c      end if
!c      xnohi0(i)=goh(i)/R(i)
	!write (1,*) "xnohi0parte2(",i,")=", xnohi0(i) 
!c   32 continue
!c      do 7 jc=1,2
!c      if(jc.eq.1) then
!c      i=n1
!c      else if(jc.eq.2) then
!c      i=n2
!c      end if
!c	!alterar
!c      xgam=xgam+R(i)*x(i)
	!write (1,*) "xgamparte2=", xgam
!c    7 continue
      !alterar
!c     xnoh=xnoh1/xgam
	!write (1,*) "xnohparte2=", xnoh 
!c------
      Q1=Q(N2)/Q(N1)                                                    
      R1=R(N2)/R(N1)                                                    
      QR=R1/Q1                                                          
      GAM=F(N2)+Q(N2)*(1.-DLOG(Q1))-R1+DLOG(R1)-5.D0*Q(N2)*(1.-QR+DLOG(R1)-DLOG(Q1))                                                      
      DO 10 I=1,NG                                                      
   10 GAM=GAM-S(I,N2)/S(I,N1)*QT(I,N1)-QT(I,N2)*DLOG(S(I,N1))           
!c------
!c      deloh=p4*(dexp(p5/t)-1.)
!c      if((xnohi0(n2).eq.0.0).or.(deloh.eq.0.0)) then
!c      xohi0(n2)=1.0
!c      else
!c      xohi0(n2)=(-1.+dsqrt(1.+4.*xnohi0(n2)*deloh))/(2.*xnohi0(n2)*
!c     *deloh)
!c      end if
!c      if((xnoh.eq.0.0).or.(deloh.eq.0.0)) then
!c      xoh=1.0
!c      else
!c      xoh=(-1.+dsqrt(1.+4.*xnoh*deloh))/(2.*xnoh*deloh)
!c      end if
!	!alterar
!c      gam=gam+goh(n2)*(2.*dlog(xoh/xohi0(n2))+xohi0(n2)-xoh)-(2.*xoh-
!c     *xoh**2.)*deloh*xnoh*(goh(n2)-R(n2)*xnoh)/(1.+2.*xnoh*xoh*deloh)
!c------
      RETURN                                                            
      END                                                               
      SUBROUTINE STIG(Y,S)                                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
      COMMON/CUFAC/N,NG,P(10,10),T                                      
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AA(10),DA(10,10),XM(10,4)                                       
      DIMENSION Y(10),V(10),YGEM(10)                                    
      common/nga/nga,mass(12)
      JPUR=0                                                            
      AMAX=0.                                                           
      DO 10 I=1,N                                                       
      IF(A(I).LT.AMAX) GOTO 10                                          
      JPUR=I                                                            
      AMAX=A(I)                                                         
   10 CONTINUE                                                          
      RMAX=1.D5                                                         
      NGN=N                                                             
      IF(NF.GT.1) NGN=N+NF                                              
      NEG=0                                                             
      DO 100 KK=1,NGN                                                   
      JM=KK                                                             
      IF(JPUR.NE.0) JM=JPUR                                             
      IF(JM.LE.N) GOTO 30                                               
      DO 20 I=1,N                                                       
   20 Y(I)=Z(I)*(2+XVL(I,JM-N)/SFAS(JM-N))/3                            
      GOTO 40                                                           
   30 SUM=0.                                                            
      DO 35 I=1,N                                                       
      GG=A(I)-GAM(JM,I)                                                 
      IF(GG.LT.-50.D0) GG=-50.D0                                        
      Y(I)=DEXP(GG)                                                     
   35 SUM=SUM+Y(I)                                                      
   40 NA=3                                                              
      DO 43 K=1,NA                                                      
      DO 36 I=1,N                                                       
   36 Y(I)=Y(I)/SUM                                                     
      CALL unifac(1,Y,AA,DA,PACT)                                       
      IF(K.EQ.NA) GOTO 44                                               
      DO 41 I=1,N                                                       
   41 Y(I)=DEXP(A(I)-AA(I))                                             
   42 SUM=0.                                                            
      DO 43 I=1,N                                                       
   43 SUM=SUM+Y(I)                                                      
   44 CONTINUE                                                          
      YV1=0.                                                            
      DO 50 J=1,NF                                                      
   50 V(J)=0.                                                           
      DO 60 I=1,N                                                       
      GD=DLOG(Y(I))+AA(I)-A(I)                                          
      YV1=YV1+Y(I)*GD                                                   
      DO 60 J=1,NF                                                      
      K=J                                                               
      VV=XVL(I,K)*Z(I)/SFAS(K)                                          
      D=GD*(Y(I)-VV)                                                    
   60 V(J)=V(J)+D                                                       
      YV2=V(1)                                                          
      DO 70 J=1,NF                                                      
      IF(V(J).LT.YV2) YV2=V(J)                                          
   70 CONTINUE                                                          
      RT1=YV1                                                           
      IF(YV2.GT.0.) RT1=RT1-YV2/2                                       
      IF(NEG.EQ.0.AND.YV1.GT.0.) GOTO 80                                
      RT1=YV1                                                           
      IF(NEG.EQ.0) RMAX=0.                                              
      NEG=1                                                             
   80 IF(RT1.GT.RMAX) GOTO 100                                          
      S=YV1                                                             
      RMAX=RT1                                                          
      CC=DEXP(-YV1)                                                     
      DO 90 I=1,N                                                       
   90 YGEM(I)=Y(I)*CC                                                   
      IF(JPUR.NE.0) GOTO 110                                            
  100 CONTINUE                                                          
  110 DO 120 I=1,N                                                      
  120 Y(I)=YGEM(I)                                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE GMIX(NARG,NDIF,FUN,GRAD,XMAT,YVAL)                     
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),XX(10),DA(10,10),XM(10,4)                                       
      DIMENSION YVAL(30),GRAD(30),X(10),TM(10),FG(10),XMAT(30,30)       
      NG=NF-1                                                           
      N=NARG/NG                                                         
      JD=1                                                              
      IF(NDIF.EQ.2) JD=2                                                
      IF(NDIF.EQ.2) CALL CHECK(N,YVAL)                                  
      IF(NF.NE.NG) GOTO 20                                              
      NARG=NARG-N                                                       
      NG=NG-1                                                           
      DO 10 I=1,N                                                       
      DO 10 J=1,NG                                                      
   10 YVAL(I+(J-1)*N)=XVL(I,J)                                          
   20 FUN=-GNUL                                                         
      DO 50 I=1,N                                                       
      XVL(I,NF)=1.                                                      
      DO 30 J=1,NG                                                      
      XVL(I,J)=YVAL(I+(J-1)*N)                                          
   30 XVL(I,NF)=XVL(I,NF)-XVL(I,J)                                      
      ! write(*,*)xvl
      DO 40 J=1,NF                                                      
      IF(XVL(I,J).GT.0.) GOTO 40                                        
      FUN=0.                                                            
      GOTO 1000                                                         
   40 CONTINUE                                                          
   50 CONTINUE                                                          
      DO 200 J=1,NF                                                     
      SFAS(J)=0.                                                        
      DO 60 I=1,N                                                       
      X(I)=XVL(I,J)*Z(I)                                                
   60 SFAS(J)=SFAS(J)+X(I)                                              
      DO 65 I=1,N                                                       
      XX(I)=X(I)/SFAS(J)                                                
   65 XM(I,J)=XX(I)                                                     
      CALL unifac(JD,XX,FG,DA,PACT)                                     
      IDUM(J)=NDUM                                                      
      DO 70 I=1,N                                                       
      TM(I)=DLOG(XVL(I,J)/SFAS(J))+FG(I)                                
   70 FUN=FUN+X(I)*TM(I)                                                
      IF(NDIF.EQ.0) GOTO 200                                            
      DO 80 I=1,N                                                       
      S=Z(I)*TM(I)                                                      
      IF(J.EQ.NF) GOTO 75                                               
      GRAD(I+(J-1)*N)=S                                                 
      GOTO 80                                                           
   75 DO 76 K=1,NG                                                      
      NK=I+(K-1)*N                                                      
   76 GRAD(NK)=GRAD(NK)-S                                               
   80 CONTINUE                                                          
      IF(NDIF.EQ.1) GOTO 200                                            
      DO 100 I=1,N                                                      
      ST=Z(I)/SFAS(J)                                                   
      DO 100 L=1,N                                                      
      S=ST*(DA(I,L)-1.)*Z(L)                                            
      IF(L.EQ.I)S=S+Z(I)/XVL(I,J)                                       
      IF(J.EQ.NF) GOTO 90                                               
      XMAT(I+(J-1)*N,L+(J-1)*N)=S                                       
      GOTO 95                                                           
   90 DO 92 K=1,NG                                                      
      DO 92 M=1,K                                                       
      NK=I+(K-1)*N                                                      
      NM=L+(M-1)*N                                                      
      IF(K.NE.M) XMAT(NK,NM)=S                                          
      IF(K.EQ.M) XMAT(NK,NM)=XMAT(NK,NM)+S                              
   92 CONTINUE                                                          
   95 CONTINUE                                                          
  100 CONTINUE                                                          
  200 CONTINUE                                                          
 1000 RETURN                                                            
      END                                                               
      SUBROUTINE STABIL(N,NDIF,FUN,GRAD,XMAT,Y)                         
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
      DIMENSION GRAD(30),XMAT(30,30),Y(30),YEX(10),XEX(10)              
      common/nga/nga,mass(12)

      SUM=0.                                                            
      DO 10 I=1,N                                                       
      YEX(I)=0.                                                         
      IF(Y(I).GT.-40.) YEX(I)=DEXP(Y(I))                                
   10 SUM=SUM+YEX(I)                                                    
      DO 15 I=1,N                                                       
   15 XEX(I)=YEX(I)/SUM                                                 
      JD=1                                                              
      IF(NDIF.EQ.2) JD=2                                                
      CALL unifac(JD,XEX,AL,DA,PACT)                                    
      FUN=1.                                                            
      DO 20 I=1,N                                                       
      S=Y(I)+AL(I)-A(I)                                                 
      IF(NDIF.EQ.0) GOTO 20                                             
      GRAD(I)=YEX(I)*S                                                  
   20 FUN=FUN+YEX(I)*(S-1)                                              
      IF(NDIF.LT.2) GOTO 50                                             
      DO 30 I=1,N                                                       
      S=XEX(I)                                                          
      DO 40 J=1,I                                                       
      XMAT(I,J)=S*YEX(J)*DA(I,J)                                        
   40 XMAT(J,I)=XMAT(I,J)                                               
   30 XMAT(I,I)=XMAT(I,I)+YEX(I)+GRAD(I)                                

  
  
   50 RETURN                                                            
      END                                                               
      SUBROUTINE CHECK(N,YVAL)                                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
      COMMON/CIPR/IPR                                                   
      COMMON/COUT/IOUT                                                  
      DIMENSION YVAL(30)                                                
      JMAX=NF                                                           
      SMAX=SFAS(NF)                                                     
      NT=NF-1                                                           
      DO 5 J=1,NT                                                       
      IF(SFAS(J).LT.SMAX) GOTO 5                                        
      SMAX=SFAS(J)                                                      
      JMAX=J                                                            
    5 CONTINUE                                                          
      IF(JMAX.EQ.NF) GOTO 20                                            
      DO 10 I=1,N                                                       
      DS=1.                                                             
      DO 15 J=1,NT                                                      
   15 DS=DS-YVAL(I+(J-1)*N)                                             
   10 YVAL(I+(JMAX-1)*N)=DS                                             
   20 IF(NF.EQ.2) GOTO 100                                              
      DO 21 I=1,N                                                       
      XVL(I,NF)=1.                                                      
      DO 21 J=1,NT                                                      
      XVL(I,J)=YVAL(I+(J-1)*N)                                          
   21 XVL(I,NF)=XVL(I,NF)-XVL(I,J)                                      
      DO 23 J=1,NF                                                      
      SFAS(J)=0.                                                        
      DO 22 I=1,N                                                       
      AL(I)=XVL(I,J)*Z(I)                                               
   22 SFAS(J)=SFAS(J)+AL(I)                                             
      DO 23 I=1,N                                                       
   23 XM(I,J)=AL(I)/SFAS(J)                                             
      DO 30 I=1,NT                                                      
      JN=I+1                                                            
      DO 30 J=JN,NF                                                     
      IF(DABS(XM(MAXZ,I)-XM(MAXZ,J)).LT.5.D-3) GOTO 40                  
   30 CONTINUE                                                          
      GOTO 100                                                          
   40 IV=I                                                              
      JV=J                                                              
      DMAX=0.                                                           
      DO 50 I=1,N                                                       
      R1=DABS(XM(I,IV)-XM(I,JV))                                        
      IF(R1.GT.DMAX) DMAX=R1                                            
   50 CONTINUE                                                          
      IF(DMAX.GT.2.5D-2) GOTO 100                                       
      WRITE(6,601)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,601)                                     
      NF=NF-1                                                           
      NT=NT-1                                                           
      DO 60 I=1,N                                                       
      XVL(I,IV)=XVL(I,IV)+XVL(I,JV)                                     
   60 XVL(I,JV)=XVL(I,NF+1)                                             
  100 CONTINUE                                                          
      IF(NF.LT.3) GOTO 250                                              
      MINF=0                                                            
      DO 200 I=1,NF                                                     
      IF(SFAS(I).LT.1.D-12) MINF=I                                      
  200 CONTINUE                                                          
      IF(MINF.EQ.0) GOTO 250                                            
      WRITE(6,602) MINF,SFAS(MINF)                                      
      IF(IOUT.NE.6) WRITE(IOUT,602) MINF,SFAS(MINF)                     
      DO 220 I=1,NF                                                     
      IF(I.EQ.MINF) GOTO 220                                            
      DO 210 J=1,N                                                      
  210 XVL(J,I)=XVL(J,I)+SFAS(I)*XVL(J,MINF)                             
  220 CONTINUE                                                          
      IF(MINF.EQ.NF) GOTO 250                                           
      NNF=NF-1                                                          
      DO 230 I=1,NNF                                                    
      IF(I.LT.MINF) GOTO 230                                            
      DO 240 J=1,N                                                      
  240 XVL(J,I)=XVL(J,I+1)                                               
  230 CONTINUE                                                          
      NF=NF-1                                                           
  250 CONTINUE                                                          
  601 FORMAT(//,' * * * TWO PHASES ARE IDENTICAL * * *'//)              
  602 FORMAT(/,' * THE AMOUNT OF PHASE',I2,' IS',D12.4,'. THE PHASE IS ELIMINATED *',/)                                                   
      RETURN                                                            
      END                                                               
      SUBROUTINE SPLIT(ND,N,IDI,BETA,DEL,G,W)                           
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION G(ND,ND),W(ND,5)                                        
      IDI=0                                                             
      DO 20 I=1,N                                                       
      I1=I-1                                                            
      I2=I+1                                                            
      ID=0                                                              
      TV=G(I,I)                                                         
      SV=TV                                                             
      IF (I1.EQ.0) GO TO 25                                             
      DO 30 J=1,I1                                                      
   30 SV=SV-G(I,J)**2                                                   
   25 IF (SV.LT. DEL*DEL) ID=1                                          
      SVR=DSQRT(DABS(SV))                                               
      IF (SVR.LT. DEL) SVR=DEL                                          
      XM=0.                                                             
      IF (I.EQ.N) GO TO 35                                              
      DO 40 J=I2,N                                                      
      S=G(J,I)                                                          
      IF (I1.EQ.0) GO TO 45                                             
      DO 50 K=1,I1                                                      
   50 S=S-G(I,K)*G(J,K)                                                 
   45 S=S/SVR                                                           
      IF (DABS(S).GT. XM) XM=DABS(S)                                    
   40 G(J,I)=S                                                          
   35 IF (XM.LT. BETA) GO TO 55                                         
      ID=1                                                              
      XM=XM/BETA                                                        
      DO 60 J=I,N                                                       
   60 G(J,I)=G(J,I)/XM                                                  
      SVR=SVR*XM                                                        
   55 IF (ID.EQ.1) W(I,1)=SVR**2-SV                                     
      G(I,I)=SVR                                                        
      IDI=IDI+ID                                                        
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE LINE(ND,N,XLAM,GD,G,W)                                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION GD(ND),G(ND,ND),W(ND,5)                                 
   65 W(1,2)=-GD(1)/(G(1,1)+XLAM)                                       
      IF (N.EQ.1) GO TO 75                                              
      DO 70 I=2,N                                                       
      I1=I-1                                                            
      S=-GD(I)                                                          
      DO 80 J=1,I1                                                      
   80 S=S-G(I,J)*W(J,2)                                                 
   70 W(I,2)=S/(G(I,I)+XLAM)                                            
   75 W(N,3)=W(N,2)/(G(N,N)+XLAM)                                       
      IF (N.EQ.1) GO TO 85                                              
      DO 90 II=2,N                                                      
      I=N+1-II                                                          
      I1=I+1                                                            
      S=W(I,2)                                                          
      DO 100  J=I1,N                                                    
  100 S=S-G(J,I)*W(J,3)                                                 
   90 W(I,3)=S/(G(I,I)+XLAM)                                            
   85 RETURN                                                            
      END                                                               
      SUBROUTINE TMSSJ(ND,N,IPR,NMAX,XLAM,GLIM,F   ,X,GD,G,W,ifunc)     
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/COUT/IOUT                                                  
      integer::ifunc
      DIMENSION X(ND),GD(ND),G(ND,ND),W(ND,5)                           
      NTM=0                                                             
      DEL=1.D-7                                                         
  150 NTM=NTM+1                                                         
      GNM=0.                                                            
      
      if(ifunc==1)then
        CALL STABIL(N,2,F,GD,G,X)                                       
      elseif(ifunc==2)then
        CALL GMIX(N,2,F,GD,G,X)  
      endif
      DO 5 I=1,N                                                        
    5 GNM=GNM+GD(I)**2                                                  
      IF (GNM.LT. GLIM) GO TO 200                                       
      BETA=0.                                                           
      DO 10 I=1,N                                                       
      W(I,4)=X(I)                                                       
      W(I,5)=GD(I)                                                      
      DO 10 J=1,I                                                       
      S=DABS(G(I,J))                                                    
      IF (I.EQ.J) S=S*N                                                 
      IF (S.GT. BETA) BETA=S                                            
   10 W(I,1)=0.                                                         
      BETA=DSQRT(BETA/N)                                                
      CALL SPLIT(ND,N,IDI,BETA,DEL,G,W)                                 
      XLM=0.                                                            
      NTS=0                                                             
  350 NTS=NTS+1                                                         
      CALL LINE(ND,N,XLM,GD,G,W)                                        
      SMAX=0.                                                           
      GP=0.                                                             
      DP=0.                                                             
      DO 205 I=1,N                                                      
      S=W(I,3)                                                          
      IF (DABS(S).GT. SMAX) SMAX=DABS(S)                                
      GP=GP+S*GD(I)                                                     
  205 DP=DP+S*S*W(I,1)                                                  
      FAC=1.                                                            
      IF (SMAX.GT. XLAM) FAC=XLAM/SMAX                                  
      DER=FAC*GP                                                        
      ALFA=1.  
      cont = 0
  210 FF=ALFA*FAC                                                       
      DO 215 I=1,N                                                      
  215 X(I)=W(I,4)+FF*W(I,3)                                             
 
       if(ifunc==1)then
        CALL STABIL(N,1,FNEW,GD,G,X)                           
      elseif(ifunc==2)then
        CALL GMIX(N,1,FNEW,GD,G,X)   
      endif
 
      cont = cont + 1
      IF(cont > 10000)THEN 
      GOTO 230
      ENDIF
 
                                             
      IF(FNEW.NE.0.) GOTO 220                                           
      ALFA=ALFA/2.                                                      
      GOTO 210                                                          
  220 DELS=FNEW-F                                                       
       IF (FNEW.NE. 0.  .AND. DELS.LT. 1.D-10) GO TO 230                
      IF (NTS.GT. 1) GO TO 125                                          
      DE2=(DELS-ALFA*DER)/ALFA**2/2                                     
      GT=-DER/DE2                                                       
      IF (GT.GT. ALFA) ALFA=ALFA/3.                                     
      IF (GT.LE. ALFA) ALFA=GT                                          
      GO TO 210                                                         
  230 PRED=FF*GP-FF**2/2*(GP+DP)                                        
      CALPRE=DELS/PRED                                                  
      F=FNEW                                                            
      DO 345 I=1,N                                                      
  345 W(I,4)=X(I)                                                       
      IF (NTS.GE.3) GO TO 125                                           
      IF (ALFA.EQ. 1.  .AND. CALPRE.GT. .8) GO TO 350                   
  125 DO 126 I=1,N                                                      
  126 X(I)=W(I,4)                                                       
      IF (IPR.NE.0) WRITE(6,130) NTM,F,CALPRE,GNM,ALFA                  
      IF(IOUT.NE.6.AND.IPR.NE.0) WRITE(IOUT,130) NTM,F,CALPRE,GNM,ALFA  
  130 FORMAT(1X,I2,':  F = ',D16.8,'  CALPRE = ',F8.3,'  GRAD = ',D16.8,'  ALFA = ',D10.2)                                                
      IF (IPR.GT.1) WRITE(6,131) (X(I),I=1,N)                           
      IF(IOUT.NE.6.AND.IPR.GT.1) WRITE(IOUT,131) (X(I),I=1,N)           
  131 FORMAT(' X-VECTOR',10F11.5,(/,9X,10F11.5))                        
      IF (IPR.GT. 1) WRITE(6,132)                                       
      IF(IOUT.NE.6.AND.IPR.GT.1) WRITE(IOUT,132)                        
  132 FORMAT(' ')                                                       
      IF (NTM.LT. NMAX) GO TO 150                                       
  200 IF(IPR.GT.0) WRITE(6,133) NTM,GNM,F                               
      IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,133) NTM,GNM,F              
  133 FORMAT(/,'  NUMBER OF ITERATIONS = ',I2,', NORM OF GRADIENT = ',  D12.5,', F = ',D14.6)                                             
      IF(IPR.GT.0) WRITE(6,131) (X(I),I=1,N)                            
      IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,131) (X(I),I=1,N)           
      RETURN                                                            
      END                                                               
      SUBROUTINE BINOD                                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      COMMON/CY/Y13,Y21,STEP                                            
      COMMON/COUT/IOUT                                                  
      DIMENSION YMAT(70,6),NITE3(70)                                    
      COMMON/CBISO/NIC1,NIC2,IC1(120),IC2(120)                          
      DIMENSION Y(4),YOLD(4),DY(4),DYOLD(4),DMAT(4,4)                   
      DIMENSION K(6)                                                    
!c-----
      dimension xandrea1(10),xandrea2(10),act1(10),act2(10),dact(10,10),pact(2,2),acttot(70,6)
      common/ig/ig
!c-----      
      DATA K/1,3,2,4,5,6/                                               
      WRITE(6,601)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,601)                                     
  601 FORMAT(///,10X,' ** BINODAL CURVE, CONCENTRATIONS IN MOLE PERCENT **',//,18X,'COMP. 1',13X,'COMP. 2',13X,'COMP. 3')                 
!c-----
      write(6,1601)
      if(iout.ne.6) write(iout,1601)
 1601 format(18x,'gamma1',14x,'gamma2',14x,'gamma3',/)
!c-----     
      DO 1 I=1,4                                                        
      DO 1 J=1,4                                                        
1     DMAT(I,J)=0.D0                                                    
      IRUND=200                                                         
      NIC1=0                                                            
      NIC2=0                                                            
      ICOND=0                                                           
      NOLD=2                                                            
      N=0                                                               
      Y(1)=1.D0-Y13/100.D0                                              
      Y(2)=0.D0                                                         
      Y(3)=Y21/100.D0                                                   
      Y(4)=0.D0                                                         
   12 CALL SOLVE(Y,DY,NOLD,NEW,NITER,N,NT)                              
      IF(NITER.LE.NT) GOTO 16                                           
      IF(N.GT.0) GOTO 19                                                
      ICOND=-10                                                         
      WRITE(6,902)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,902)                                     
  902 FORMAT(/,' THE BASE LINE CALCULATION DID NOT CONVERGE IN 10 ITERATIONS'/)                                                           
      GOTO 3                                                            
   19 IF(IHALF.LT.5) GOTO 20                                            
      ICOND=-3                                                          
      GOTO 3                                                            
   20 IHALF=IHALF+1                                                     
      ST=ST/2.D0                                                        
      GOTO 17                                                           
   16 IF(DABS(Y(1)-Y(3))+DABS(Y(2)-Y(4)).GT.1.D-8) GOTO 21              
      IF(N.GT.0) GOTO 19                                                
      WRITE(6,903)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,903)                                     
  903 FORMAT(/,' THE CONCENTRATIONS ON THE BASE LINE ARE IDENTICAL IN THE TWO PHASES'/)                                                   
      GOTO 3                                                            
   21 N=N+1                                                             
      NITE3(N)=NITER                                                    
      IHALF=0                                                           
      DO 2 I=1,4                                                        
2     YMAT(N,I)=Y(I)                                                    
      IF(ICOND.EQ.2.AND.Y(1).LT.1.D-10) GOTO 3                          
      DYMAX=DABS(DY(NEW))                                               
      DO 4 I=1,4                                                        
4     DY(I)=DY(I)/DYMAX                                                 
      IF(N.EQ.1)GOTO 5                                                  
      STAP=DABS(Y(NEW)-YOLD(NEW))                                       
      IF(DY(NEW)*DYOLD(NEW).GT.0.D0)GOTO 6                              
      DO 7 I=1,4                                                        
7     DY(I)=-DY(I)                                                      
6     IF(NEW.EQ.NOLD)GOTO 8                                             
      RR=DY(NEW)/DYOLD(NEW)                                             
      DO 9 I=1,4                                                        
9     DYOLD(I)=DYOLD(I)*RR                                              
8     DO 10 I=1,4                                                       
      Z=(YOLD(I)-Y(I))/STAP                                             
      DMAT(I,3)=(3.D0*Z+2.D0*DY(I)+DYOLD(I))/STAP                       
10    DMAT(I,4)=(2.D0*Z+DY(I)+DYOLD(I))/STAP**2                         
5     ST=RUND(Y(NEW),DY(NEW),STEP,IRUND)                                
      DO 18 I=1,4                                                       
      DMAT(I,1)=Y(I)                                                    
      DMAT(I,2)=DY(I)                                                   
      YOLD(I)=Y(I)                                                      
18    DYOLD(I)=DY(I)                                                    
17    DO 11 I=1,4                                                       
      Y(I)=DMAT(I,4)                                                    
      DO 11 J=1,3                                                       
11    Y(I)=ST*Y(I)+DMAT(I,4-J)                                          
      IF(IHALF.GT.0)GOTO 12                                             
      CALL TERM(Y,DMAT,ICOND,NEW)                                       
      NOLD=NEW                                                          
      IF(ICOND.EQ.0.OR.ICOND.EQ.2) GOTO 12                              
3     NGIT=N                                                            
      IF(N.EQ.0) RETURN                                                 
      IF(ICOND.NE.1)GOTO 13                                             
      N=N+1                                                             
      DO 14 I=1,4                                                       
14    YMAT(N,I)=Y(I)                                                    
      NITE3(N)=0                                                        
13    IF(N.EQ.0)RETURN                                                  
      DO 15 I=1,N                                                       
      YMAT(I,5)=1.D0-YMAT(I,1)-YMAT(I,2)                                
15    YMAT(I,6)=1.D0-YMAT(I,3)-YMAT(I,4)                                
      DO 30 I=1,N                                                       

      xandrea1(1)=ymat(i,k(1))
      xandrea1(2)=ymat(i,k(3))
      xandrea1(3)=ymat(i,k(5))
      xandrea2(1)=ymat(i,k(2))
      xandrea2(2)=ymat(i,k(4))
      xandrea2(3)=ymat(i,k(6))
      call unifac(0,xandrea1,act1,dact,pact)
      call unifac(0,xandrea2,act2,dact,pact)
      acttot(i,1)=dexp(act1(1))
      acttot(i,2)=dexp(act2(1))
      acttot(i,3)=dexp(act1(2))
      acttot(i,4)=dexp(act2(2))
      acttot(i,5)=dexp(act1(3))
      acttot(i,6)=dexp(act2(3))

      WRITE(6,600) I,NITE3(I),(YMAT(I,K(J)),J=1,6)   
	WRITE(3,611) I,NITE3(I),(YMAT(I,K(J)),J=1,6)   
	
      write(8,47) YMAT (i,1),  YMAT (i,2), YMAT (i,5)    !Alfonsina
	 write(8,47) YMAT (i,3),  YMAT (i,4), YMAT (i,6)  !Alfonsina
	 write(8,*) 
	 
   47	FORMAT (8X,F12.6 , 8X,F12.6 , 8X,F12.6) !Alfonsina                  

      if(ig.eq.1) write(6,1600) (acttot(i,j),j=1,6)
   30 continue
 
  600 FORMAT(5X,2I3,2X,3(2P2F9.4,2X))     
  611 FORMAT(5X,2I3,2X,3(2P2F9.4,2X))                                   
      IF(NIC1.GT.0) WRITE(6,900) IC1(1),IC1(NIC1)                       
 1600 format(12x,6(f9.3,1x))
  900 FORMAT(/,' FALSE SOLUTION IN PHASE 1 FOR CALCULATED TIE LINES',I3,' TO',I3,/)                                                       
      IF(NIC2.GT.0) WRITE(6,901) IC2(1),IC2(NIC2)                       
  901 FORMAT(/,' FALSE SOLUTION IN PHASE 2 FOR CALCULATED TIE LINES',I3,' TO',I3,/)                                                       
      IF(IOUT.EQ.6) GOTO 50                                             
      DO 35 I=1,N                                                       
      WRITE(IOUT,600) I,NITE3(I),(YMAT(I,K(J)),J=1,6)                   
    


      if(ig.eq.1) write(iout,1600) (acttot(i,j),j=1,6)
   35 continue

      IF(NIC1.GT.0) WRITE(IOUT,900) IC1(1),IC1(NIC1)                    
      IF(NIC2.GT.0) WRITE(IOUT,901) IC2(1),IC2(NIC2)                    
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SOLVE(Y,DY,NOLD,NEW,NITER,N,NT)                        
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CBISO/NIC1,NIC2,IC1(120),IC2(120)                          
      DIMENSION Y(4),DY(4),AMAT(3,5)                                    
      DIMENSION INO(3)                                                  
      common/nga/nga,mass(12)
      NK=1000  !cambi� de 3 a 100 !Alfonsina                            
      NITER=0                                                           
      NT=1000  !cambi� de 3 a 100 !Alfonsina                            
      DO 1 I=1,3                                                        
    1 IF(Y(I).LT.1.D-14) NT=1000    !Cambie de 10 iteraciones a100 !Alfonsina 0
11    NITER=NITER+1                                                     
      IF(NITER.GT.NT) RETURN                                            
      DO 2 I=1,4                                                        
2     IF(Y(I).LT.0.D0)Y(I)=0.D0                                         
      DO 3 I=1,2                                                        
      Y1(I)=Y(I)                                                        
3     Y2(I)=Y(I+2)                                                      
      Y1(3)=1.D0-Y1(1)-Y1(2)                                            
      Y2(3)=1.D0-Y2(1)-Y2(2)                                            
      IF(Y1(3).LT.0.)Y1(3)=0.                                           
      IF(Y2(3).LT.0.) Y2(3)=0.                                          
      CALL unifac(3,Y1,ACT1,DACT1,PACT)                                 
      CALL unifac(3,Y2,ACT2,DACT2,PACT)                                 
      J=0                                                               
      DO 6 I=1,4                                                        
      IF(I.EQ.NOLD)GOTO 6                                               
      J=J+1                                                             
      INO(J)=I                                                          
6     CONTINUE                                                          
      DO 7 I=1,3                                                        
      DO 7 J=1,2                                                        
      AMAT(I,J)=DACT1(I,J)-DACT1(I,3)                                   
7     AMAT(I,J+2)=DACT2(I,3)-DACT2(I,J)                                 
      DO 8 I=1,3                                                        
      AMAT(I,5)=AMAT(I,NOLD)                                            
      DO 9 J=1,3                                                        
9     AMAT(I,J)=AMAT(I,INO(J))                                          
8     AMAT(I,4)=ACT1(I)-ACT2(I)                                         
      CALL GAUSL(3,5,3,2,AMAT)                                          
      RES=0.D0                                                          
      DO 10 I=1,3                                                       
      Y(INO(I))=Y(INO(I))-AMAT(I,4)                                     
      DY(INO(I))=-AMAT(I,5)                                             
10    RES=RES+AMAT(I,4)**2                                              
      IF(RES.GT.1.D-10)GOTO 11                                          
      IZ=0                                                              
      DO 14 I=1,2                                                       
      IF(Y1(I).LT.1.D-14) IZ=1                                          
   14 IF(Y2(I).LT.1.D-14) IZ=1                                          
      IF(IZ.EQ.1) GOTO 13                                               
      CALL GCON(3,Y1,ACT1,DACT1,ICVEX)                                  
      IF(ICVEX.EQ.1) GOTO 15                                            
      NIC1=NIC1+1                                                       
      IC1(NIC1)=N+1                                                     
   15 CALL GCON(3,Y2,ACT2,DACT2,ICVEX)                                  
      IF(ICVEX.EQ.1) GOTO 13                                            
      NIC2=NIC2+1                                                       
      IC2(NIC2)=N+1                                                     
13    DY(NOLD)=1.D0                                                     
      NEW=NOLD                                                          
      DYMAX=1.D0                                                        
      DO 12 I=1,4                                                       
      IF(DABS(DY(I)).LE.DYMAX)GOTO 12                                   
      NEW=I                                                             
      DYMAX=DABS(DY(I))                                                 
12    CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE GCON(NK,X,ACT,DACT,ICVEX)                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION X(3),DG(2),DDG(2,2),ACT(3),DACT(10,10)                  
      ICVEX=1                                                           
      DO 1 I=1,NK                                                       
    1 IF(ACT(I).LT.1.D-10) ACT(I)=1.D-10                                
      DO 5 I=1,NK                                                       
      DO 5 J=1,NK                                                       
    5 DACT(I,J)=DACT(I,J)/ACT(I)                                        
      IF(NK.EQ.3) GOTO 9                                                
      DDG(2,2)=DACT(2,2)-DACT(1,2)-DACT(2,1)+DACT(1,1)                  
      GOTO 30                                                           
9     DO 20 I=2,NK                                                      
      II=I-1                                                            
      DO 20 J=2,NK                                                      
      JJ=J-1                                                            
   20 DDG(II,JJ)=DACT(I,J)-DACT(1,J)-DACT(I,1)+DACT(1,1)                
      IF(X(1).LE.1.D-12.OR.X(2).LE.1.D-12) GOTO 30                      
      DET=DDG(1,1)*DDG(2,2)-DDG(2,1)*DDG(2,1)                           
      IF(DET.LE.0.D0.OR.DDG(1,1).LE.0.D0.OR.DDG(2,2).LE.0.D0) ICVEX=-1  
      GOTO 100                                                          
   30 CONTINUE                                                          
      IF(DDG(2,2).LE.0.D0) ICVEX=-1                                     
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      FUNCTION RUND(Y,DY,S,IRUND)                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      X=Y+S*DY+1.D-8*DY**2                                              
      IX=IRUND*X                                                        
      Z=DFLOAT(IX)/IRUND-Y                                              
      RUND=DABS(Z)                                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE TERM(Y,DMAT,ICOND,NEW)                                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION Y(4),DMAT(4,4),A(4)                                     
      IF(Y(1).LT.1.D-14.OR.Y(3).LT.1.D-14) GOTO 1                       
      IF(Y(1)+Y(2).GT.1.D0.OR.Y(3)+Y(4).GT.1.D0) GOTO 2                 
      IF(Y(1)+Y(2)-.01D0.LT.Y(3)+Y(4).AND.Y(1)-.01D0.LT.Y(3))GOTO 3     
      RETURN                                                            
1     ICOND=2                                                           
      DS=DMAT(1,1)/(DMAT(1,1)-Y(1))                                     
      DO 5 I=1,4                                                        
5     Y(I)=DMAT(I,1)+DS*(Y(I)-DMAT(I,1))                                
      Y(1)=0.D0                                                         
      NEW=1                                                             
      RETURN                                                            
2     ICOND=-2                                                          
      RETURN                                                            
3     ICOND=1                                                           
      ND=2+NEW                                                          
      IF(ND.GT.4)ND=ND-4                                                
      DO 6 I=1,4                                                        
6     A(I)=DMAT(NEW,I)-DMAT(ND,I)                                       
      DS=0.D0                                                           
      NITER=0                                                           
7     NITER=NITER+1                                                     
      IF(NITER.LE.10)GOTO 8                                             
      ICOND=-1                                                          
      RETURN                                                            
8     F=((A(4)*DS+A(3))*DS+A(2))*DS+A(1)                                
      DF=(3.D0*A(4)*DS+2.D0*A(3))*DS+A(2)                               
      DF=-F/DF                                                          
      DS=DS+DF                                                          
      IF(DABS(DF).GT.1.D-6)GOTO 7                                       
      DO 9 I=1,4                                                        
      Y(I)=DMAT(I,4)                                                    
      DO 9 J=1,3                                                        
9     Y(I)=Y(I)*DS+DMAT(I,4-J)                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE SOLBIN                                                 
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CACT/X1(10),X2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CY/Y13,Y21,STEP                                            
      COMMON/COUT/IOUT                                                  
      DIMENSION DMAT(2,3)                                               
      common/nga/nga,mass(12)

      NITER=0                                                           
      X1(1)=1.D0-Y13/100.D0                                             
      X2(1)=Y21/100.D0                                                  
   10 NITER=NITER+1                                                     
      IF(NITER.GT.10) GOTO 50                                           
      IF(X1(1).LT.0.D0) X1(1)=0.D0                                      
      IF(X2(1).LT.0.D0) X2(1)=0.D0                                      
      X1(2)=1.D0-X1(1)                                                  
      X2(2)=1.D0-X2(1)                                                  
      CALL unifac(3,X1,ACT1,DACT1,PACT)                                 
      CALL unifac(3,X2,ACT2,DACT2,PACT)                                 
      DO 20 I=1,2                                                       
      DMAT(I,1)=DACT1(I,1)-DACT1(I,2)                                   
      DMAT(I,2)=DACT2(I,2)-DACT2(I,1)                                   
   20 DMAT(I,3)=ACT1(I)-ACT2(I)                                         
      CALL GAUSL(2,3,2,1,DMAT)                                          
      RES=DMAT(1,3)**2+DMAT(2,3)**2                                     
      X1(1)=X1(1)-DMAT(1,3)                                             
      X2(1)=X2(1)-DMAT(2,3)                                             
      IF(RES.GT.1.D-20) GOTO 10                                         
   50 CONTINUE                                                          
      WRITE(6,603)                                                      
      IF(IOUT.NE.6) WRITE(IOUT,603)                                     
      IF(IOUT.NE.6) WRITE(IOUT,604) X1(1),X2(1),X1(2),X2(2)             
      WRITE(6,604) X1(1),X2(1),X1(2),X2(2)                              
  603 FORMAT(///,5X,'** BINARY SOLUBILITIES IN MOLE FRACTIONS **',//,11X,'COMPONENT 1',15X,'COMPONENT 2',/)                               
  604 FORMAT(2(2X,2P2D12.2)//)                                          
      CALL GCON(2,X1,ACT1,DACT1,ICVEX)                                  
      IF(IOUT.NE.6.AND.ICVEX.EQ.-1) WRITE(IOUT,601)                     
      IF(ICVEX.EQ.-1) WRITE(6,601)                                      
      CALL GCON(2,X2,ACT2,DACT2,ICVEX)                                  
      IF(IOUT.NE.6.AND.ICVEX.EQ.-1) WRITE(IOUT,602)                     
      IF(ICVEX.EQ.-1) WRITE(6,602)                                      
  601 FORMAT(' FALSE SOLUTION IN PHASE 1')                              
  602 FORMAT(' FALSE SOLUTION IN PHASE 2')                              
      RETURN                                                            
      END                                                               
      SUBROUTINE FUNC(N,M,NDIF,X,SSQ)                                   
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
      COMMON/CMARQ/GRAD(2),XJTJ(2,2)                                    
      COMMON/CUFAC/NK,NG,P(10,10),T                                     
      COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
      DIMENSION X(2),F(2)                                               
      JD=0                                                              
      IF(NDIF.EQ.1) JD=4                                                
      P(1,2)=X(1)*300.D0                                                
      P(2,1)=X(2)*300.D0                                                
      CALL PARAM2                                                       
      SSQ=0.                                                            
      IF(NDIF.EQ.0) GOTO 11                                             
      DO 10 I=1,2                                                       
      GRAD(I)=0.                                                        
      DO 10 J=1,2                                                       
   10 XJTJ(I,J)=0.                                                      
   11 CONTINUE                                                          
      DO 21 L=1,5                                                       
      Y1(1)=XC(L)                                                       
      Y1(2)=1.D0-XC(L)                                                  
      CALL unifac(JD,Y1,ACT1,DACT1,PACT)                                
      DO 17 I=1,2                                                       
      F(I)=ACT1(I)-GE(L,I)                                              
      GC(L,I)=ACT1(I)                                                   
   17 SSQ=SSQ+F(I)*F(I)                                                 
      IF(JD.EQ.0) GOTO 21                                               
      DO 19 I=1,2                                                       
      DO 20 J=1,2                                                       
      GRAD(J)=GRAD(J)+F(I)*PACT(I,J)                                    
      DO 20 K=1,2                                                       
   20 XJTJ(J,K)=XJTJ(J,K)+PACT(I,J)*PACT(I,K)                           
   19 CONTINUE                                                          
   21 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE MARQ(FUNC,N,M,X,XLAMB,FAC,EPSG,MAXF)                   
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/COUT/IOUT                                                  
      COMMON/CMARQ/GRAD(2),XJTJ(2,2)                                    
      COMMON/CIPR/IPR                                                   
      DIMENSION X(2),Y(2),XNY(2),A(2,2),DX(2)                           
      IEVAL=0                                                           
      ISTOP=0                                                           
      IER=0                                                             
      XMAXL=XLAMB*1.D+4                                                 
      ITER=0                                                            
      IF(IOUT.NE.6.AND.IPR.EQ.1) WRITE(IOUT,603)                        
      IF(IPR.EQ.1) WRITE(6,603)                                         
      CALL FUNC(N,M,1,X,SRES)                                           
      IEVAL=IEVAL+1                                                     
      SSQ=SRES                                                          
      IF(IPR.EQ.1.AND.IOUT.NE.6) WRITE(IOUT,601) ITER,SSQ               
      IF(IPR.EQ.1) WRITE(6,601) ITER,SSQ                                
   10 CONTINUE                                                          
      IF(IEVAL.NE.1) CALL FUNC(N,M,1,X,SRES)                            
      GNORM=0.D0                                                        
      DO 140 I=1,N                                                      
  140 GNORM=GNORM+GRAD(I)**2                                            
      GNORM=DSQRT(GNORM)                                                
      IF(GNORM.LT.EPSG) ISTOP=ISTOP+1                                   
      IF(ISTOP.GT.0) GOTO 1000                                          
      ITER=ITER+1                                                       
   49 CONTINUE                                                          
      IF(IEVAL.GT.MAXF) GOTO 998                                        
      DO 41 I=1,N                                                       
      DO 40 J=1,N                                                       
   40 A(I,J)=XJTJ(I,J)                                                  
   41 A(I,I)=A(I,I)+XLAMB                                               
      CALL CHOL(N,A)                                                    
      Y(1)=-GRAD(1)/A(1,1)                                              
      DO 81 I=2,N                                                       
      SUM=0.D0                                                          
      II=I-1                                                            
      DO 80 J=1,II                                                      
   80 SUM=SUM+A(I,J)*Y(J)                                               
   81 Y(I)=(-GRAD(I)-SUM)/A(I,I)                                        
      DX(N)=Y(N)/A(N,N)                                                 
      DO 85 I=2,N                                                       
      II=N-I+1                                                          
      SUM=0.D0                                                          
      III=II+1                                                          
      DO 84 J=III,N                                                     
   84 SUM=SUM+A(J,II)*DX(J)                                             
   85 DX(II)=(Y(II)-SUM)/A(II,II)                                       
      DO 90 I=1,N                                                       
   90 XNY(I)=X(I)+DX(I)                                                 
      CALL FUNC(N,M,0,XNY,SRES)                                         
      IEVAL=IEVAL+1                                                     
      SSQNY=SRES                                                        
      SQ1=0.D0                                                          
      SQ2=0.D0                                                          
      DO 110 I=1,N                                                      
      Y(I)=XNY(I)*300.D0                                                
      SQ1=SQ1+DX(I)**2                                                  
      SQ2=SQ2-DX(I)*GRAD(I)                                             
  110 CONTINUE                                                          
      CCAL=SSQ-SSQNY                                                    
      CPRE=SQ2+XLAMB*SQ1                                                
      CALPRE=CCAL/(CPRE+1.D-14)                                         
      IF(IPR.EQ.1) WRITE(6,601) ITER,SSQNY,Y(1),Y(2),GNORM,XLAMB,CALPRE 
      IF(IOUT.NE.6.AND.IPR.EQ.1) WRITE(IOUT,601) ITER,SSQNY,Y(1),Y(2),  GNORM,XLAMB,CALPRE                                                
      IF(SSQ-SSQNY.GT..75*CPRE) XLAMB=XLAMB/FAC                         
      IF(SSQ-SSQNY.LT..25*CPRE) XLAMB=XLAMB*FAC                         
      IF(XLAMB.GT.XMAXL) GOTO 999                                       
      IF(SSQNY-SSQ) 120,120,119                                         
  119 CONTINUE                                                          
      IF(SSQNY-SSQ.GT.DABS(SSQ)*.5D0) XLAMB=XLAMB*FAC                   
      GOTO 49                                                           
  120 CONTINUE                                                          
      IF(SSQ-SSQNY.GT.DABS(SSQ)*.8D0) XLAMB=XLAMB/FAC                   
      DO 130 I=1,N                                                      
  130 X(I)=XNY(I)                                                       
      SSQ=SSQNY                                                         
      GOTO 10                                                           
  998 IER=1                                                             
      GOTO 1000                                                         
  999 IER=2                                                             
      GOTO 1000                                                         
  601 FORMAT(1X,I3,2X,D12.4,2F10.3,5X,2D12.4,3X,F10.2)                  
  602 FORMAT(//,' ISTOP= ',I2,5X,'IER = ',I2,5X,'IEVAL = ',I5,//)       
  603 FORMAT(///,'  ** ITERATIONCOURSE, UNIQUAC-PARAMETERS FROM UNIFAC **'/)                                                              
 1000 CONTINUE                                                          
      IF(IOUT.NE.6.AND.IPR.EQ.1) WRITE(IOUT,602) ISTOP,IER,IEVAL        
      IF(IPR.EQ.1) WRITE(6,602)ISTOP, IER, IEVAL                        
      RETURN                                                            
      END                                                               
      SUBROUTINE CHOL(N,A)                                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION A(2,2)                                                  
      DO 50 I=1,N                                                       
      I1=I-1                                                            
      IF(I1.EQ.0) GOTO 30                                               
      DO 20 J=I,N                                                       
      DO 20 K=1,I1                                                      
   20 A(I,J)=A(I,J)-A(I,K)*A(J,K)                                       
   30 IF(A(I,I).LT.1.D-14) A(I,I)=1.D-14                                
      A(I,I)=DSQRT(A(I,I))                                              
      IF(I.EQ.N) GOTO 100                                               
      J1=I+1                                                            
      DO 50 J=J1,N                                                      
   50 A(J,I)=A(I,J)/A(I,I)                                              
  100 RETURN                                                            
      END                                                               
!C  SUBROUTINE GAUSL SOLVES N LINEAR ALGEBRAIC EQUATIONS BY GAUSS        
!C  ELIMINATION WITH ROW PIVOTING                                        
!C  TO SOLVE THE PROBLEM QX=U, WHERE Q IS A NXN MATRIX AND U IS NXNS,    
!C  ONE PLACES Q IN THE FIRST N COLUMNS OF A AND U IS PLACED IN THE      
!C  FOLLOWING NS COLUMNS.                                                
!C  THE PROGRAM RETURNS X=Q**(-1)*U AT THE PREVIOUS POSITION OF U.       
!C  *                                                                    
!C  ND IS THE ROW DIMENSION AND NCOL IS THE COLUMN DIMENSION OF A.       
!C  BOTH MUST BE TRANSFERRED TO THE SUBROUTINE.                          
!C  *****************                                                    
!C                                                                       
      SUBROUTINE GAUSL(ND,NCOL,N,NS,A)                                  
                                                                     
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(ND,NCOL)                                              
      N1=N+1                                                            
      NT=N+NS                                                           
      IF (N .EQ. 1) GO TO 50                                            
!C      START ELIMINATION                                                
!C                                                                       
!C                                                                       
      DO 10 I=2,N                                                       
      IP=I-1                                                            
      I1=IP                                                             
      X=DABS(A(I1,I1))                                                  
      DO 11 J=I,N                                                       
      IF (DABS(A(J,I1)) .LT. X) GO TO 11                                
      X=DABS(A(J,I1))                                                   
      IP=J                                                              
   11 CONTINUE                                                          
      IF (IP .EQ. I1) GO TO 13                                          
!C                                                                       
!C     ROW INTERCHANGE                                                   
!C                                                                       
      DO 12 J=I1,NT                                                     
      X=A(I1,J)                                                         
      A(I1,J)=A(IP,J)                                                   
   12 A(IP,J)=X                                                         
   13 DO 10 J=I,N                                                       
      X=A(J,I1)/A(I1,I1)                                                
      DO 10 K=I,NT                                                      
   10 A(J,K)=A(J,K) - X*A(I1,K)                                         
!C                                                                       
!C      ELIMINATION FINISHED, NOW BACKSUBSTITUTION                       
!C                                                                       
   50 DO 20 IP=1,N                                                      
      I=N1-IP                                                           
      DO 20 K=N1,NT                                                     
      A(I,K) = A(I,K)/A(I,I)                                            
      IF (I .EQ. 1) GO TO 20                                            
      I1=I-1                                                            
      DO 25 J=1,I1                                                      
   25 A(J,K) = A(J,K) - A(I,K)*A(J,I)                                   
   20 CONTINUE                                                          
      RETURN                                                            
      END                                            
      

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      double precision a(np,np),b(n)
      INTEGER i,ii,j,ll
      double precision sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n

          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d1=-d1
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)

          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END      
                         
!C******************************* F I N ***************************************
!C
