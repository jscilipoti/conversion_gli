integer function mainsgfunc (i1,ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el numero de grupo principal correspon-
!c      diente al subgrupo i1 de la tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
    integer:: il,ipareq
	irec1=i1+(ipareq-1)*150
      read (14,10,rec=irec1) main
      mainsgfunc = main
  10  format (i4)
      return
endfunction mainsgfunc


SUBROUTINE LEEPAR (J,IREC1,IPAREQ,NGRUPA,ENASST,RKASST) 
!c-----------------------------------------------------------------
!c     Lee los par�metros de las bases de datos PARVOLAS.MDS 
!c     (volumen de asociaci�n, UNIT=15) y PARENEAS.MDS (energ�a de 
!c     asociaci�n, UNIT=16)  
!c-----------------------------------------------------------------
IMPLICIT real*8 (A-H,O-Z)
PARAMETER (NMG=150)
INTEGER J,IREC1,IPAREQ
DIMENSION NGRUPA(NMG)
      READ(16,502,REC=IREC1)FS,(ENASST,I=1,MAINSGfunc(NGRUPA(J),IPAREQ))
      READ(15,502,REC=IREC1)FS,(RKASST,I=1,MAINSGfunc(NGRUPA(J),IPAREQ))     
 502  FORMAT(a8,70d12.5)      
END

subroutine ab_ban1(model)
!-------------------------------------------------------------------
!      Esta subrunina ABRE los bancos de datos siguientes:
!     - INTRCN.MDS      (UNIT=13)
!     - INTRCNAS.MDS    (UNIT=13)
!     - GRUPOSRAM.MDS   (UNIT=14)
!     - PARVOLAS.MDS    (UNIT=15)  
!     - PARENEAS.MDS    (UNIT=16)     
!-------------------------------------------------------------------
!
    !use input_data, only:model
      COMMON/AS/ASOC
      LOGICAL ASOC
      integer::model, mod
   
      
    !CALL MODEL(mod)
    if(model /= 3)then
        if(model==1)then
            open (unit=13,file='src/database/intrcn.mds',status='old',&
                  access='direct',form='formatted',recl=850)   
        else     
            open (unit=13,file='src/database/intrcnas.mds',status='old',&
                  access='direct',form='formatted',recl=850)        
        endif
        open (unit=14,file='src/database/gruposram.mds',status='old',&
     	      access='direct',form='formatted',recl=300)
        open (unit=15,file='src/database/parvolas.mds',status='old',&
              access='direct',form='formatted',recl=850)
        open (unit=16,file='src/database/pareneas.mds',status='old',&
              access='direct',form='formatted',recl=850)    
    else

        open (unit=14,file='src/database/gruposramgc.mds',status='old',&
     	      access='direct',form='formatted',recl=263)
        open (unit=13,file='src/database/intrcngcalpha.mds',status='old',&
              access='direct',form='formatted',recl=730)
        open (unit=16,file='src/database/intrcngckapa.mds',status='old',&
              access='direct',form='formatted',recl=730)    
        
    endif
      

      return
endsubroutine ab_ban1



SUBROUTINE Model(m)
  !  use Input
      COMMON/AS/ASOC
      LOGICAL ASOC
      INTEGER M
      
      IDEV=6
      ASOC=.FALSE.
      write (idev,21)  ! Choose the MODEL
61    write (idev,30) !"> "
      READ (5,981,ERR=61) M
      iF (M<1.OR.M>3) GOTO 61	
      IF(M.EQ.2)ASOC=.TRUE.
  !  InputProblem%model = M
 21   format (//' Choose the model:',//,&
             15x,'UNIFAC  : 1',/,&
             15x,'A-UNIFAC: 2',/,&
             15x,'GC-EOS  : 3')
 30   format (1x,/,60x,'> ',$)
 981  FORMAT (I1)
endsubroutine


SUBROUTINE BUSCARAS (N,VECTOR,tam,VL)
!c---------------------------------------------------------------------
!c     Comprueba la existencia del grupo N en VECTOR
!c--------------------------------------------------------------------
    implicit none
    integer::tam,n,i
    INTEGER VECTOR(tam)
    LOGICAL VL

    VL=.FALSE.
    I=0
    DO WHILE (.NOT.VL)
        I=I+1
	  IF (N.EQ.VECTOR(I)) VL=.TRUE.
	  IF (I.EQ.tam) EXIT
	ENDDO
endsubroutine

SUBROUTINE BUSCAR (I,IPAREQL,FF,NUM)
      CHARACTER*8 FS,FF,MG
      LOGICAL IPAREQL
!C
      IPAREQL=.FALSE.
      num=(I-1)*150
      fs= ''
      MG= ''
      DO WHILE (MG.NE.'fin')
	  num=num+1
	  read (14,10,rec=num)MG,fs
	  IF (FS.EQ.FF) THEN
	      IPAREQL=.TRUE.
	      EXIT
	  ENDIF
	enddo
!c
!c
  10  format (4x,2a8) 
      return
endsubroutine buscar