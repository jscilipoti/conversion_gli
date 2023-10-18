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