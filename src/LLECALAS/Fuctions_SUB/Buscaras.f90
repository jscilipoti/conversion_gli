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