SUBROUTINE CR_PUNTF (IGRUP,NPUNT,NGRUP,NPINT,NINTT)
    !--------------------------------------------------------------------------
    !   ESTA SUBRUTINA CREA EL VECTOR NGRUP, QUE CONTIENE TODOS LOS GRUPOS 
    !   PRESENTES EN EL CASO QUE SE ESTA CORRIENDO.
    !	NPUNT es el vector que contiene las posiciones que los distintos
    !	grupos ocupan en el vector NGRUP. ej. NPUNT(9)=1 ; NPUNT(12)=2 ...
    !--------------------------------------------------------------------------
    !
        use InputData
        implicit none
    
    !Variables de ENTRADA
        integer,intent(in)::igrup
    !Variables de ENTRADA/SALIDA      
        integer,dimension(NMG)::NPUNT,NGRUP
        integer,dimension(NINT):: NPINT,NINTT
    !Variables INTERNAS      
        integer::j,k,maingrup
    !COMMONS      
        integer,dimension(NMG)::MAIN
        real*8,dimension(NMG)::R,Q
    !FUNCTIONS    
        integer,external::mainsgfunc
    
    !SENTENCIAS
        if(npunt(igrup)==0)then             !PREGUNTA SI EST� CARGADO EL SUBGRUPO
          do j=1,NMG                        !BUSCA PRIMER POSICI�N VAC�A
              if(ngrup(J)==0)exit
          enddo
          npunt(igrup) = J
          ngrup(J) = igrup
          maingrup = mainsgfunc(ngrup(J),ipareq)!BUSCA A QU� GRUPO PPAL PERTENECE
          if(npint(maingrup)==0)then        !PREGUNTA SI EST� CARGADO EL GRUPO PPAL
              do K=1,NINT                   !BUSCA PRIMER POSICI�N VAC�A
                  if(nintt(K)==0)exit
              enddo
              npint(maingrup) = k
              nintt(K) = maingrup
              main(k) = maingrup
          endif      
        endif
        
    endsubroutine CR_PUNTF