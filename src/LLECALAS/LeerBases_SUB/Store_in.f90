SUBROUTINE Store_In (model,nintt)
    !C-----------------------------------------------------------------------
    !C      ESTA SUBRUTINA CARGA TODOS LOS PARAMETROS DEL BANCO INTRCN QUE SE
    !C      USAN EN MOLDES. LAS PROPIEDADES SE ALMACENAN EN COMMONS.
    !C-----------------------------------------------------------------------
    !C
       use InputData
       !PARAMETER (NG=70,NMG=150)
        implicit none
    
    !INTEGER
        integer::i,j,id_g1,id_g2,model
    
    !COMMONS    
    
        integer::NPINT(NINT),NINTT(NINT),NUMINT
    
       ! real*8,dimension(NMG,NMG)::A !,kstr(NMG,NMG),kdot(NMG,NMG),alpha(NMG,NMG)
    
    !FUNCIONES
        real*8,external::Leer_In !Leer_Alpha, Leer_Kapa,
    
       ! if(InputProblem%Model /= 3)then
        I=1
        do while (NINTT(I).ne.0)
           J=1
           do while (NINTT(J).ne.0)
    
              Aint1(I,J) = Leer_In(NINTT(I),NINTT(J),IPAREQ)
    
              J=J+1
           enddo
           I=I+1
        enddo
    
    !    else !para la GC
    !
    !    I=1
    !    do while (NINTT(I).ne.0)
    !       J=I
    !       do while (NINTT(J).ne.0)
    !
    !             if(Nintt(i) > Nintt(j))then
    !                 id_g1 = Nintt(j)
    !                 id_g2 = Nintt(i)
    !             else
    !                 id_g1 = Nintt(i)
    !                 id_g2 = Nintt(j)
    !             endif
    !             
    !             kstr(id_g1,id_g2) = Leer_Kapa(id_g1,id_g2) !; kstr(id_g2,id_g1) = param(1)
    !             if(isnan(kstr(id_g1,id_g2))) kstr(id_g1,id_g2) = 1
    !             kstr(id_g2,id_g1) = kstr(id_g1,id_g2)
    !             
    !             kdot(id_g1,id_g2) = Leer_Kapa(id_g2,id_g1) !; kdot(id_g2,id_g1) = param(2)
    !             if(isnan(kdot(id_g1,id_g2))) kdot(id_g1,id_g2) = 0
    !             kdot(id_g2,id_g1) = kdot(id_g1,id_g2) 
    !             
    !             alpha(id_g1,id_g2) = Leer_Alpha(id_g1,id_g2)
    !             if(isnan(alpha(id_g1,id_g2))) alpha(id_g1,id_g2) = 0
    !             alpha(id_g2,id_g1) = Leer_Alpha(id_g2,id_g1)
    !             if(isnan(alpha(id_g2,id_g1))) alpha(id_g2,id_g1) = 0
    !
    !          J=J+1
    !       enddo
    !       I=I+1
    !    enddo
    !    
    !    endif    
        
        RETURN
        ENDSUBROUTINE Store_In