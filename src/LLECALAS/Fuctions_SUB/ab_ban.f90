subroutine ab_ban1(model)
    !-------------------------------------------------------------------
    !   This subroutine opens each database:
    !   - INTRCN.MDS      (UNIT=13)
    !   - INTRCNAS.MDS    (UNIT=13)
    !   - GRUPOSRAM.MDS   (UNIT=14)
    !   - PARVOLAS.MDS    (UNIT=15)  
    !   - PARENEAS.MDS    (UNIT=16)     
    !-------------------------------------------------------------------
    !
    use iso_fortran_env, only: int8

    implicit none
    
    integer(kind=int8), intent(in) :: model

    if(model /= 3) then
        if(model == 1) then !if model = 1(UNIQUAC)
            open (unit=13, file='src/database/intrcn.mds', status='old',&
            access='direct', form='formatted', recl=850)   
        else !if model = 0 (UNIFAC) or 2(A-UNIFAC)    
            open (unit=13, file='src/database/intrcnas.mds', status='old',&
            access='direct',form='formatted',recl=850)        
        endif

        !if model = 0(A-UNIFAC) or 1(UNIQUAC) or 2(A-UNIFAC)
        open (unit=14, file='src/database/gruposram.mds', status='old', &
            access='direct',form='formatted',recl=300)
        open (unit=15, file='src/database/parvolas.mds', status='old', &
            access='direct', form='formatted', recl=850)
        open (unit=16, file='src/database/pareneas.mds', status='old', &
            access='direct', form='formatted', recl=850)

    else !if model = 3(GC)
        open (unit=14, file='src/database/gruposramgc.mds', status='old',&
            access='direct', form='formatted', recl=263)
        open (unit=13, file='src/database/intrcngcalpha.mds', status='old',&
            access='direct', form='formatted', recl=730)
        open (unit=16, file='src/database/intrcngckapa.mds', status='old',&
            access='direct', form='formatted', recl=730)    
    endif
    return
    
    end subroutine ab_ban1