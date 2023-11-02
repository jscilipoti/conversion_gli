program check
    !ANT, NGG, Pxx, Txx: No fueron checkeados
    use do_tests
    use InputData, only:&
    &name, name_maxlen, ICALC, modelo, IPRm, IOUTm, NOVAPm, igm, ipareq,&
    &ANT,NTEXT
    use flash, only: P, T, Z
    use CUFAC, only: NKK, NGG, Pxx, Txx
    
    implicit none
    real*8,dimension(10) :: z_check
    integer :: i 

    print *,""
    print *, test_run//"leer_input_flash-Test"
    if (leer_input_flash_check) then
        if (pause_test) pause
        call leer_input_flash("test/name.dat")
        !Check if everything went OK
        if (name /= "test/llecalas2.dat")& 
            &ERROR STOP "'name' was not read correctly."
        if (abs(ICALC - 0) > 1E-8)&
            &ERROR STOP "'ICALC' was not read correctly."
        if (abs(modelo - 0) > 1E-8)&
            &ERROR STOP "'modelo' was not read correctly."
        if (abs(IPRm - 0) > 1E-8)&
            &ERROR STOP "'IPRm' was not read correctly."
        if (abs(IOUTm - 1) > 1E-8)&
            &ERROR STOP "'IOUTm' was not read correctly."
        if (abs(NOVAPm - 0) > 1E-8)&
            &ERROR STOP "'NOVAPm' was not read correctly."
        if (abs(igm - 1) > 1E-8)&
            &ERROR STOP "'igm' was not read correctly."
        if (abs(ipareq - 2) > 1E-8)&
            &ERROR STOP "'ipareq' was not read correctly."
        if (abs(T - 358.15) > 1E-5)&
            &ERROR STOP "'T' was not read correctly."
        if (abs(P - 1.) > 1E-8)&
            &ERROR STOP "'P' was not read correctly."
        z_check=Z-(/0.65,0.13,1E-8,1E-8,1E-8,1E-8,1E-8,0.,0.,0./)
        do i = 1 , size(z_check)
            if (abs(z_check(i)) > 1E-5) &
            &ERROR STOP "'Z' was not read correctly."
        enddo
        if (abs(NKK - 7) > 1E-8) ERROR STOP "'NKK' was not read correctly."
        if (NTEXT /= "llecalas2") ERROR STOP "'NTEXT' was not read correctly."


        print *, test_ok
    else 
        print *, test_disabled
    endif
end program check
    