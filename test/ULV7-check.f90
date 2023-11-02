program check
    use do_tests
    use ULV7_values
    use flash 
    use fobjtype
    use CUFAC
    
    implicit none
    
    integer::i
    integer::variables
    real*8::fmin
    real*8,allocatable,dimension(:)::x
    real*4,dimension(10)::zf_check
    real*4,dimension(3)::ksalida_check
    character(len=*),parameter :: name_filename =  "test/ULV7.dat"
    double precision,external::conversion
    double precision,external::praxis_n,f_n ,newton 

    print *,""
    print *, test_run//"ULV7-Test"
    if (ULV7_check) then
        continue
    else 
        print *, test_disabled
        goto 999
    endif
    if (pause_test) pause
    OPEN (unit=333,file='salida.OUT')
    call leer_input_flash(name_filename)
    !call leer_input_flash_ULV7()
    variables = 1
    allocate(x(variables))
    x(1) = 0.0000001
    if (.false.) then
        agextr = .true.
        aglim = 1.D-4
    else
        agextr = .false.
    endif
    call genDatExtr_ULV7(x,variables,2,10,.false.,zf_check,ksalida_check)
    close (unit=333)
    
    !Check if everything went OK
    if (abs(zf_check(1) - ULV7_data(1))>1E-8 .or. isnan(zf_check(1)))&
        &ERROR STOP "Zf of component 1 has changed"
    if (abs(zf_check(2) - ULV7_data(2))>1E-8 .or. isnan(zf_check(2)))&
        &ERROR STOP "Zf of component 2 has changed"
    if (abs(zf_check(3) - ULV7_data(3))>1E-8 .or. isnan(zf_check(3)))&
        &ERROR STOP "Zf of component 3 has changed"
    if (abs(zf_check(4) - ULV7_data(4))>1E-8 .or. isnan(zf_check(4)))& 
        &ERROR STOP "Zf of component 4 has changed"
    if (abs(zf_check(5) - ULV7_data(5))>1E-8 .or. isnan(zf_check(5)))& 
        &ERROR STOP "Zf of component 5 has changed"
    if (abs(zf_check(6) - ULV7_data(6))>1E-8 .or. isnan(zf_check(6)))& 
        &ERROR STOP "Zf of component 6 has changed"
    if (abs(zf_check(7) - ULV7_data(7))>1E-8 .or. isnan(zf_check(7)))& 
        &ERROR STOP "Zf of component 7 has changed"
    if (abs(ksalida_check(1) - ULV7_data(8))>1E-8 .or. isnan(ksalida_check(1)))& 
        &ERROR STOP "K1 has changed"
    if (abs(ksalida_check(2) - ULV7_data(9))>1E-8 .or. isnan(ksalida_check(2)))& 
        &ERROR STOP "K2 has changed"
    if (abs(ksalida_check(3) - ULV7_data(10))>1E-8.or. isnan(ksalida_check(3)))& 
        &ERROR STOP "K3 has changed"

    print *, test_ok
    999 continue
end program check

subroutine leer_input_flash_ULV7()
    use InputData
    use flash
    use CUFAC 
    implicit none    
    integer::i,j,k
    
    
    call open_file_name_ULV7()
    
    OPEN (UNIT=2,FILE=flashInput_name,status='OLD',FORM='FORMATTED')
    READ(2,501) NTEXT   
    501 FORMAT(36A2)  
    
    READ(2,*) ICALC,modelo,IPRm,IOUTm,NOVAPm,igm, ipareq     
    
    call open_database(modelo)
    CALL PARIN2(NKK,NGG,Pxx,Txx)
    
    IF(NOVAPm/=0) then                                             
        DO 6 J=1,NKK                                                        
                                      
    6       READ(2,*) (ANT(K,J),K=1,3)                                        
        DO 7 J=1,NKK                                                        
            ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
    7       ANT(2,J)=2.302585*ANT(2,J)                                        
    endif   
    READ(2,502) T,P
    502 FORMAT(2F10.2)  
    READ(2,*) (Z(I),I=1,NKK)  
    close(unit=2)
endsubroutine leer_input_flash_ULV7

subroutine open_file_name_ULV7()
    use InputData
    implicit none
    
    integer::parameters
    
    !apertura de bases de datos generadas desde excel
    OPEN (UNIT=1,FILE='test/ULV7-name.dat',status='OLD',FORM='FORMATTED')
    read(1,*)parameters 
    read(1,"(A36)") flashInput_name
    flashInput_name = flashInput_name(2:len_trim(flashInput_name)-1)
    flashInput_name = "test/ULV7-"//flashInput_name

    if (parameters==1)then
        call leerBases()
        stop
    endif    
    CLOSE (UNIT=1)  

endsubroutine open_file_name_ULV7

subroutine genDatExtr_ULV7(x,n,extI,extS,temp,zf_check,ksalida_check)
    use Flash
    use flashout

    implicit none
    logical::temp,proc
    integer::extI,extS,camb
    integer,intent(in)::n
    real*8::x(n)
    real*8,dimension(size(z))::zf
    doubleprecision::obj
    real*4,intent(out),dimension(3)::ksalida_check
    real*4,intent(out),dimension(10)::zf_check
    
    integer::i
    
    i = extI
    camb = 1
    proc = .true.
    x(1) = i
    call conversion(zf)
    !Escribe en el documento de salida las condiciones
    !write(*,'(A,F10.8)') "Acido: ",zf(1)
    !write(*,'(A,F10.8)') "Gly: ", zf(2)
    !write(*,'(A,F10.8)') "Mono: ", zf(3)
    !write(*,'(A,F10.8)') "Di: ", zf(4)
    !write(*,'(A,F10.8)') "Tri: ", zf(5)
    !write(*,'(A,F10.8)') "Agua: ", zf(6)
    !write(*,'(A,F10.8)') "Tol: ", zf(7)
    i = i + camb
    !write(*,'(A,F10.8)') "K1: ", ksalida(1)
    !write(*,'(A,F10.8)') "K2: ", ksalida(2)
    !write(*,'(A,F10.8)') "K3: ", ksalida(3)
    ksalida_check(1)=ksalida(1)
    ksalida_check(2)=ksalida(2)
    ksalida_check(3)=ksalida(3)
    zf_check = zf

endsubroutine genDatExtr_ULV7