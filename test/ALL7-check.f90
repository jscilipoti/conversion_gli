program check
    use do_tests
    use flash 
    use fobjtype
    use CUFAC
    
    implicit none
    
    integer::i
    integer::variables
    real*8::fmin
    real*8,allocatable,dimension(:)::x
    
    double precision,external::conversion
    double precision,external::praxis_n,f_n ,newton 
    
    print *,""
    print *, test_run//"ALL7-Test"
    
    if (ALL7_check) then
        continue
    else 
        print *, test_disabled
        goto 999
    endif
    pause

    OPEN (unit=333,file='salida.OUT')
    call leer_input_flash_ALL7()

    
    variables = 1
    allocate(x(variables))
    
    x(1) = 0.0000001
    
    !Ve en el caso de la extraccion del agua    
    if (.false.) then
        agextr = .true.
        aglim = 1.D-4
    else
        agextr = .false.
    endif
    
    !Optimiza la funcion f_n en el documento Fobj
    !FMIN = praxis_n ( 3.D-4, 1.D1, variables, 3, x, f_n )
    
    !Genera una lista de concentraciones segun la funcion genDat que esta en el documento Conversion_f     
    !call genDat(x,variables)
    
    call genDatExtr(x,variables,2,10,.false.)
    
    !Genera una espacion de datos de la funcion objetivo en todos los valores posibles de X
    !call genDatosGraf(40)
    
    
    if (pause_test) pause

    close (unit=333)
    print *, test_ok
    999 continue
end program check

subroutine leer_input_flash_ALL7()
    use InputData
    use flash
    use CUFAC 
    implicit none    
    integer::i,j,k
    
    
    call open_file_name_ALL7()
    
    OPEN (UNIT=2,FILE=name,status='OLD',FORM='FORMATTED')
    READ(2,501) NTEXT   
    501 FORMAT(36A2)  
    
    READ(2,*) ICALC,modelo,IPRm,IOUTm,NOVAPm,igm, ipareq     
    
    call ab_ban1(modelo)
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
endsubroutine leer_input_flash_ALL7

subroutine open_file_name_ALL7()
    use InputData
    implicit none
    
    integer::parameters
    
    !apertura de bases de datos generadas desde excel
    OPEN (UNIT=1,FILE='test/ALL7-name.dat',status='OLD',FORM='FORMATTED')
    read(1,*)parameters 
    read(1,"(A36)") name
    name = name(2:len_trim(name)-1)
    name = "test/ALL7-"//name
    !name = "test/llecalas2.dat"
    if (parameters==1)then
        call LeerBases()
        stop
    endif    
    CLOSE (UNIT=1)  

endsubroutine open_file_name_ALL7