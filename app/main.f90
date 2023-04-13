program main 
    
    use flash 
    use fobjtype
    
    implicit none
    
    integer::i
    integer::variables
    real*8::fmin
    real*8,allocatable,dimension(:)::x
    
    double precision,external::conversion
    double precision,external::praxis_n,f_n ,newton 
    
    OPEN (unit=333,file='salida.OUT')
    call leer_input_flash()

    
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
    call genDatExtr(x,variables,2,98,.false.)
    
    !Genera una espacion de datos de la funcion objetivo en todos los valores posibles de X
    !call genDatosGraf(40)
    
    
    pause

    close (unit=333)
endprogram main
    
subroutine leer_input_flash()
    use InputData
    use flash 
    implicit none    
    integer::N,i,j,k,ng
    real*8::Tx,px
    COMMON/CUFAC/N,NG,Px(10,10),Tx

    
    
    call open_file_name()
    OPEN (UNIT=2,FILE=name,status='OLD',FORM='FORMATTED')
    READ(2,501) NTEXT   
    501 FORMAT(36A2)  
    READ(2,*) ICALC,modelo,IPRm,IOUTm,NOVAPm,igm, ipareq     
    call ab_ban1(modelo)
    CALL PARIN2
    IF(NOVAPm/=0) then                                             
        DO 6 J=1,N                                                        
!C   6 READ(2,502) (ANT(K,J),K=1,3)                                      
    6       READ(2,*) (ANT(K,J),K=1,3)                                        
        DO 7 J=1,N                                                        
            ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
    7       ANT(2,J)=2.302585*ANT(2,J)                                        
    endif   
    READ(2,502) T,P
    502 FORMAT(2F10.2)  
    READ(2,*) (Z(I),I=1,N)  
    close(unit=2)
endsubroutine leer_input_flash
    
    
subroutine open_file_name()
    use InputData
    implicit none
    
    integer::parameters
    
    !apertura de bases de datos generadas desde excel
    OPEN (UNIT=1,FILE='name.dat',status='OLD',FORM='FORMATTED')
    read(1,*)parameters 
    read(1,"(A36)") name
    name = name(2:len_trim(name)-1)
    if (parameters==1)then
        call LeerBases()
        stop
    endif    
    CLOSE (UNIT=1)  

endsubroutine open_file_name