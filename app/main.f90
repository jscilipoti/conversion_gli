program main 
    
    use flash 
    use fobjtype
    
    implicit none
    
    integer::i
    integer::variables
    real(kind=8) ::fmin
    real(kind=8) ,allocatable,dimension(:)::x
    
    double precision,external::conversion
    double precision,external::praxis_n,f_n ,newton 
    
    OPEN (unit=333,file='salida.OUT')
    call leer_input_flash()

    
    variables = 1
    allocate(x(variables))
    
    x(1) = 0.0000001
    
    if (.false.) then
        agextr = .true.
        aglim = 1.D-4
    else
        agextr = .false.
    endif
    
    !Genera una lista de concentraciones segun la funcion genDat que esta en el documento Conversion_f     
    call genDatExtr(x,variables,2,10,.false.)
    
    
    pause

    close (unit=333)
endprogram main
    
subroutine leer_input_flash()
    use InputData
    use flash 
    implicit none   
   
    integer :: N, i, j, k, ng
    real(kind=8) :: Tx, px
    COMMON/CUFAC/N,NG,Px(10,10),Tx

    call open_file_name()
    open(UNIT=2, file=name, status='OLD', form='FORMATTED')
    read(2, '(36A2)') NTEXT
    read(2, *) ICALC, modelo, IPRm, IOUTm, NOVAPm, igm, ipareq
    call ab_ban1(modelo)
    CALL PARIN2
    IF(NOVAPm/=0) then                                            
        DO J=1,N                                                                                          
            READ(2,*) (ANT(K,J),K=1,3)                                       
        ENDDO 
        DO J=1,N                                                       
            ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                            
            ANT(2,J)=2.302585*ANT(2,J)                                       
        ENDDO 
    endif  
    READ(2,502) T,P
    502 FORMAT(2F10.2) 
    READ(2,*) (Z(I),I=1,N)

end subroutine leer_input_flash
    
    
subroutine open_file_name()
    use InputData
    implicit none
    
    integer::parameters
    
    ! Open file 'name.dat'
    OPEN (UNIT=1,FILE='name.dat',status='OLD',FORM='FORMATTED')
    
    ! Read parameters from file
    read(1,*)parameters 
    
    ! Read name from file
    read(1,"(A36)") name
    
    ! Remove extra characters from name
    name = name(2:len_trim(name)-1)
    
    ! Check if parameters is 1
    if (parameters==1)then
        call LeerBases()
        stop
    endif   
    
    ! Close file
    CLOSE (UNIT=1) 

end subroutine open_file_name