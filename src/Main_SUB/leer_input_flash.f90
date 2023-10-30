subroutine leer_input_flash(name_filename)
    ! This Fortran subroutine opens a file and reads its contents. 
    ! It first reads the number of parameters and then reads the name of the 
    ! file which contains the data for the flash calculation. 
    ! The subroutine then removes the leading and trailing quotes from the name.
    ! If the number of parameters is 1, the subroutine reads the bases from the 
    ! file using the `LeerBases()` subroutine and stops the program. Finally, 
    ! the file is closed.
    use InputData
    use flash
    use CUFAC 
    implicit none    
    integer::i,j,k

    integer :: parameters = 0
    character(len=*), intent(in) :: name_filename
    character(len=name_maxlen), dimension(2) :: file_data
    
    
    ! Open the file for reading
    call open_textfile(name_filename,file_data,2,name_maxlen)

     ! Read the number of parameters and the name from the file
    parameters = ichar(trim(file_data(1)))
    name = file_data(2)

    ! Remove the leading and trailing quotes from the name
    name = name(2:len_trim(name)-1)

    ! If the number of parameters is 1, read the bases from the file
    if (parameters==1)then
        call LeerBases()
        stop
    endif
    
    Open (UNIT=2,FILE=name,status='OLD',FORM='FORMATTED')
    READ(2,501) NTEXT   
    501 FORMAT(36A2)  
    
    READ(2,*) ICALC,modelo,IPRm,IOUTm,NOVAPm,igm, ipareq     
    
    call ab_ban1(modelo)
    CALL PARIN2(NKK,NGG,Pxx,Txx)
    
    IF(NOVAPm/=0) then                                             
        DO 6 J=1,NKK                                                        
!C   6 READ(2,502) (ANT(K,J),K=1,3)                                      
    6       READ(2,*) (ANT(K,J),K=1,3)                                        
        DO 7 J=1,NKK                                                        
            ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
    7       ANT(2,J)=2.302585*ANT(2,J)                                        
    endif   
    READ(2,502) T,P
    502 FORMAT(2F10.2)  
    READ(2,*) (Z(I),I=1,NKK)  
    close(unit=2)
endsubroutine leer_input_flash
