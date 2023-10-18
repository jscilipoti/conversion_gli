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
