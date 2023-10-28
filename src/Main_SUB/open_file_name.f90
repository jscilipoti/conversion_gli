subroutine open_file_name
    ! This Fortran subroutine opens a file named 'name.dat' and reads its 
    ! contents. It first reads the number of parameters and then reads the name 
    ! of the file. 
    ! The subroutine then removes the leading and trailing quotes from the name.
    ! If the number of parameters is 1, the subroutine reads the bases from the 
    ! file using the `LeerBases()` subroutine and stops the program. Finally, 
    ! the file is closed.
    use InputData
    implicit none
    
    integer::parameters = 0
    character(len=name_maxlen), dimension(2) :: file_data
    
    ! Open the file for reading
    call open_file(name_filename,file_data,2,name_maxlen)

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
end subroutine open_file_name


