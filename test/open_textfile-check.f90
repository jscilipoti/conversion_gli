program check
    use do_tests
    implicit none
    character(len=8) :: filename = "name.dat"
    character(len=16) :: flash_filename = "llecalas2.dat"
    character(len=36),dimension(2) :: array = (/"",""/)
    print *,""
    print *, test_run//"open_file-check"
    if (pause_test) pause
    call open_textfile(filename,array,2,36)
    if (open_file_check) then
        if (filename/="name.dat")&
            ERROR STOP "Filename has been changed."
        if (size(array)/=2)&
            ERROR STOP "Max. lines has been changed."
        if (len(array(1))/=36)&
            ERROR STOP "Max. char number has been changed."
        if (ichar(trim(array(1)))*0/=0)&
            ERROR STOP "First line read was not a number"
        if (len_trim(array(2))/=15)&
            ERROR STOP "Second line has been changed(15)."
        if (trim(array(2))/='"llecalas2.dat"')&
            ERROR STOP "Second line has been changed."
    print *, test_ok
    else
        print *, test_disabled
    end if
end program check