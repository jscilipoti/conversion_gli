subroutine open_textfile(filename,array,max_lines,max_chars)
    ! A subroutine that reads each line of a file and saves its values in an 
    ! array. There is a maximum number of lines in the file and
    ! maximum number of characters per line allowed to read.
    ! The OPEN statement can connect an existing external file to a unit, 
    ! create a file and connect it to a unit, or change some specifiers of the
    ! connection.The following provides a detailed list of the OPEN specifier 
    ! keywords:
    ! - UNIT=u: 
    !       u is an interger expression that specifies the unit number where
    !       the external file is connected.
    ! - FILE=fin:
    !       fin is a character string containing the name and path of the 
    !       file being opened.     
    ! - ACCESS=acc: 
    !       The "ACCESS=acc" clause is optional. "acc" is a character 
    !       expression. If ACCESS='DIRECT': RECL must also be given, 
    !       since all I/O transfers are done in multiples of fixed-size 
    !       records. A direct-access file contains a number of records that 
    !       are written to or read from by referring to the record number. 
    !       Direct access is also called random access.
    !       If FORM is not given, formatted transfer is assumed. 
    ! - FORM=fm:
    !       The FORM=fm clause is optional. fm is a character expression.
    !       The default is 'FORMATTED'.
    !       If FORM='FORMATTED', each record is terminated with a newline 
    !       (\n) character; that is, each record actually has one extra 
    !        character.
    ! - RECL=rl:
    !       The RECL=rl clause is required if ACCESS='DIRECT' and 
    !       ignored otherwise. It indicates the length of each record in a 
    !       file connected for direct access, or the maximum length of a 
    !       record in a file connected for sequential access. (recl=850)
    ! - ACTION=act:
    !       This specifier denotes file permissions. Possible values are: 
    !       READ, WRITE, and READWRITE.
    !       If act is READ, it specifies that the file is opened read-only.
    !       If act is WRITE, it specifies that the file is opened 
    !       write-only. You cannot execute a BACKSPACE statement on a 
    !       write-only file.
    !       If act is READWRITE, it specifies that the file is opened with 
    !       both read and write permissions.
    ! - ERR=s
    !   The ERR=s clause is optional. s is a statement label of a statement to 
    !   branch to if an error occurs during execution of the OPEN statement.
    ! - IOSTAT=ios:
    !       The IOSTAT=ios clause is optional. ios is an integer variable 
    !       that receives the error status from an OPEN. After the execution of 
    !       the OPEN, if no error condition exists, then ios is zero; otherwise,
    !       it is some positive number.
    !       If you want to avoid aborting the program when an error occurs on 
    !       an OPEN, include ERR=s or IOSTAT=ios.
    ! - STATUS=sta:
    !       The STATUS=sta clause is optional. 
    !       sta is a character expression. Possible values are: 
    !       'OLD', 'NEW', 'UNKNOWN', or 'SCRATCH'.
    !       If sta is 'OLD', the file already exists (nonexistence is an error).

    use iso_fortran_env, only: int16
    implicit none
    
    ! Maximum number of lines in the file (100 suggested)
    integer(kind=int16), intent(in)  :: max_lines 
    ! Maximum number of characters per line (256 suggested)
    integer(kind=int16), intent(in):: max_chars 
    
    ! Array to store the lines
    character(len=max_chars), dimension(max_lines):: array 
    
    integer(kind=int16) :: i, n, stat
    integer(kind=int16) :: file_unit
    
    ! Name of the file to open
    character(len=*),intent(in) :: filename
    
    open(unit=file_unit, file=filename,&
        &status='old', action='read', iostat=stat) ! open the file for reading
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', filename
        stop
    end if
    
    n = 0 ! number of lines read
    do i = 1, max_lines ! loop over the lines
        read(file_unit, '(A)', iostat=stat) array(i) ! read a line into the array
        if (stat == -1) exit ! end of file reached
        if (stat /= 0) then ! check for errors
        print *, 'Error reading file ', filename
        stop
        end if
        n = n + 1 ! increment the number of lines read
    end do
    
    close(file_unit) ! close the file

end subroutine open_textfile