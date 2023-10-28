subroutine open_file(filename,array,max_lines,max_chars)
    ! A subroutine that reads each line of a file and saves its values in an 
    ! array. There is a maximum number of lines in the file and
    ! maximum number of characters per line allowed to read.

    implicit none
    
    ! Maximum number of lines in the file (100 suggested)
    integer, intent(in)  :: max_lines 
    ! Maximum number of characters per line (256 suggested)
    integer, intent(in):: max_chars 
    
    ! Array to store the lines
    character(len=max_chars), dimension(max_lines):: array 
    
    integer :: i, n, stat
    integer :: file_unit = 1
    
    ! Name of the file to open
    character(len=*),intent(in) :: filename
    
    open(unit=file_unit, file=filename,&
        status='old', action='read', iostat=stat) ! open the file for reading
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', filename
        stop
    end if
    
    n = 0 ! number of lines read
    do i = 1, max_lines ! loop over the lines
        read(file_unit, '(A)', iostat=stat) array(i) ! read a line into the array
        if (stat == -1) exit ! end of file reached
        if (stat /= 0) then ! check for errors
        !print *, 'Error reading file ', filename
        stop
        end if
        n = n + 1 ! increment the number of lines read
    end do
    
    close(file_unit) ! close the file

end subroutine open_file