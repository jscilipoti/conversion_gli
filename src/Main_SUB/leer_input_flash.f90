subroutine leer_input_flash(name_filename)
    ! This Fortran subroutine opens a file and reads its contents. 
    ! It first reads the number of parameters and then reads the name of the 
    ! file which contains the data for the flash calculation. 
    ! The subroutine then removes the leading and trailing quotes from the name.
    ! If the number of parameters is 1, the subroutine reads the bases from the 
    ! file using the `leerBases()` subroutine. Finally, 
    ! the file is closed.

    use CUFAC, only: NKK, NGG, Pxx, Txx
    use flash, only: P, T, Z
    use iso_fortran_env, only: int16, int8
    use InputData, only:&
        &flashInput_name, name_maxlen,&
        &ICALC, modelo, IPRm, IOUTm, NOVAPm, igm, ipareq,&
        &ANT,NTEXT
    use openunits, only: flashInput_unit
    
    
    implicit none    
    
    integer(kind=int16):: i, j, k, stat
    integer(kind=int8) :: doLeerBases = 0
    character(len=*), intent(in) :: name_filename
    character(len=name_maxlen), dimension(2) :: file_data
    
    ! Open the file for reading
    call open_textfile(name_filename, file_data, 2, name_maxlen)

     ! Read the number of parameters and the name from the file
    doLeerBases = ichar(trim(file_data(1)))
    flashInput_name = file_data(2)

    ! Remove the leading and trailing quotes from the name
    flashInput_name = flashInput_name(2:len_trim(flashInput_name)-1)

    ! If the number of doLeerBases is 1, read the bases from the file
    if (doLeerBases == 1) then
        call leerBases()
        stop
    end if
    
    !Old way to open flashInput file
    open(unit=2, file=flashInput_name, status='old', form='formatted',&
        &action='read', iostat=stat)
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', flashInput_name
        ERROR STOP 'Error opening file.'
    end if
    !New way to open flashInput file
    open(newunit=flashInput_unit, file=flashInput_name, status='old',&
        &form='formatted', action='read', iostat=stat)
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', flashInput_name
        ERROR STOP 'Error opening file.'
    end if

    ! Read the first line of the flash calc input file where the name is.
    read(2,*) NTEXT

    ! Read some info related to the calculation (see InputData module...)   
    read(2,*) ICALC, modelo, IPRm, IOUTm, NOVAPm, igm, ipareq     
    
    call open_database(modelo)
    call PARIN2(NKK, NGG, Pxx, Txx)
    
    ! The following code is run whether vapor phase is not included in flash
    ! calculations.
    if(NOVAPm /= 0) then                                             
        do j = 1, NKK                                                                                           
            read(2,*) (ANT(k,j), k = 1, 3)                                        
        end do
        do j = 1, NKK                                                        
            ANT(1,j) = 2.302585 * (ANT(1,j) - 2.880814)                             
            ANT(2,j) = 2.302585 * ANT(2,j)
        end do                                        
    end if

    ! Read the Temperature and the Pressure at the end of the flash input file
    read(2,'(2F10.2)') T, P
    
    ! Read the composition of each component of the system.
    read(2,*) (Z(i), i = 1, NKK) 
    
    close(flashInput_unit)
    close(unit=2)

end subroutine leer_input_flash
