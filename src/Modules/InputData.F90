module InputData
    
    use iso_fortran_env
    
    implicit none
    
    ! Max site number per group.
    integer(kind=int32),parameter :: NMS = 2
    ! Max size of the arrays which contains subgroups' data .
    integer(kind=int32),parameter :: NMG = 150
    ! Max size of the arrays which contains the interaction parameters.
    integer(kind=int32),parameter :: NINT = 70
    ! The numbers of the second line of the flash-parameters-file.
    !   icalc:  0-' **** FLASH CALCULATION ****'
    !           1-' **** BINODAL CURVE CALCULATION ****'
    !           2-' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** '
    !   model:  0-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'
    !           1-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'
    !           2-' MODEL USED FOR LIQUID PHASE NON-IDEALITY: A-UNIFAC'
    !   iprm:    1-' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND 
    !                UNIQUAC, RESPECTIVELY **'
    !   ioutm:   1-'open 'lleasoccuzada.OUT''
    !   novapm:  0-'VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS'
    !   igm:     0-'write compositions'
    !           1-'write compositions and activities'
    !   ipareq: 1-'liquid-liquid parameters table (UNIFAC)'
    !           2-'liquid-vapor parameters table (UNIFAC)'
    !           3-'infinite dilution parameters table (UNIFAC)'
    !           4-'GC-EOS parameters'
    integer(kind=int32) :: ipareq,icalc,modelo,iprm,ioutm,novapm,igm

    integer(kind=int32) :: output
    ! The name of the file with the filename of the flash-parameters-file.
    !character(len=8) :: name_filename = "name.dat"
    ! The max character lenght of the filename with the flash parameters.
    integer(kind=int32), parameter :: name_maxlen = 36
    ! The filename of the flash-parameters-file.
    character(len=name_maxlen) :: name
    ! A matrix with the size of (the max site number per group x
    ! Max size of the arrays which contains subgroups' data)
    real(kind=real64), dimension(NMG,NMG) :: aint1

    ! A string with the first line of the flash calc input file. It is usually
    ! the same name in the namefile
    integer(kind=int32), dimension(name_maxlen) :: NTEXT
    
    real(kind=real64), dimension(10,3):: ANT
    
endmodule InputData