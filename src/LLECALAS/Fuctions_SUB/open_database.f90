subroutine open_database(model)
    !---------------------------------------------------------------------------
    !   This subroutine opens each database according to UNIFAC, UNIQUAC,
    !   A-UNIFAC and Group-Contribution (GC) models:
    !
    !   UNIQUAC & UNIFAC:
    !   - INTRCN.MDS        (UNIT=13) Interaction parameters.
    !   - GRUPOSRAM.MDS     (UNIT=14) Molecular groups data.
    !
    !   A-UNIFAC:
    !   - INTRCNAS.MDS      (UNIT=13) Associative interaction parameters.
    !   - GRUPOSRAM.MDS     (UNIT=14) Molecular groups data.
    !   - PARVOLAS.MDS      (UNIT=15) Associative volume  parameters.
    !   - PARENEAS.MDS      (UNIT=16) Associative energy  parameters.
    !
    !   GC:
    !   - INTRCNGCALPHA.MDS (UNIT=13) Alpha interaction parameters for GC.
    !   - GRUPOSRAMGC.MDS   (UNIT=14) Molecular groups data for GC.
    !   - INTRCNGCKAPA.MDS  (UNIT=16) Kappa interaction  parameters for GC.
    !
    !   The OPEN statement can connect an existing external file to a unit, 
    !   create a file and connect it to a unit, or change some specifiers of the
    !   connection.The following provides a detailed list of the OPEN specifier 
    !   keywords:
    !       - UNIT=u: 
    !           u is an interger expression that specifies the unit number where
    !           the external file is connected.
    !       - FILE=fin:
    !           fin is a character string containing the name and path of the 
    !           file being opened.     
    !       - ACCESS=acc: 
    !           The "ACCESS=acc" clause is optional. "acc" is a character 
    !           expression. If ACCESS='DIRECT': RECL must also be given, 
    !           since all I/O transfers are done in multiples of fixed-size 
    !           records. A direct-access file contains a number of records that 
    !           are written to or read from by referring to the record number. 
    !           Direct access is also called random access.
    !           If FORM is not given, formatted transfer is assumed. 
    !       - FORM=fm:
    !           The FORM=fm clause is optional. fm is a character expression.
    !           The default is 'FORMATTED'.
    !           If FORM='FORMATTED', each record is terminated with a newline 
    !           (\n) character; that is, each record actually has one extra 
    !           character.
    !       - RECL=rl:
    !           The RECL=rl clause is required if ACCESS='DIRECT' and 
    !           ignored otherwise. It indicates the length of each record in a 
    !           file connected for direct access, or the maximum length of a 
    !           record in a file connected for sequential access. (recl=850)
    !       - ERR=s
    !           The ERR=s clause is optional. s is a statement label of a 
    !           statement to branch to if an error occurs during execution of 
    !           the OPEN statement.
    !       - IOSTAT=ios:
    !           The IOSTAT=ios clause is optional. ios is an integer variable 
    !           that receives the error status from an OPEN. After the execution
    !           of the OPEN, if no error condition exists, then ios is zero; 
    !           otherwise, it is some positive number.
    !           If you want to avoid aborting the program when an error occurs  
    !           on an OPEN, include ERR=s or IOSTAT=ios.
    !       - ACTION=act:
    !           This specifier denotes file permissions. Possible values are: 
    !           READ, WRITE, and READWRITE.
    !           If act is READ, it specifies that the file is opened read-only.
    !           If act is WRITE, it specifies that the file is opened 
    !           write-only. You cannot execute a BACKSPACE statement on a 
    !           write-only file.
    !           If act is READWRITE, it specifies that the file is opened with 
    !           both read and write permissions.
    !       - STATUS=sta:
    !           The STATUS=sta clause is optional. 
    !           sta is a character expression. Possible values are: 
    !           'OLD', 'NEW', 'UNKNOWN', or 'SCRATCH'.
    !           If sta is 'OLD', the file already exists 
    !           (nonexistence is an error).
    !---------------------------------------------------------------------------

    use iso_fortran_env, only: int8, int16

    implicit none
    
    ! The model which is currently used:
    integer(kind=int8), intent(in) :: model
    ! An error-handling variable related to OPEN:
    integer(kind=int16) :: stat

    ! These are the PATH and files to be opened:
    character(len=*), parameter :: path = "src/database/"
    character(len=*), parameter :: intrcn_mds = "intrcn.mds"
    character(len=*), parameter :: intrcnas_mds = "intrcnas.mds"
    character(len=*), parameter :: parvolas_mds = "parvolas.mds"
    character(len=*), parameter :: pareneas_mds = "pareneas.mds"
    character(len=*), parameter :: gruposram_mds = "gruposram.mds"
    character(len=*), parameter :: intrcngcalpha_mds = "intrcngcalpha.mds"
    character(len=*), parameter :: gruposramgc_mds = "gruposramgc.mds"
    character(len=*), parameter :: intrcngckapa_mds = "intrcngckapa.mds"

    if(model /= 3) then ! When model is not GC (3)...    
        !if(model == 1) then ! When model is UNIQUAC (1)... 
        ! [The correct way to proceed was to use only A-UNIFAC for the 
        ! associative interaction parameters, as classic UNIFAC is not 
        ! suitable for this purpose. The previous statement was erroneous and 
        ! should be disregarded.]
        if(model /= 2) then ! When model is either UNIFAC (0) or UNIQUAC (1)...
            open (unit=13, file=path//intrcn_mds, status='old',&
            access='direct', form='formatted', recl=850, action='read',&
            iostat=stat)
            if (stat /= 0) then ! check for errors
                print *, 'Error opening file ', intrcn_mds
                stop
            end if   
        else ! When model is A-UNIFAC (2)...    
            open (unit=13, file=path//intrcnas_mds, status='old',&
            access='direct',form='formatted',recl=850, action='read', &
            iostat=stat)
            if (stat /= 0) then ! check for errors
                print *, 'Error opening file ', intrcnas_mds
                stop
            end if   
            open (unit=15, file=path//parvolas_mds, status='old', &
            access='direct', form='formatted', recl=850, action='read', &
            iostat=stat)
            if (stat /= 0) then ! check for errors
                print *, 'Error opening file ', parvolas_mds
                stop
            end if   
            open (unit=16, file=path//pareneas_mds, status='old', &
            access='direct', form='formatted', recl=850, action='read', &
            iostat=stat)
            if (stat /= 0) then ! check for errors
                print *, 'Error opening file ', pareneas_mds
                stop
            end if          
        endif

        ! When model is either A-UNIFAC (0), UNIQUAC (1) or A-UNIFAC (2)...
        open (unit=14, file=path//gruposram_mds, status='old', &
            access='direct',form='formatted',recl=300, action='read', &
            iostat=stat)
        if (stat /= 0) then ! check for errors
            print *, 'Error opening file ', gruposram_mds
            stop
        end if  
        !open (unit=15, file='src/database/parvolas.mds', status='old', &
        !   access='direct', form='formatted', recl=850)
        !open (unit=16, file='src/database/pareneas.mds', status='old', &
        !   access='direct', form='formatted', recl=850)
        ! [The four previous lines were misplaced in the code, as they are only 
        ! relevant for the A-UNIFAC model (model = 2), which uses associative 
        ! interaction parameters. For the classic UNIFAC model (model = 0), 
        ! these parameters are not required. The lines have been moved to the 
        ! appropriate section of the code.]
           

    else ! When model is GC (3)...
        open (unit=13, file=path//intrcngcalpha_mds, status='old',&
        access='direct', form='formatted', recl=730, action='read', &
        iostat=stat)
        if (stat /= 0) then ! check for errors
            print *, 'Error opening file ', intrcngcalpha_mds
            stop
        end if  
        open (unit=14, file=path//gruposramgc_mds, status='old',&
            access='direct', form='formatted', recl=263, action='read', &
            iostat=stat)
        if (stat /= 0) then ! check for errors
            print *, 'Error opening file ', gruposramgc_mds
            stop
        end if  
        open (unit=16, file=path//intrcngckapa_mds, status='old',&
            access='direct', form='formatted', recl=730, action='read', &
            iostat=stat)
        if (stat /= 0) then ! check for errors
            print *, 'Error opening file ', intrcngckapa_mds
            stop
        end if      
    endif
    
    ! Since "ab_ban" is not a function, it does not need to return.
    !return 
    
    end subroutine open_database