module do_tests
    implicit none
    private
    public :: &
    pause_test,&
    test_ok,&
    test_disabled,&
    test_run,&
    test_error,&
    ALL7_check,&
    ULL7_check,&
    ALV7_check,&
    ULV7_check,&
    open_file_check,&
    leer_input_flash_check
  
 
    ! These are module variables that can be used by any program unit that uses this module
    character(len=*),parameter :: test_run = "Running test: "
    character(len=*),parameter :: test_ok = "Ok! "
    character(len=*),parameter :: test_disabled = "This test is disabled."
    character(len=*),parameter :: test_error = "ERROR! "
    
    logical,parameter :: pause_test = .false.

    logical,parameter :: ALL7_check = .false.
    logical,parameter :: ULL7_check = .false.
    logical,parameter :: ALV7_check = .true.
    logical,parameter :: ULV7_check = .true.
    logical,parameter :: open_file_check = .true.
    logical,parameter :: leer_input_flash_check = .true.

    
  
  end module do_tests