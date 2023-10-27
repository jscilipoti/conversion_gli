module do_tests
    implicit none
    private
    public :: test_ok, test_disabled, test_run, test_error, ALL7_check, ULL7_check, ALV7_check, ULV7_check
  
 
    ! These are module variables that can be used by any program unit that uses this module
    character*14,parameter :: test_run = "Running test: "
    character*24,parameter :: test_ok = "Ok! "
    character*24,parameter :: test_disabled = "This test is disabled."
    character*24,parameter :: test_error = "ERROR! "
    logical,parameter:: ALL7_check = .false.
    logical,parameter:: ULL7_check = .false.
    logical,parameter:: ALV7_check = .true.
    logical,parameter:: ULV7_check = .true.

    
  
  end module do_tests