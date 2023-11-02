program check
    use do_tests
    
    implicit none


    print *,""
    print *, test_run//"open_database-check"
    
    if (open_database_check) then
        if (pause_test) pause
        call open_database(0)
        close(unit=13)
        close(unit=14)
        call open_database(1)
        close(unit=13)
        close(unit=14)
        call open_database(2)
        close(unit=13)
        close(unit=14)
        close(unit=15)
        close(unit=16)
        call open_database(3)
        close(unit=13)
        close(unit=14)
        close(unit=16)
        !Check if everything went OK
        ! if ("string" /= "string")& 
        !     &ERROR STOP ""
        ! if (abs(0 - 0) > 1E-8)&
        !     &ERROR STOP ""
            
        print *, test_ok
    else 
        print *, test_disabled
    endif

!print *, "Put some tests in here!"

end program check
