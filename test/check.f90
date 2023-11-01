program check
    use do_tests
    
    implicit none


    print *,""
    print *, test_run//"check"
    
    if (.true.) then
        if (pause_test) pause
 
        !Check if everything went OK
        if ("string" /= "string")& 
            ERROR STOP ""
        if (abs(0 - 0) > 1E-8)&
            ERROR STOP ""
            
        print *, test_ok
    else 
        print *, test_disabled
    endif

!print *, "Put some tests in here!"

end program check
