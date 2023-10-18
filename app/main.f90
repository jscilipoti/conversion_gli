program main

    use flash
    use fobjtype
    !use iso_fortran_env, only: 8

    implicit none

    integer :: i
    integer :: variables
    real(8) :: fmin
    real(8), allocatable :: x(:)

    interface
        function conversion(x) result(y)
            real(8), intent(in) :: x
            real(8) :: y
        end function conversion

        function praxis_n(x, n) result(f)
            real(8), intent(in) :: x(n)
            integer, intent(in) :: n
            real(8) :: f
        end function praxis_n

        function f_n(x, n) result(f)
            real(8), intent(in) :: x(n)
            integer, intent(in) :: n
            real(8) :: f
        end function f_n

        function newton(x, n) result(f)
            real(8), intent(in) :: x(n)
            integer, intent(in) :: n
            real(8) :: f
        end function newton
    end interface

    open(333,file='salida.OUT')
    call leer_input_flash()

    variables = 1
    allocate(x(variables))

    x(1) = 0.0000001

    if (.false.) then
        agextr = .true.
        aglim = 1.0e-4
    else
        agextr = .false.
    endif

    !Genera una lista de concentraciones segun la funcion genDat que esta en el documento Conversion_f     
    call genDatExtr(x,variables,2,10,.false.)

    stop

    close(333)
end program main
    

    
    
subroutine open_file_name()
    use InputData
    implicit none
    
    integer::parameters
    
    ! Open file 'name.dat'
    OPEN (UNIT=1,FILE='name.dat',status='OLD',FORM='FORMATTED')
    
    ! Read parameters from file
    read(1,*)parameters 
    
    ! Read name from file
    read(1,"(A36)") name
    
    ! Remove extra characters from name
    name = name(2:len_trim(name)-1)
    
    ! Check if parameters is 1
    if (parameters==1)then
        call LeerBases()
        stop
    endif   
    
    ! Close file
    CLOSE (UNIT=1) 

end subroutine open_file_name