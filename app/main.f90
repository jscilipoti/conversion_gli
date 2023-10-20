program main 
    
    use flash 
    use fobjtype
    use CUFAC
    
    implicit none
    
    integer::i
    integer::variables
    real*8::fmin
    real*8,allocatable,dimension(:)::x
    
    double precision,external::conversion
    double precision,external::praxis_n,f_n ,newton 
    
    OPEN (unit=333,file='salida.OUT')
    call leer_input_flash()

    
    variables = 1
    allocate(x(variables))
    
    x(1) = 0.0000001
    
    !Ve en el caso de la extraccion del agua    
    if (.false.) then
        agextr = .true.
        aglim = 1.D-4
    else
        agextr = .false.
    endif
    
    !Optimiza la funcion f_n en el documento Fobj
    !FMIN = praxis_n ( 3.D-4, 1.D1, variables, 3, x, f_n )
    
    !Genera una lista de concentraciones segun la funcion genDat que esta en el documento Conversion_f     
    !call genDat(x,variables)
    call genDatExtr(x,variables,2,10,.false.)
    
    !Genera una espacion de datos de la funcion objetivo en todos los valores posibles de X
    !call genDatosGraf(40)
    
    
    !pause

    close (unit=333)
endprogram main
