subroutine conversion(zf)
    use InputData
    use Flash
    use fobjtype
    implicit none
    integer::variables
    real*8::fmin,esq,minmar
    real*8,allocatable,dimension(:)::x
    real*8,allocatable,dimension(:)::xcic
    real*8,dimension(size(z))::zf
    logical::espos
    integer::i
    interface
    function anexcep(esp,minmar,x)

        use InputData
        use Flash
        use fobjtype
        implicit none
        integer::variables
        real*8::fmin,esq,minmar,limI,limS,pas
        real*8,dimension(:)::x
        real*8,allocatable,dimension(:)::xcic
        real*8,dimension(size(z))::zf
        logical::espos
        integer::i
        real*8::anexcep
        integer::esp
        integer, dimension(3,3)::datx
        
        double precision,external::praxis_n,f,newton
    endfunction
    endinterface
    integer::esp
    
    double precision,external::praxis_n,f,newton
    
    variables = 3
    allocate(x(variables))
    allocate(xcic(variables))
    x(1) = 0.1926 !z(2)/(z(2) + z(3) + z(4) +z(5))
    x(2) = 0.2965 !z(3)/(z(2) + z(3) + z(4) +z(5))
    x(3) = 0.4394 !z(4)/(z(2) + z(3) + z(4) +z(5))
    
    fobj = .false. !se llama desde conversion
    
    write(*,*) x
    FMIN = praxis_n(1.D-5,5.D-2,variables,0,x,F)
    write(*,*) FMIN
    write(*,*)
    
    esp = 23
    
    minmar = 0.0001
    
    if (FMIN > minmar) then
        FMIN = anexcep(esp,minmar,x)
        if (FMIN > minmar) then
            esp = 59
            FMIN = anexcep(esp,minmar,x)
        endif
    endif
    
    if (FMIN > minmar) then
        x(1) = 0
        x(1) = 0/x(1)
        x(2) = 0/x(1)
        x(3) = 0/x(1)
    endif
    
    call calcZ(x,variables,zf)
    
endsubroutine conversion
    

    

    

    
