doubleprecision function f_n(x,n)
    use Flash
    use flashout

    implicit none
    integer,intent(in)::n
    real*8::x(n)
    real*8,dimension(size(z))::zf
    doubleprecision::obj
    
    integer::i
    
    !Sentencias
    
    !Variables que se modifican para optimizar
    
    T = x(1)
    
    contador=0
    call conversion(zf)
    f_n = obj(zf)
    write(333,*) zf,f_n

endfunction f_n