subroutine calcZ(x,variables,zf)
   
    use Flash
    implicit none
    real*8::esq
    integer::variables
    real*8,dimension(variables)::x
    real*8,dimension(size(z))::zf
    
    esq = z(2) + z(3) + z(4) + z(5)
    zf(:) = z(:)
    zf(2) = esq * x(1)
    zf(3) = esq * x(2)
    zf(4) = esq * x(3)
    zf(5) = esq * (1 - sum(x(:)))
    zf(1) = z(1) + (z(3) - zf(3)) * 1 + (z(4) - zf(4)) * 2 + (z(5) - zf(5)) * 3
    if (agextr) then
        zf(6) = aglim
    else
        zf(6) = z(6) + (z(1) - zf(1))
    endif
    
endsubroutine calcZ