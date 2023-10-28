subroutine genDat(x,n)
    use Flash
    use flashout

    implicit none
    integer,intent(in)::n
    real*8::x(n)
    real*8,dimension(size(z))::zf
    doubleprecision::obj
    
    integer::i
    
    do i = 1,101
        
        x(1) = (i - 1)/50. + 0.000000001 !+ 393.15
        
        !verbose!write(*,*)
        !verbose!write(*,*)
        !verbose!write(*,*) i
        !verbose!write(*,*)
        !verbose!write(*,*)
    
        !z(1) = x(1)
        !aglim = x(1)
!        z(1) = x(1)
!        z(2) = (101 - i)/50. + 0.000000001
        
        contador=0
        call conversion(zf)
        
        !Escribe en el documento de salida las condiciones
        write(333,"(i0,6(es11.4))") i,zf(1),zf(2),zf(3),zf(4),zf(5),zf(6)
    
    enddo

endsubroutine