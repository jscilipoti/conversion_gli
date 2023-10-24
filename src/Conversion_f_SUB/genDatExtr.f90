subroutine genDatExtr(x,n,extI,extS,temp)
    use Flash
    use flashout

    implicit none
    logical::temp,proc
    integer::extI,extS,camb
    integer,intent(in)::n
    real*8::x(n)
    real*8,dimension(size(z))::zf
    doubleprecision::obj
    
    integer::i
    
    i = extI
    camb = 1
    proc = .true.
    !do while ((proc) .and. (i <= extS) .and. (i >= extI))
        
        write(*,*)
        write(*,*)
        write(*,*) i
        write(*,*)
        write(*,*)
        
        x(1) = i
        !T = x(1) + 273.15
    
    !    x(1) = (i - 1)/50. + 0.000000001 !+ 393.15
    !    z(1) = x(1)
    !    z(2) = (101 - i)/50. + 0.000000001
        
        contador=0
        call conversion(zf)
        
   !     if (isnan(zf(1))) then
   !         if (camb == 1) then
   !             camb = -1
   !             i = extS
   !         else
   !             proc = .false.
   !         endif
   !     else
            !Escribe en el documento de salida las condiciones
        open(unit=3333,file="salida3.out",status="old")
            write(3333,*) i,";",zf(1),";",zf(2),";",zf(3),";",zf(4),";",zf(5),";",zf(6),";",zf(7)
            write(*,*) i,zf(1),zf(2),zf(3),zf(4),zf(5),zf(6),zf(7)
            write(*,*) "Acido: ",zf(1)
            write(*,*) "Gly: ", zf(2)
            write(*,*) "Mono: ", zf(3)
            write(*,*) "Di: ", zf(4)
            write(*,*) "Tri: ", zf(5)
            write(*,*) "Agua: ", zf(6)
            write(*,*) "Tol: ", zf(7)
            close(unit=3333)
            i = i + camb
   !     endif
 !   enddo
    write(*,*) "K1", ksalida(1)
    write(*,*) "K2", ksalida(2)
    write(*,*) "K3", ksalida(3)

endsubroutine genDatExtr