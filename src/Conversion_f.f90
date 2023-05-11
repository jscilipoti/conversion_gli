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
        
        write(*,*)
        write(*,*)
        write(*,*) i
        write(*,*)
        write(*,*)
    
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
    
    variables = 3
    allocate(xcic(variables))
    datx = reshape((/ 1, 2, 3, 2, 3, 1, 3, 1, 2 /), shape(datx))
    
    limI = (3*z(2)+2*z(3)+1*z(4)-z(1))/(z(2)+z(3)+z(4)+z(5))
    limS = (3*z(2)+2*z(3)+1*z(4)+z(6))/(z(2)+z(3)+z(4)+z(5))
    pas = (limS-limI)/esp

    xcic(1) = 0.
    
    FMIN = 3.
    
    do while ((FMIN > minmar) .and. (xcic(1) < limS) .and. (xcic(1)/3 < 1))
        xcic(2) = xcic(1)
        do while ((FMIN > minmar) .and. ((xcic(1) + xcic(2)) < limS) .and. ((xcic(2)/3 + xcic(1)/2) < 1))
            if ((limI - (xcic(1) + xcic(2))) > (xcic(2))) then
                xcic(3) = limI - (xcic(1) + xcic(2))
            else
                xcic(3) = xcic(2)
            endif
            do while ((FMIN > minmar) .and. (xcic(1) + xcic(2) + xcic(3) < limS) .and. ((xcic(3)/3. + xcic(2)/2. + xcic(1)/1.) < 1))
                write(*,*) xcic
                !write(*,*) xcic(3)/3 + xcic(2)/2 + xcic(1)/1
                do i=1,3
                    if (FMIN > minmar) then
                        x(1) = xcic(datx(i,1))/3.
                        x(2) = xcic(datx(i,2))/2.
                        x(3) = xcic(datx(i,3))/1.
                        if (x(1)+x(2)+x(3) < 1) then
                            FMIN = praxis_n(1.D-5,5.D-2,variables,0,x,F)
                            write(*,*) FMIN
                            if (sum(xcic(:)) < 1.D-10) goto 10
                        endif
                    endif
                enddo
10              xcic(3) = xcic(3) + pas
                write(*,*)
            enddo
            xcic(2) = xcic(2) + pas
        enddo
        xcic(1) = xcic(1) + pas
    enddo
    
    anexcep = FMIN

endfunction