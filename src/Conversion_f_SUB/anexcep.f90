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
                !verbose!write(*,*) xcic
                !!verbose!write(*,*) xcic(3)/3 + xcic(2)/2 + xcic(1)/1
                do i=1,3
                    if (FMIN > minmar) then
                        x(1) = xcic(datx(i,1))/3.
                        x(2) = xcic(datx(i,2))/2.
                        x(3) = xcic(datx(i,3))/1.
                        if (x(1)+x(2)+x(3) < 1) then
                            FMIN = praxis_n(1.D-5,5.D-2,variables,0,x,F)
                            !verbose!write(*,*) FMIN
                            if (sum(xcic(:)) < 1.D-10) goto 10
                        endif
                    endif
                enddo
10              xcic(3) = xcic(3) + pas
                !verbose!write(*,*)
            enddo
            xcic(2) = xcic(2) + pas
        enddo
        xcic(1) = xcic(1) + pas
    enddo
    
    anexcep = FMIN

endfunction