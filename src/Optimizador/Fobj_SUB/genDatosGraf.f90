subroutine genDatosGraf(valLim)
    use InputData
    use Flash
    use fobjtype
    
    implicit none
    
    integer::valLim
    real*8,allocatable,dimension(:)::x
    integer::i,j,k
    double precision::f,limi1,limi2,limi3,lims1,lims2,lims3
    
    allocate(x(3))
    
    !Limites de los valores de X para la generacion del espacio de datos Fobj
    limi1 = 0.07
    limi2 = 0.07
    limi3 = 0.34
    lims1 = 0.13
    lims2 = 0.13
    lims3 = 0.5

    do i=0,valLim
        do j=0,valLim
            do k=0,valLim
                x(1) = limi1+(i*(lims1-limi1))/valLim
                x(2) = limi2+(j*(lims2-limi2))/valLim
                x(3) = limi3+(k*(lims3-limi3))/valLim
                write(333,*) x,f(x,3)
            enddo
        enddo
    enddo


endsubroutine