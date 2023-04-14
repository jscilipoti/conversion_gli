doubleprecision function f(x,n)
    use Flash
    use flashout

    implicit none
    integer,intent(in)::n
    real*8::x(n)
    integer::i
    real*8::act(10),Keqcalc(3),Kexp(3),zf(10),esq,errneg
    double precision,external::conversion
    logical::esneg
    !Sentencias


        !Calculo composici�n final a partir de la nueva conversi�n
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
        
        !Analizo si los valores son positivos
        esneg = .false.
        errneg = 1
        do i=1,size(zf)
            if (zf(i) < 0) then
                esneg = .true.
                errneg = errneg*zf(i)*(-1)
            endif
        enddo
        if (esneg) then
            f = errneg*10D50
            contador = contador + 1
            return
        endif
        
        !normalizo Z
        zf(:) = zf(:)/sum(zf(:))
        
        !flash
        call llecalas(T,P,zf)
        
        !c�lculo de actividades
        do i=1,size(z)
            if (compfases(i,1) > 1.0D-10) then
                act(i) = compfases(i,1)*exp(agam(i,1))
            else
                act(i) = compfases(i,2)*exp(agam(i,2))
            endif
        enddo
        
        !c�lculo de Keq(clac)
        Keqcalc(1) = act(3)*act(6)/(act(2)*act(1))
        Keqcalc(2) = act(4)*act(6)/(act(3)*act(1))
        Keqcalc(3) = act(5)*act(6)/(act(4)*act(1))
        
        !Datos de Kexp        
        
        !Acido Oleico 1 (T = 413.15-453.15 K) da Silva

        !Kexp(1) = 0.000124912*exp((-8296.4)*((1/T)-(1/298.15)))
        !Kexp(2) = 2.43536E-07*exp((-15050)*((1/T)-(1/298.15)))
        !Kexp(3) = 0.001392994*exp((-5976.2)*((1/T)-(1/298.15)))
        
        !Acido Oleico 2 (T = 393.15-433.15 K) Kong
        
        !Kexp(1) = 2.25306E-05*exp((-11250.90807)*((1/T)-(1/298.15)))
        !Kexp(2) = 1.90293E-05*exp((-11817.11444)*((1/T)-(1/298.15)))
        !Kexp(3) = 0.000145242*exp((-8818.098831)*((1/T)-(1/298.15)))
        
        !Acido Oleico 3 (T = 413.15-473.15 K) Kotwel
        
        !Kexp(1) = 1.71818E-09*exp((-17855)*((1/T)-(1/298.15)))
        !Kexp(2) = 7.88412E-10*exp((-18325)*((1/T)-(1/298.15)))
        !Kexp(3) = 8.23053E-10*exp((-17246)*((1/T)-(1/298.15)))
        
        !Acido Oleico 4 (T = 443.15-483.15 K) Z. Zhang
        
        !Kexp(1) = 0.003923387*exp((-5230.7)*((1/T)-(1/298.15)))
        !Kexp(2) = 0.001259551*exp((-7017.1)*((1/T)-(1/298.15)))
        !Kexp(3) = 6.08916E-06*exp((-10229)*((1/T)-(1/298.15)))
        
        !Acido Oleico 5 (T = 393.15-453.15 K) Ratchadapiban
        
        !Kexp(1) = 0.000511379*exp((-7416)*((1/T)-(1/298.15)))
        !Kexp(2) = 0.000178236*exp((-8902.5)*((1/T)-(1/298.15)))
        !Kexp(3) = 2.98947E-07*exp((-13842)*((1/T)-(1/298.15)))
        
        !Acido Oleico 6 (T = 393.15-483.15 K) Promedio

        Kexp(1) = 4.30465!0.002157910*exp((-5861.4)*((1/T)-(1/298.15)))
        Kexp(2) = 1.2!0.000339003*exp((-8253.8)*((1/T)-(1/298.15)))
        Kexp(3) = 0.239!0.000206217*exp((-7568.0)*((1/T)-(1/298.15)))
        
        f = 0.
        do i=1,size(Keqcalc)
            !f = f + ((Keqcalc(i) - Kexp(i))/Kexp(i))**2
            if (Keqcalc(i) > Kexp(i)) then
                f = f + (1 - (kexp(i)/keqcalc(i)))**2
            else
                f = f + (1 - (keqcalc(i)/kexp(i)))**2
            endif
        enddo
        contador = contador + 1
        
    endfunction f
    
    
    
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
    
doubleprecision function obj(zf)
    use Flash
    use flashout
    
    implicit none
    real*8,dimension(size(z))::zf
    double precision::totz
    
    obj = 0.
    
    !Variables que se buscan optimizar
    
    !Maximizar la cantidad de MG
    !obj = obj + (1 - zf(3)/sum(zf(:)))**2
    
    !Maximizar la cantidad de DG
    !obj = obj + (1 - zf(4)/sum(zf(:)))**2
    obj = obj + (1-zf(4)/(zf(3)+zf(4)+zf(5)))**2
    
    !Maximizar la cantidad de TG
    !obj = obj + (1 - zf(5)/sum(zf(:)))**2
    
endfunction
    
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