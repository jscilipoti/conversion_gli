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
        !verbose!write(*,*) compfases(1,1),compfases(2,1),compfases(3,1),compfases(4,1),&
        !compfases(5,1),compfases(6,1),compfases(7,1)
        !verbose!write(*,*) agam(1,1),agam(2,1),agam(3,1),agam(4,1),&
        !agam(5,1),agam(6,1),agam(7,1)
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

        Kexp(1) = 0.33274673!0.002157910*exp((-5861.4)*((1/T)-(1/298.15)))
        Kexp(2) = 1.373980!0.000339003*exp((-8253.8)*((1/T)-(1/298.15)))
        Kexp(3) = 0.249791!0.000206217*exp((-7568.0)*((1/T)-(1/298.15)))
        
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
        ksalida(:) = keqcalc(:)
        
    endfunction f