subroutine conversion(zf)
    use InputData
    use Flash
    use fobjtype
    implicit none ! Declara explícitamente todos los tipos de variables
    integer :: variables ! Número de variables
    real*8 :: fmin, esq, minmar ! Variables reales de doble precisión
    real*8, allocatable, dimension(:) :: x, xcic 
    ! Vectores reales de doble precisión asignables
    real*8, dimension(size(z)) :: zf 
    ! Vector real de doble precisión del mismo tamaño que z
    logical :: espos ! Variable lógica
    integer :: i, esp ! Variables enteras
    interface ! Define una interfaz para la función anexcep
        function anexcep(esp, minmar, x) 
            ! Esta función calcula el valor mínimo de una excepción
            use InputData ! Importa el módulo InputData
            use Flash     ! Importa el módulo Flash
            use fobjtype  ! Importa el módulo fobjtype
            implicit none ! Declara explícitamente todos los tipos de variables
            integer :: variables ! Número de variables
            real*8 :: fmin, esq, minmar, limI, limS, pas 
            ! Variables reales de doble precisión
            real*8, dimension(:) :: x ! Vector real de doble precisión
            real*8, allocatable, dimension(:) :: xcic 
            ! Vector real de doble precisión asignable
            real*8, dimension(size(z)) :: zf 
            ! Vector real de doble precisión del mismo tamaño que z
            logical :: espos ! Variable lógica
            integer :: i, esp ! Variables enteras
            real*8 :: anexcep ! Valor de retorno de la función
            integer, dimension(3,3) :: datx ! Matriz entera
            double precision, external :: praxis_n, f, newton 
            ! Funciones externas de doble precisión
        end function anexcep
    end interface
    
    double precision, external :: praxis_n, f, newton 
    ! Funciones externas de doble precisión
    
    variables = 3 ! Asigna el valor 3 a la variable variables
    allocate(x(variables), xcic(variables)) 
    ! Asigna memoria para los vectores x y xcic
    
    ! Asigna valores iniciales a los elementos de x basados en los 
    ! elementos de z
    x(1) = 0.1926d0 ! z(2)/(z(2) + z(3) + z(4) + z(5))
    x(2) = 0.2965d0 ! z(3)/(z(2) + z(3) + z(4) + z(5))
    x(3) = 0.4394d0 ! z(4)/(z(2) + z(3) + z(4) + z(5))
    
    fobj = .false. !
    
    write(*,*) x   ! Escribe el vector x en la salida estándar
    fmin = praxis_n(1.d-5, 5.d-2, variables, 0, x, f) 
    ! Llama a la función praxis_n para obtener 
    ! el valor mínimo de f con respecto a x
    write(*,*) fmin! Escribe el valor mínimo en la salida estándar
    
    esp = 23   ! Asigna el valor 23 a la variable esp
    minmar = 0.0001d0! Asigna el valor 0.0001 a la variable minmar
    
    if (fmin > minmar) then! Si el valor mínimo es mayor que el margen mínimo
        fmin = anexcep(esp,minmar,x)
        ! Llama a la función anexcep con el primer valor de esp para obtener 
        ! un nuevo valor mínimo
        if (fmin > minmar) then
            ! Si el nuevo valor mínimo sigue siendo mayor que el margen mínimo
            esp = 59   
            fmin = anexcep(esp,minmar,x)
            ! Llama a la función anexcep con el segundo valor de esp para 
            ! obtener otro valor mínimo
        end if
    end if
    
    if (fmin > minmar) then
        ! Si el valor mínimo final es mayor que el margen mínimo
        x(1) = 0d0! Asigna el valor 0 al primer elemento de x
        x(1) = 0d0/x(1)! Divide el primer elemento de x por sí mismo, lo que produce un error de división por cero
        x(2) = 0d0/x(1)! Divide el segundo elemento de x por el primero, lo que produce un error de división por cero
        x(3) = 0d0/x(1)! Divide el tercer elemento de x por el primero, lo que produce un error de división por cero
    end if
    
    call calcZ(x,variables,zf)
    ! Llama a la subrutina calcZ para calcular el vector zf 
    ! a partir de x y variables
end subroutine conversion
    
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
    
end subroutine calcZ
    
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
           
        contador=0
        call conversion(zf)
        
        !Escribe en el documento de salida las condiciones
        write(333,"(i0,6(es11.4))") i,zf(1),zf(2),zf(3),zf(4),zf(5),zf(6)
    
    enddo

end subroutine genDat
    
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
        
        write(*,*)
        write(*,*)
        write(*,*) i
        write(*,*)
        write(*,*)
        
        x(1) = i
        
        contador=0
        call conversion(zf)
        
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

    write(*,*) "K1", ksalida(1)
    write(*,*) "K2", ksalida(2)
    write(*,*) "K3", ksalida(3)

end subroutine genDatExtr
    
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

end function anexcep