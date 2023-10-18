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