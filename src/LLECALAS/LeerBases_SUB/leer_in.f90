real*8 function Leer_In(i1,i2,ipareq)
!c-------------------------------------------------------------------
!c      Esta funcion devuelve el parametro de interaccion entre los 
!c      grupos i1 y i2 de la tabla de parametros ipareq:
!c                     1: liquido-liquido
!c                     2: liquido-vapor
!c                     3: dilucion-infinita
!c-------------------------------------------------------------------
!c
      implicit real*8 (a-h,o-z)
      character*8 fs

      if (i1.eq.i2) then
          aint1 = 0.0
      else
	    irec4=i1+(ipareq-1)*70
	    read (13,20,rec=irec4) fs,(aint1,i=1,i2)
      end if
      
      Leer_In = aint1

  20  format (a8,70d12.5)
      return
end function Leer_In