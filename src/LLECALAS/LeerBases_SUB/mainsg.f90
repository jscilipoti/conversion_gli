!integer function mainsg (i1,ipareq)
!!c-------------------------------------------------------------------
!!c      Esta funcion devuelve el numero de grupo principal correspon-
!!c      diente al subgrupo i1 de la tabla de parametros ipareq:
!!c                     1: liquido-liquido
!!c                     2: liquido-vapor
!!c                     3: dilucion-infinita
!!c-------------------------------------------------------------------
!!c
!    integer:: il,ipareq
!	irec1=i1+(ipareq-1)*150
!      read (14,10,rec=irec1) main
!      mainsg = main
!  10  format (i4)
!      return
!endfunction mainsg