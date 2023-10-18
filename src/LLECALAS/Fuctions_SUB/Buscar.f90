SUBROUTINE BUSCAR (I,IPAREQL,FF,NUM)
    CHARACTER*8 FS,FF,MG
    LOGICAL IPAREQL
!C
    IPAREQL=.FALSE.
    num=(I-1)*150
    fs= ''
    MG= ''
    DO WHILE (MG.NE.'fin')
    num=num+1
    read (14,10,rec=num)MG,fs
    IF (FS.EQ.FF) THEN
        IPAREQL=.TRUE.
        EXIT
    ENDIF
  enddo
  !c
!c
  10  format (4x,2a8) 
      return
endsubroutine buscar