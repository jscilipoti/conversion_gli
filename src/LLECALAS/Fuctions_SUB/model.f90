SUBROUTINE Model(m)
    !  use Input
        COMMON/AS/ASOC
        LOGICAL ASOC
        INTEGER M
        
        IDEV=6
        ASOC=.FALSE.
        write (idev,21)  ! Choose the MODEL
  61    write (idev,30) !"> "
        READ (5,981,ERR=61) M
        iF (M<1.OR.M>3) GOTO 61	
        IF(M.EQ.2)ASOC=.TRUE.
    !  InputProblem%model = M
   21   format (//' Choose the model:',//,&
               15x,'UNIFAC  : 1',/,&
               15x,'A-UNIFAC: 2',/,&
               15x,'GC-EOS  : 3')
   30   format (1x,/,60x,'> ',$)
   981  FORMAT (I1)
  endsubroutine