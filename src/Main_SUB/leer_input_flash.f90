subroutine leer_input_flash()
    use InputData
    use flash 
    implicit none    
    integer::N,i,j,k,ng
    !real*8::Tx,px
    real*8::Tx
    real*8,dimension(10,10):: Px
    !COMMON/CUFAC/N,NG,Px(10,10),Tx
    COMMON/CUFAC/N,NG,Px,Tx

    
    
    call open_file_name()
    
    OPEN (UNIT=2,FILE=name,status='OLD',FORM='FORMATTED')
    READ(2,501) NTEXT   
    501 FORMAT(36A2)  
    
    READ(2,*) ICALC,modelo,IPRm,IOUTm,NOVAPm,igm, ipareq     
    
    call ab_ban1(modelo)
    CALL PARIN2(N,NG,Px,Tx)
    
    IF(NOVAPm/=0) then                                             
        DO 6 J=1,N                                                        
!C   6 READ(2,502) (ANT(K,J),K=1,3)                                      
    6       READ(2,*) (ANT(K,J),K=1,3)                                        
        DO 7 J=1,N                                                        
            ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
    7       ANT(2,J)=2.302585*ANT(2,J)                                        
    endif   
    READ(2,502) T,P
    502 FORMAT(2F10.2)  
    READ(2,*) (Z(I),I=1,N)  
    close(unit=2)
endsubroutine leer_input_flash
