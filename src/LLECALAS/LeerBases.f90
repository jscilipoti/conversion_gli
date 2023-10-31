subroutine leerBases()
    use InputData
    implicit none
    integer::model,  Ncomp, i, j, k,cant,parameters, irec, mgr,ICS,ITS1,ITS2,num,nga,ng,aa,bb
    integer,dimension(:,:,:):: ms(10,10,2)
    integer,dimension(NMG)::NPUNT,NGRUP,ngrupa,cs
    integer,dimension(NINT):: NPINT,NINTT,npunta
    integer,external::mainsgfunc
    integer::TS(NMG,NMS)    
    real*8::r,q,rrt,ENASST,RKASST
    logical::esta,vl

    open (unit=2,FILE=name,status='OLD',FORM='FORMATTED')
    open (unit=3,FILE="parameters.dat",FORM='FORMATTED')
    
    read(2,"(3i1)") model, ipareq, Ncomp
    call ab_ban1(model)

    ms(:,:,:)=0
    do i=1,Ncomp
        read(2,"(20i2)") (ms(I,J,1),ms(I,J,2),j=1,size(ms(1,:,1)))   
        j=1
        do while (ms(i,j,1)/=0)
            call cr_puntf (ms(i,j,1),NPUNT,NGRUP,NPINT,NINTT)
            j=j+1
        enddo 
    enddo
    
    call store_in(model,nintt)
 !Escritura de grupos   
    j=1
    do while (ngrup(j)/=0)
        write(3,"(i3)") ngrup(j)
        j=j+1
    enddo
    write(3,*)"end"
!Escritura de par�metros grupales
    j=1
    
    do while (ngrup(j)/=0)
        irec=ngrup(j)+(ipareq-1)*150
        read(14,"(42x,2d15.8)",rec=irec)r,q
        write(3,"(f10.3)") r
        write(3,"(f10.3)") q
        j=j+1
    enddo      
    write(3,*)"end"    
!Escritura par�metros interacci�n    
    I=1
    do while (ngrup(I)/=0)
       J=1
       do while (ngrup(J)/=0)

          write(3,"(f10.3)") Aint1(npint(mainsgfunc(ngrup(i),ipareq)),npint(mainsgfunc(ngrup(j),ipareq)))

          J=J+1
       enddo
       write(3,*)"end"
       I=I+1
    enddo
    write(3,*)"endint"   
!Escritura composici�n grupal
    do i=1,Ncomp
        k=1
        do while (ngrup(k)/=0)
            j=1
            esta=.False.
            do while (ms(i,j,1)/=0)
                if(ngrup(k)==ms(i,j,1))then
                    esta=.True.
                    exit
                endif    
                j=j+1
            enddo
            if(esta)then
                write(3,*)ms(i,j,2)
            else
                write(3,*)0
            endif              
            k=k+1
        enddo
        write(3,*)"end"
    enddo  
    
!Averigua cu�les grupos son asociativos    
    if (model==2) then
        NUM = 0
        NGA = 0
	    DO J=1,Ncomp
	        DO I=1,size(ms(1,:,1))
	            IF(MS(J,I,1).EQ.0)CYCLE
	            IREC=MS(J,I,1)+(ipareq-1)*150      
      	        READ(14,"(2x,I2,37X,D15.8,120X,3I2)",REC=IREC)MGR,RRT,ICS,ITS1,ITS2 !500  FORMAT(2x,I2,37X,D15.8,120X,3I2)
                IF (ICS.NE.0) THEN !GRUPOS ASOCIATIVOS
                    NG = MS(J,I,1)
                    CALL BUSCARAS (NG,NGRUPA,NMG,VL)
                    IF(.not.VL)THEN
                        NGA=NGA+1
                        NPUNTA(NG) = NGA
                        NGRUPA(NGA) = NG
                        CS(NGA) = ICS
                        TS(NGA,1) = ITS1
                        TS(NGA,2) = ITS2
                    ENDIF
                ENDIF      	    
            ENDDO
        ENDDO    
    

 !Escritura de grupos asociativos 
    j=1
    do while (ngrupa(j)/=0)
        write(3,"(i3)") ngrupa(j)
        write(3,"(i3)") cs(j)
        j=j+1
    enddo
    write(3,*)"end"
    
!      DO J=1,nga
!        DO K=1,nga
!            DO BB=1,NMS
!                DO AA=1,NMS
!                    if(AA>CS(k).or.bb>cs(j))then
!                        write(3,*)"-",AA,K,BB,J
!                        cycle
!                    endif
!                    IREC = MAINSGfunc(NGRUPA(K),IPAREQ)+(ipareq-1)*70
!                    IF(TS(J,BB).EQ.1.OR.TS(K,AA).EQ.1)THEN !Si alguno de ambos sitios es del tipo 1
!                       CALL LEEPAR (J,IREC,IPAREQ,NGRUPA,ENASST,RKASST)
!                       write(3,*)ENASST,AA,K,BB,J
!                       !RKASS(AA,K,BB,J)=RKASST                        
!                    ELSEIF (TS(J,BB).NE.TS(K,AA)) THEN
!                       CALL LEEPAR (J,IREC,IPAREQ,NGRUPA,ENASST,RKASST)
!                       write(3,*)ENASST,AA,K,BB,J
!!                       !RKASS(AA,K,BB,J)=RKASST
!                    ELSE !IF ((TS(J,BB).EQ.TS(K,AA)).and.(TS(J,BB)/=1).and.(TS(K,AA)/=1)) THEN
!                       write(3,*)0.0,AA,K,BB,J
!                       !RKASS(AA,K,BB,J)=0.0     
!                    ENDIF   
!                ENDDO
!            ENDDO
!        ENDDO
!      ENDDO  
      
      DO BB=1,NMS
      DO J=1,nga
        
            
            
                DO AA=1,NMS
                DO K=1,nga
                    if(AA>CS(k).or.bb>cs(j))then
                        write(3,*)"-" !,AA,K,BB,J
                        cycle
                    endif
                    IREC = MAINSGfunc(NGRUPA(K),IPAREQ)+(ipareq-1)*70
                    IF(TS(J,BB).EQ.1.OR.TS(K,AA).EQ.1)THEN !Si alguno de ambos sitios es del tipo 1
                       CALL LEEPAR (J,IREC,IPAREQ,NGRUPA,ENASST,RKASST)
                       write(3,*)ENASST !,AA,K,BB,J
                       !RKASS(AA,K,BB,J)=RKASST                        
                    ELSEIF (TS(J,BB).NE.TS(K,AA)) THEN
                       CALL LEEPAR (J,IREC,IPAREQ,NGRUPA,ENASST,RKASST)
                       write(3,*)ENASST !,AA,K,BB,J
!                       !RKASS(AA,K,BB,J)=RKASST
                    ELSE !IF ((TS(J,BB).EQ.TS(K,AA)).and.(TS(J,BB)/=1).and.(TS(K,AA)/=1)) THEN
                       write(3,*)0.0 !,AA,K,BB,J
                       !RKASS(AA,K,BB,J)=0.0     
                    ENDIF   
                ENDDO
            ENDDO
            write(3,*)"end"
        ENDDO
      ENDDO  
      write(3,*)"enden"         
    
       DO BB=1,NMS
      DO J=1,nga
        
            
            
                DO AA=1,NMS
                DO K=1,nga
                    if(AA>CS(k).or.bb>cs(j))then
                        write(3,*)"-" !,AA,K,BB,J
                        cycle
                    endif
                    IREC = MAINSGfunc(NGRUPA(K),IPAREQ)+(ipareq-1)*70
                    IF(TS(J,BB).EQ.1.OR.TS(K,AA).EQ.1)THEN !Si alguno de ambos sitios es del tipo 1
                       CALL LEEPAR (J,IREC,IPAREQ,NGRUPA,ENASST,RKASST)
                       write(3,*)RKASST !,AA,K,BB,J
                       !RKASS(AA,K,BB,J)=RKASST                        
                    ELSEIF (TS(J,BB).NE.TS(K,AA)) THEN
                       CALL LEEPAR (J,IREC,IPAREQ,NGRUPA,ENASST,RKASST)
                       write(3,*)RKASST !,AA,K,BB,J
!                       !RKASS(AA,K,BB,J)=RKASST
                    ELSE !IF ((TS(J,BB).EQ.TS(K,AA)).and.(TS(J,BB)/=1).and.(TS(K,AA)/=1)) THEN
                       write(3,*)0.0 !,AA,K,BB,J
                       !RKASS(AA,K,BB,J)=0.0     
                    ENDIF   
                ENDDO
            ENDDO
            write(3,*)"end"
        ENDDO
      ENDDO  
      write(3,*)"enden"      
    endif !if model == 2  

    close (unit=2)
    close (unit=3)
    pause
endsubroutine leerBases
