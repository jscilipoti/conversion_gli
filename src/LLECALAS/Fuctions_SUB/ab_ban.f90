subroutine ab_ban1(model)
    !-------------------------------------------------------------------
    !      Esta subrunina ABRE los bancos de datos siguientes:
    !     - INTRCN.MDS      (UNIT=13)
    !     - INTRCNAS.MDS    (UNIT=13)
    !     - GRUPOSRAM.MDS   (UNIT=14)
    !     - PARVOLAS.MDS    (UNIT=15)  
    !     - PARENEAS.MDS    (UNIT=16)     
    !-------------------------------------------------------------------
    !
        !use input_data, only:model
          COMMON/AS/ASOC
          LOGICAL ASOC
          integer::model, mod
       
          
        !CALL MODEL(mod)
        if(model /= 3)then
            if(model==1)then
                open (unit=13,file='src/database/intrcn.mds',status='old',&
                      access='direct',form='formatted',recl=850)   
            else     
                open (unit=13,file='src/database/intrcnas.mds',status='old',&
                      access='direct',form='formatted',recl=850)        
            endif
            open (unit=14,file='src/database/gruposram.mds',status='old',&
                   access='direct',form='formatted',recl=300)
            open (unit=15,file='src/database/parvolas.mds',status='old',&
                  access='direct',form='formatted',recl=850)
            open (unit=16,file='src/database/pareneas.mds',status='old',&
                  access='direct',form='formatted',recl=850)    
        else
    
            open (unit=14,file='src/database/gruposramgc.mds',status='old',&
                   access='direct',form='formatted',recl=263)
            open (unit=13,file='src/database/intrcngcalpha.mds',status='old',&
                  access='direct',form='formatted',recl=730)
            open (unit=16,file='src/database/intrcngckapa.mds',status='old',&
                  access='direct',form='formatted',recl=730)    
            
        endif
          
    
          return
    endsubroutine ab_ban1