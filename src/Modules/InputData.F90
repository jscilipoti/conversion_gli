module InputData
    
    integer,parameter::NMS = 2 !número máximo de sitios por grupo
    integer,parameter::NMG = 150 !(dimensión máxima para vectores que guardan inf. sobre subgrupos)
    integer,parameter::NINT = 70 !Dimensión para vectores que guardan inf. sobre parámetros de interacción
    integer::ipareq,icalc,modelo,iprm,ioutm,novapm,igm
    integer::output
    character(len=36)::name
    real*8,dimension(NMG,NMG)::aint1
    integer:: NTEXT(36)
    real*8::ANT(10,3) 
    
endmodule InputData