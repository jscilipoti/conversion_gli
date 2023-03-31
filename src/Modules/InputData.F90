module InputData
    
    integer,parameter::NMS = 2 !n�mero m�ximo de sitios por grupo
    integer,parameter::NMG = 150 !(dimensi�n m�xima para vectores que guardan inf. sobre subgrupos)
    integer,parameter::NINT = 70 !Dimensi�n para vectores que guardan inf. sobre par�metros de interacci�n
    integer::ipareq,icalc,modelo,iprm,ioutm,novapm,igm
    integer::output
    character(len=36)::name
    real*8,dimension(NMG,NMG)::aint1
    integer:: NTEXT(36)
    real*8::ANT(10,3) 
    
endmodule InputData