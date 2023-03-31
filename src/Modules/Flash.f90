module flash
    
    real*8::z(10),P,T,xguess=0.6,aglim
    integer::contador=1
    logical::agextr
    
endmodule flash
    
module flashout
    real*8::agam(10,4),compfases(10,4)
    integer::cont
endmodule flashout