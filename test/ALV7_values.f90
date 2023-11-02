module ALV7_values
    implicit none
    private
    public :: ALV7_data
    integer, parameter :: QR_K = selected_real_kind (4)
  
 
    ! These are module variables that can be used by any program unit that uses this module
    real(kind=QR_K), dimension(10) :: ALV7_data =&
    &(/0.44344967596404450,&
    &2.2377817584918511E-002,&
    &2.6720123199777176E-002,&
    &6.2876006809734569E-002,&
    &1.8026082405569736E-002,&
    &0.20655033403595552,&
    &1.0000000000000000E-8,&
    &0.33274674363719842,&
    &1.3739800562019173,&
    &0.24979099579751418/)
  end module ALV7_values