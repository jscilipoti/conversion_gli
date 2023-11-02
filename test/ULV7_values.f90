module ULV7_values
    implicit none
    private
    public :: ULV7_data
    integer, parameter :: QR_K = selected_real_kind (4)
  
 
    ! These are module variables that can be used by any program unit that uses this module
    real(kind=QR_K), dimension(10) :: ULV7_data =&
    (/0.47108774710957019,&
    &2.9794662252076217E-2,&
    &3.1367147395115078E-2,&
    &5.8969495563111381E-2,&
    &9.8687247896973255E-3,&
    &0.17891226289042983,&
    &.0000000000000000E-8,&
    &0.33274694860814902,&
    &1.3739778518466321,&
    &0.24979108614114404/)
  end module ULV7_values