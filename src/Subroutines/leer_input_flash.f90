subroutine leer_input_flash()
    use InputData
    use flash

    implicit none

    integer :: N, i, j, k, ng
    real(8) :: Tx, px
    common /CUFAC/ N, NG, Px(10,10), Tx

    call open_file_name()
    open(2, file=name, status='old', form='formatted')
    read(2, '(36A2)') NTEXT
    read(2, *) ICALC, modelo, IPRm, IOUTm, NOVAPm, igm, ipareq
    call ab_ban1(modelo)
    call PARIN2
    if (NOVAPm /= 0) then
        do j = 1, N
            read(2, *) (ANT(k,j), k = 1, 3)
        end do
        do j = 1, N
            ANT(1,j) = 2.302585 * (ANT(1,j) - 2.880814)
            ANT(2,j) = 2.302585 * ANT(2,j)
        end do
    endif
    read(2,'(2F10.2)') T,P
    read(2, *) (Z(i), i = 1, N)
end subroutine leer_input_flash
