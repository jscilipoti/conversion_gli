subroutine r8mat_print_some_n ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! r8mat_print_some_n prints some of an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
    !
    !    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
    !
    !    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
      implicit none
    
      integer ( kind = 4 ), parameter :: incx = 5
      integer ( kind = 4 ) m
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) a(m,n)
      character ( len = 14 ) ctemp(incx)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) i2hi
      integer ( kind = 4 ) i2lo
      integer ( kind = 4 ) ihi
      integer ( kind = 4 ) ilo
      integer ( kind = 4 ) inc
      integer ( kind = 4 ) j
      integer ( kind = 4 ) j2
      integer ( kind = 4 ) j2hi
      integer ( kind = 4 ) j2lo
      integer ( kind = 4 ) jhi
      integer ( kind = 4 ) jlo
      character ( len = * ) title
    
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
    
      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx
    
        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )
    
        inc = j2hi + 1 - j2lo
    
        write ( *, '(a)' ) ' '
    
        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i8,6x)' ) j
        end do
    
        write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '
    
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )
    
        do i = i2lo, i2hi
    
          do j2 = 1, inc
    
            j = j2lo - 1 + j2
    
            if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
            else
              write ( ctemp(j2), '(g14.6)' ) a(i,j)
            end if
    
          end do
    
          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )
    
        end do
    
      end do
    
      return
    end