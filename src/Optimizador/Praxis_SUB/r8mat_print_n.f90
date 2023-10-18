subroutine r8mat_print_n ( m, n, a, title )

    !*****************************************************************************80
    !
    !! r8mat_print_n prints an R8MAT.
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
    !    12 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of rows in A.
    !
    !    Input, integer ( kind = 4 ) N, the number of columns in A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
      implicit none
    
      integer ( kind = 4 ) m
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) a(m,n)
      character ( len = * ) title
    
      call r8mat_print_some_n ( m, n, a, 1, 1, m, n, title )
    
      return
    end