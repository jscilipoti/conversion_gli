function r8_hypot ( x, y )

    !*****************************************************************************80
    !
    !! R8_HYPOT returns the value of sqrt ( X^2 + Y^2).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, Y, the arguments.
    !
    !    Output, real ( kind = 8 ) R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
    !
      implicit none
    
      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) r8_hypot
      real ( kind = 8 ) x
      real ( kind = 8 ) y
    
      if ( abs ( x ) < abs ( y ) ) then
        a = abs ( y )
        b = abs ( x )
      else
        a = abs ( x )
        b = abs ( y )
      end if
    !
    !  A contains the larger value.
    !
      if ( a == 0.0D+00 ) then
        c = 0.0D+00
      else
        c = a * sqrt ( 1.0D+00 + ( b / a ) ** 2 )
      end if
    
      r8_hypot = c
    
      return
    end