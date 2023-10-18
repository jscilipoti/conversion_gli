subroutine print2_n ( n, x, prin, fx, nf, nl )

    !*****************************************************************************80
    !
    !! print2_n prints certain data about the progress of the iteration.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 May 2006
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Richard Brent.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Richard Brent,
    !    Algorithms for Minimization with Derivatives,
    !    Prentice Hall, 1973,
    !    Reprinted by Dover, 2002.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of variables.
    !
    !    Input, real ( kind = 8 ) X(N), the current estimate of the minimizer.
    !
    !    Input, integer ( kind = 4 ) PRIN, the user-specifed print level.
    !    0, nothing is printed.
    !    1, f_n is printed after every n+1 or n+2 linear minimizations.  
    !       final X is printed, but intermediate X is printed only 
    !       if N is at most 4.
    !    2, the scale factors and the principal values of the approximating 
    !       quadratic form are also printed.
    !    3, X is also printed after every few linear minimizations.
    !    4, the principal vectors of the approximating quadratic form are 
    !       also printed.
    !
    !    Input, real ( kind = 8 ) FX, the smallest value of f_n(X) found so far.
    !
    !    Input, integer ( kind = 4 ) NF, the number of function evaluations.
    !
    !    Input, integer ( kind = 4 ) NL, the number of linear searches.
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) fx
      integer ( kind = 4 ) nf
      integer ( kind = 4 ) nl
      integer ( kind = 4 ) prin
      real ( kind = 8 ) x(n)
    
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Linear searches      ', nl
      write ( *, '(a,i8)' ) '  Function evaluations ', nf 
    
      if ( n <= 4 .or. 2 < prin ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'X:'
        write ( *, '(5g14.6)' ) x(1:n)
      end if
    
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The function value FX: ', fx
    
      return
    end