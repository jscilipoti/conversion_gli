subroutine quad_n ( n, f_n, x, t, h, v, q0, q1, nl, nf, dmin, ldt, fx, qf1, &
    qa, qb, qc, qd0, qd1 )
  
  !*****************************************************************************80
  !
  !! quad_n seeks to minimize the scalar function f_n along a particular curve.
  !
  !  Discussion:
  !
  !    The minimizer to be sought is required to lie on a curve defined
  !    by Q0, Q1 and X.
  !
  !    This function was modified by removing the common blocks,
  !    28 July 2016.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    30 July 2016
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
  !    Input, external real ( kind = 8 ) f_n, is the name of the function to 
  !    be minimized.  The function should have the form 
  !      function f_n(x,n)
  !      integer ( kind = 4 ) n
  !      real ( kind = 8 ) f_n
  !      real ( kind = 8 ) x(n)
  !    and accepts X and N as input, returning in f_n the function value.
  !
  !    Input/output, real ( kind = 8 ) X(N), ?
  !
  !    Input, real ( kind = 8 ) T, ?
  !
  !    Input, real ( kind = 8 ) H, ?
  !
  !    Input, real ( kind = 8 ) V(N,N), the matrix of search directions.
  !
  !    Input/output, real ( kind = 8 ) Q0(N), Q1(N), auxiliary points used to \
  !    define a curve through X.
  !
  !    Input/output, integer ( kind = 4 ) NL, the number of linear searches.
  !
  !    Input/output, integer ( kind = 4 ) NF, the number of function evaluations.
  !
  !    Input, real ( kind = 8 ) DMIN, an estimate for the smallest eigenvalue.
  !
  !    Input, real ( kind = 8 ) LDT, the length of the step.
  !
  !    Input/output, real ( kind = 8 ) FX, the value of f_n(X,N).
  !
  !    Input/output, real ( kind = 8 ) QF1, QA, QB, QC, QD0, QD1 ?
  !
    implicit none
  
    integer ( kind = 4 ) n
  
    real ( kind = 8 ) dmin
    real ( kind = 8 ), external :: f_n
    logical fk
    real ( kind = 8 ) fx
    real ( kind = 8 ) h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) jsearch
    real ( kind = 8 ) l
    real ( kind = 8 ) ldt
    integer ( kind = 4 ) nf
    integer ( kind = 4 ) nits
    integer ( kind = 4 ) nl
    real ( kind = 8 ) q0(n)
    real ( kind = 8 ) q1(n)
    real ( kind = 8 ) qa
    real ( kind = 8 ) qb
    real ( kind = 8 ) qc
    real ( kind = 8 ) qd0
    real ( kind = 8 ) qd1
    real ( kind = 8 ) qf1
    real ( kind = 8 ) s
    real ( kind = 8 ) t
    real ( kind = 8 ) temp
    real ( kind = 8 ) v(n,n)
    real ( kind = 8 ) value
    real ( kind = 8 ) x(n)
  
    temp = fx
    fx   = qf1
    qf1  = temp
  
    call r8vec_swap_n ( n, x, q1 )
  
    qd1 = sqrt ( sum ( ( x(1:n) - q1(1:n) ) ** 2 ) )
  
    l = qd1
    s = 0.0D+00
  
    if ( qd0 <= 0.0D+00 .or. qd1 <= 0.0D+00 .or. nl < 3 * n * n ) then
  
      fx = qf1
      qa = 0.0D+00
      qb = 0.0D+00
      qc = 1.0D+00
  
    else
  
      jsearch = -1
      nits = 2
      value = qf1
      fk = .true.
  
      call minny_n ( n, jsearch, nits, s, l, value, fk, f_n, x, t, &
        h, v, q0, q1, nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 )
  
      qa =                 l * ( l - qd1 )       / ( qd0 + qd1 ) / qd0
      qb = - ( l + qd0 )     * ( l - qd1 ) / qd1                 / qd0
      qc =   ( l + qd0 ) * l               / qd1 / ( qd0 + qd1 )
  
    end if
  
    qd0 = qd1
  
    do i = 1, n
      s = q0(i)
      q0(i) = x(i)
      x(i) = qa * s + qb * x(i) + qc * q1(i)
    end do
  
    return
  end