subroutine minny ( n, jsearch, nits, d2, x1, f1, fk, f, x, t, h, v, q0, q1, &
    nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 )
  
  !*****************************************************************************80
  !
  !! MINNY minimizes a scalar function of N variables along a line.
  !
  !  Discussion:
  !
  !    MINNY minimizes F along the line from X in the direction V(*,J) unless
  !    J is less than 1, when a quadratic search is made in the plane
  !    defined by Q0, Q1 and X.
  !
  !    If FK = true, then F1 is FLIN(X1).  Otherwise X1 and F1 are ignored
  !    on entry unless final FX is greater than F1.
  !
  !    This function was modified by removing the common blocks
  !    and the use of labeled statements, 28 July 2016.
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
  !    Input, integer ( kind = 4 ) JSEARCH, indicates the kind of search.
  !    If JSEARCH is a legal column index, linear search in direction of V(*,J).
  !    Otherwise, the search is parabolic, based on X, Q0 and Q1.
  !
  !    Input, integer ( kind = 4 ) NITS, the maximum number of times the interval 
  !    may be halved to retry the calculation.
  !
  !    Input/output, real ( kind = 8 ) D2, is either zero, or an approximation to 
  !    the value of (1/2) times the second derivative of F.
  !
  !    Input/output, real ( kind = 8 ) X1, on entry, an estimate of the 
  !    distance from X to the minimum along V(*,J), or, if J = 0, a curve.  
  !    On output, the distance between X and the minimizer that was found.
  !
  !    Input/output, real ( kind = 8 ) F1, ?
  !
  !    Input, logical FK; if FK is TRUE, then on input F1 contains 
  !    the value FLIN(X1).
  !
  !    Input, external real ( kind = 8 ) F, is the name of the function to 
  !    be minimized.  The function should have the form 
  !      function f(x,n)
  !      integer ( kind = 4 ) n
  !      real ( kind = 8 ) f
  !      real ( kind = 8 ) x(n)
  !    and accepts X and N as input, returning in F the function value.
  !
  !    Input/output, real ( kind = 8 ) X(N), ?
  !
  !    Input, real ( kind = 8 ) T, ?
  !
  !    Input, real ( kind = 8 ) H, ?
  !
  !    Input, real ( kind = 8 ) V(N,N), a matrix whose columns are direction
  !    vectors along which the function may be minimized.
  !
  !    Input, real ( kind = 8 ) Q0(N), an auxiliary point used to define
  !    a curve through X.
  !
  !    Input, real ( kind = 8 ) Q1(N), an auxiliary point used to define
  !    a curve through X.
  !
  !    Input/output, integer ( kind = 4 ) NL, the number of linear searches.
  !
  !    Input/output, integer ( kind = 4 ) NF, the number of function evaluations.
  !
  !    Input, real ( kind = 8 ) DMIN, an estimate for the smallest eigenvalue.
  !
  !    Input, real ( kind = 8 ) LDT, the length of the step.
  !
  !    Input/output, real ( kind = 8 ) FX, the value of F(X,N).
  !
  !    Input/output, real ( kind = 8 ) QA, QB, QC, ?
  !
  !    Input, real ( kind = 8 ) QD0, QD1, ?.
  !
    implicit none
  
    real ( kind = 8 ) d1
    real ( kind = 8 ) d2
    real ( kind = 8 ) dmin
    logical dz
    real ( kind = 8 ), external :: f
    real ( kind = 8 ) f0
    real ( kind = 8 ) f1
    real ( kind = 8 ) f2
    logical fk
    real ( kind = 8 ) flin
    real ( kind = 8 ) fm
    real ( kind = 8 ) fx
    real ( kind = 8 ) h
    integer ( kind = 4 ) jsearch
    integer ( kind = 4 ) k
    real ( kind = 8 ) ldt
    real ( kind = 8 ) m2
    real ( kind = 8 ) m4
    real ( kind = 8 ) machep
    integer ( kind = 4 ) n
    integer ( kind = 4 ) nf
    integer ( kind = 4 ) nits
    integer ( kind = 4 ) nl
    logical ok
    real ( kind = 8 ) q0(n)
    real ( kind = 8 ) q1(n)
    real ( kind = 8 ) qa
    real ( kind = 8 ) qb
    real ( kind = 8 ) qc
    real ( kind = 8 ) qd0
    real ( kind = 8 ) qd1
    real ( kind = 8 ) s
    real ( kind = 8 ) sf1
    real ( kind = 8 ) small
    real ( kind = 8 ) sx1
    real ( kind = 8 ) t
    real ( kind = 8 ) t2
    real ( kind = 8 ) temp
    real ( kind = 8 ) v(n,n)
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) x1
    real ( kind = 8 ) x2
    real ( kind = 8 ) xm
  
    machep = epsilon ( machep )
    small = machep ** 2
    m2 = sqrt ( machep )
    m4 = sqrt ( m2 )
    sf1 = f1
    sx1 = x1
    k = 0
    xm = 0.0D+00
    fm = fx
    f0 = fx
    dz = ( d2 < machep )
  !
  !  Find the step size.
  ! 
    s = sqrt ( sum ( x(1:n) ** 2 ) )
  
    if ( dz ) then
      temp = dmin
    else
      temp = d2
    end if
  
    t2 = m4 * sqrt ( abs ( fx ) / temp + s * ldt ) + m2 * ldt
    s = m4 * s + t
    if ( dz .and. s < t2 ) then
      t2 = s
    end if
  
    t2 = max ( t2, small )
    t2 = min ( t2, 0.01D+00 * h )
  
    if ( fk .and. f1 <= fm ) then
      xm = x1
      fm = f1
    end if
  
    if ( .not. fk .or. abs ( x1 ) < t2 ) then
  
      if ( 0.0D+00 <= x1 ) then
        temp = 1.0D+00
      else
        temp = - 1.0D+00
      end if
  
      x1 = temp * t2
      f1 = flin ( n, jsearch, x1, f, x, nf, v, q0, q1, qd0, qd1, qa, qb, qc )
  
    end if
  
    if ( f1 <= fm ) then
      xm = x1
      fm = f1
    end if
  !
  !  Evaluate FLIN at another point and estimate the second derivative.
  !
    do
  
      if ( dz ) then
  
        if ( f1 <= f0 ) then
          x2 = 2.0D+00 * x1
        else
          x2 = - x1
        end if
  
        f2 = flin ( n, jsearch, x2, f, x, nf, v, q0, q1, qd0, qd1, qa, qb, qc )
  
        if ( f2 <= fm ) then
          xm = x2
          fm = f2
        end if
  
        d2 = ( x2 * ( f1 - f0 ) - x1 * ( f2 - f0 ) ) &
          / ( ( x1 * x2 ) * ( x1 - x2 ) )
  
      end if
  !
  !  Estimate the first derivative at 0.
  !
      d1 = ( f1 - f0 ) / x1 - x1 * d2
      dz = .true.
  !
  !  Predict the minimum.
  !
      if ( d2 <= small ) then
  
        if ( 0.0D+00 <= d1 ) then
          x2 = - h
        else
          x2 = h
        end if
  
      else
  
        x2 = ( - 0.5D+00 * d1 ) / d2
  
      end if
  
      if ( h < abs ( x2 ) ) then
  
        if ( x2 <= 0.0D+00 ) then
          x2 = - h
        else
          x2 = h
        end if
  
      end if
  !
  !  Evaluate F at the predicted minimum.
  !
      ok = .true.
  
      do
  
        f2 = flin ( n, jsearch, x2, f, x, nf, v, q0, q1, qd0, qd1, qa, qb, qc )
  
        if ( nits <= k .or. f2 <= f0 ) then
          exit
        end if
  
        k = k + 1
  
        if ( f0 < f1 .and. 0.0D+00 < x1 * x2 ) then
          ok = .false.
          exit
        end if
  
        x2 = 0.5D+00 * x2
  
      end do
  
      if ( ok ) then
        exit
      end if
  
    end do
  !
  !  Increment the one-dimensional search counter.
  !
    nl = nl + 1
  
    if ( fm < f2 ) then
      x2 = xm
    else
      fm = f2
    end if
  !
  !  Get a new estimate of the second derivative.
  !
    if ( small < abs ( x2 * ( x2 - x1 ) ) ) then
      d2 = ( x2 * ( f1 - f0 ) - x1 * ( fm - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )
    else
      if ( 0 < k ) then
        d2 = 0.0D+00
      end if
    end if
  
    d2 = max ( d2, small )
  
    x1 = x2
    fx = fm
  
    if ( sf1 < fx ) then
      fx = sf1
      x1 = sx1
    end if
  !
  !  Update X for linear but not parabolic search.
  !
    if ( 1 <= jsearch ) then
  
      x(1:n) = x(1:n) + x1 * v(1:n,jsearch)
  
    end if
  
    return
  end