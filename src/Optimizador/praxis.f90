function praxis ( t0, h0, n, prin, x, f )

!*****************************************************************************80
!
!! PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
!
!  Discussion:
!
!    PRAXIS returns the minimum of the function F(X,N) of N variables
!    using the principal axis method.  The gradient of the function is
!    not required.
!
!    The approximating quadratic form is
!
!      Q(x') = F(x,n) + (1/2) * (x'-x)' * A * (x'-x)
!
!    where X is the best estimate of the minimum and 
!
!      A = inverse(V') * D * inverse(V)
!
!    V(*,*) is the matrix of search directions; 
!    D(*) is the array of second differences.  
!
!    If F(X) has continuous second derivatives near X0, then A will tend 
!    to the hessian of F at X0 as X approaches X0.
!
!    Thanks to Andreas Zuend for pointing out an error in the form of the
!    call to the routine r8mat_print (), 22 March 2012.
!
!    This function was modified by eliminating the use of labeled statements,
!    and removing the common blocks, 28 July 2016.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 July 2016
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
!    Input, real ( kind = 8 ) T0, is a tolerance.  PRAXIS attempts to return 
!    praxis = f(x) such that if X0 is the true local minimum near X, then
!    norm ( x - x0 ) < T0 + sqrt ( EPSILON ( X ) ) * norm ( X ),
!    where EPSILON ( X ) is the machine precision for X.
!
!    Input, real ( kind = 8 ) H0, is the maximum step size.  H0 should be 
!    set to about the maximum distance from the initial guess to the minimum.
!    If H0 is set too large or too small, the initial rate of
!    convergence may be slow.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!
!    Input, integer ( kind = 4 ) PRIN, controls printing intermediate results.
!    0, nothing is printed.
!    1, F is printed after every n+1 or n+2 linear minimizations.  
!       final X is printed, but intermediate X is printed only 
!       if N is at most 4.
!    2, the scale factors and the principal values of the approximating 
!       quadratic form are also printed.
!    3, X is also printed after every few linear minimizations.
!    4, the principal vectors of the approximating quadratic form are 
!       also printed.
!
!    Input/output, real ( kind = 8 ) X(N), is an array containing on entry a
!    guess of the point of minimum, on return the estimated point of minimum.
!
!    Input, external real ( kind = 8 ) F, is the name of the function to be
!    minimized.  The function should have the form 
!      function f(x,n)
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!    and accepts X and N as input, returning in F the function value.
!
!    Output, real ( kind = 8 ) PRAXIS, the function value at the minimizer.
!
!  Local parameters:
!
!    Local, real ( kind = 8 ) DMIN, an estimate for the smallest eigenvalue.
!
!    Local, real ( kind = 8 ) FX, the value of F(X,N).
!
!    Local, logical ILLC, is TRUE if the system is ill-conditioned.
!
!    Local, real ( kind = 8 ) LDT, the length of the step.
!
!    Local, integer ( kind = 4 ) NF, the number of function evaluations.
!
!    Local, integer ( kind = 4 ) NL, the number of linear searches.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) d(n)
  real ( kind = 8 ) d2
  real ( kind = 8 ) df
  real ( kind = 8 ) dmin
  real ( kind = 8 ) dn
  real ( kind = 8 ) dni
  real ( kind = 8 ), external :: f
  real ( kind = 8 ) f1
  logical fk
  real ( kind = 8 ) fx
  real ( kind = 8 ) h
  real ( kind = 8 ) h0
  integer ( kind = 4 ) i
  logical illc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jsearch
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k2
  integer ( kind = 4 ) kl
  integer ( kind = 4 ) kt
  integer ( kind = 4 ) ktm
  real ( kind = 8 ) large
  real ( kind = 8 ) ldfac
  real ( kind = 8 ) lds
  real ( kind = 8 ) ldt
  real ( kind = 8 ) m2
  real ( kind = 8 ) m4
  real ( kind = 8 ) machep
  integer ( kind = 4 ) nits
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nf
  real ( kind = 8 ) praxis
  integer ( kind = 4 ) prin
  real ( kind = 8 ) q0(n)
  real ( kind = 8 ) q1(n)
  real ( kind = 8 ) qa
  real ( kind = 8 ) qb
  real ( kind = 8 ) qc
  real ( kind = 8 ) qd0
  real ( kind = 8 ) qd1
  real ( kind = 8 ) qf1
  real ( kind = 8 ) r
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) s
  real ( kind = 8 ) scbd
  integer ( kind = 4 ) seed
  real ( kind = 8 ) sf
  real ( kind = 8 ) sl
  real ( kind = 8 ) small
  real ( kind = 8 ) t
  real ( kind = 8 ) t0
  real ( kind = 8 ) t2
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ) value
  real ( kind = 8 ) vlarge
  real ( kind = 8 ) vsmall
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) z(n)
  integer::loop
!
!  Initialization.
!
  machep = epsilon ( machep )
  small = machep * machep
  vsmall = small * small
  large = 1.0D+00 / small
  vlarge = 1.0D+00 / vsmall
  m2 = sqrt ( machep )
  m4 = sqrt ( m2 )
  seed = 123456789
!
!  Heuristic numbers:
!
!  If the axes may be badly scaled (which is to be avoided if
!  possible), then set SCBD = 10.  Otherwise set SCBD = 1.
!
!  If the problem is known to be ill-conditioned, initialize ILLC = true.
!
!  KTM is the number of iterations without improvement before the
!  algorithm terminates.  KTM = 4 is very cautious; usually KTM = 1
!  is satisfactory.
!
  scbd = 1.0D+00
  illc = .false.
  ktm = 1

  if ( illc ) then
    ldfac = 0.1D+00
  else
    ldfac = 0.01D+00
  end if

  kt = 0
  nl = 0
  nf = 1
  fx = f ( x, n )
  qf1 = fx
  t = small + abs ( t0 )
  t2 = t
  dmin = small
  h = h0
  h = max ( h, 100.0D+00 * t )
  ldt = h
!
!  The initial set of search directions V is the identity matrix.
!
  v(1:n,1:n) = 0.0D+00
  do i = 1, n
    v(i,i) = 1.0D+00
  end do

  d(1:n) = 0.0D+00
  qa = 0.0D+00
  qb = 0.0D+00
  qc = 0.0D+00
  qd0 = 0.0D+00
  qd1 = 0.0D+00
  q0(1:n) = x(1:n)
  q1(1:n) = x(1:n)

  if ( 0 < prin ) then
    call print2 ( n, x, prin, fx, nf, nl )
  end if
!
!  The main loop starts here.
!
  do

    sf = d(1)
    d(1) = 0.0D+00
!
!  Minimize along the first direction V(*,1).
!
    jsearch = 1
    nits = 2
    d2 = d(1)
    s = 0.0D+00
    value = fx
    fk = .false.

    call minny ( n, jsearch, nits, d2, s, value, fk, f, x, t, &
      h, v, q0, q1, nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 )

    d(1) = d2

    if ( s <= 0.0D+00 ) then
      v(1:n,1) = - v(1:n,1)
    end if

    if ( sf <= 0.9D+00 * d(1) .or. d(1) <= 0.9D+00 * sf ) then
      d(2:n) = 0.0D+00
    end if
!
!  The inner loop starts here.
!
    loop=2
    if(n==1)loop=1
    do k = 2, n

      y(1:n) = x(1:n)

      sf = fx

      if ( 0 < kt ) then
        illc = .true.
      end if

      do

        kl = k
        df = 0.0D+00
!
!  A random step follows, to avoid resolution valleys.
!
        if ( illc ) then

          do i = 1, n
            r = r8_uniform_01 ( seed )
!           call random_number ( harvest = r )
            s = ( 0.1D+00 * ldt + t2 * 10.0D+00 ** kt ) * ( r - 0.5D+00 )
            z(i) = s
            x(1:n) = x(1:n) + s * v(1:n,i)
          end do

          fx = f ( x, n )
          nf = nf + 1

        end if
!
!  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
!
        do k2 = k, n

          sl = fx

          jsearch = k2
          nits = 2
          d2 = d(k2)
          s = 0.0D+00
          value = fx
          fk = .false.

          call minny ( n, jsearch, nits, d2, s, value, fk, f, x, t, &
            h, v, q0, q1, nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 )

          d(k2) = d2

          if ( illc ) then
            s = d(k2) * ( ( s + z(k2) ) ** 2 )
          else
            s = sl - fx
          end if

          if ( df <= s ) then
            df = s
            kl = k2
          end if

        end do
!
!  If there was not much improvement on the first try, set
!  ILLC = true and start the inner loop again.
!
        if ( illc ) then
          exit
        end if

        if ( abs ( 100.0D+00 * machep * fx ) <= df ) then
          exit
        end if

        illc = .true.

      end do

      if ( k == 2 .and. 1 < prin ) then
        call r8vec_print ( n, d, '  The second difference array:' )
      end if
!
!  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
!
      do k2 = 1, k - 1

        jsearch = k2
        nits = 2
        d2 = d(k2)
        s = 0.0D+00
        value = fx
        fk = .false.

        call minny ( n, jsearch, nits, d2, s, value, fk, f, x, t, &
          h, v, q0, q1, nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 )

        d(k2) = d2

      end do

      f1 = fx
      fx = sf
      lds = 0.0D+00

      do i = 1, n
        sl = x(i)
        x(i) = y(i)
        sl = sl - y(i)
        y(i) = sl
        lds = lds + sl ** 2
      end do

      lds = sqrt ( lds )
!
!  Discard direction V(*,kl).
!
!  If no random step was taken, V(*,KL) is the "non-conjugate"
!  direction along which the greatest improvement was made.
!
      if ( small < lds ) then

        do j = kl - 1, k, -1
          v(1:n,j+1) = v(1:n,j)
          d(j+1) = d(j)
        end do

        d(k) = 0.0D+00

        v(1:n,k) = y(1:n) / lds
!
!  Minimize along the new "conjugate" direction V(*,k), which is
!  the normalized vector:  (new x) - (old x).
!
        jsearch = k
        nits = 4
        d2 = d(k)
        value = f1
        fk = .true.

        call minny ( n, jsearch, nits, d2, lds, value, fk, f, x, t, &
          h, v, q0, q1, nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 )

        d(k) = d2

        if ( lds <= 0.0D+00 ) then
          lds = - lds
          v(1:n,k) = - v(1:n,k)
        end if

      end if

      ldt = ldfac * ldt
      ldt = max ( ldt, lds )

      if ( 0 < prin ) then
        call print2 ( n, x, prin, fx, nf, nl )
      end if

      t2 = m2 * sqrt ( sum ( x(1:n) ** 2 ) ) + t
!
!  See whether the length of the step taken since starting the
!  inner loop exceeds half the tolerance.
!
      if ( 0.5D+00 * t2 < ldt ) then
        kt = - 1
      end if

      kt = kt + 1

      if ( ktm < kt ) then

        if ( 0 < prin ) then
          call r8vec_print ( n, x, '  X:' )
        end if

        praxis = fx

        return

      end if

    end do
!
!  The inner loop ends here.
!
!  Try quadratic extrapolation in case we are in a curved valley.
!
    call quad ( n, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, fx, qf1, &
      qa, qb, qc, qd0, qd1 )

    d(1:n) = 1.0D+00 / sqrt ( d(1:n) )

    dn = maxval ( d(1:n) )

    if ( 3 < prin ) then
      call r8mat_print ( n, n, v, '  The new direction vectors:' )
    end if

    do j = 1, n
      v(1:n,j) = ( d(j) / dn ) * v(1:n,j)
    end do
!
!  Scale the axes to try to reduce the condition number.
!
    if ( 1.0D+00 < scbd ) then

      do i = 1, n
        z(i) = max ( m4, sqrt ( sum ( v(i,1:n) ** 2 ) ) )
      end do

      s = minval ( z(1:n) )

      do i = 1, n

        sl = s / z(i)
        z(i) = 1.0D+00 / sl

        if ( scbd < z(i) ) then
          sl = 1.0D+00 / scbd
          z(i) = scbd
        end if

        v(i,1:n) = sl * v(i,1:n)

      end do

    end if
!
!  Calculate a new set of orthogonal directions before repeating
!  the main loop.
!
!  Transpose V for MINFIT:
!
    v(1:n,1:n) = transpose ( v(1:n,1:n) )
!
!  Call MINFIT to find the singular value decomposition of V.
!
!  This gives the principal values and principal directions of the
!  approximating quadratic form without squaring the condition number.
!
    call minfit ( n, vsmall, v, d )
!
!  Unscale the axes.
!
    if ( 1.0D+00 < scbd ) then

      do i = 1, n
        v(i,1:n) = z(i) * v(i,1:n)
      end do

      do i = 1, n

        s = sqrt ( sum ( v(1:n,i) ** 2 ) )
        d(i) = s * d(i)
        v(1:n,i) = v(1:n,i) / s

      end do

    end if

    do i = 1, n

      dni = dn * d(i)

      if ( large < dni ) then
        d(i) = vsmall
      else if ( dni < small ) then
        d(i) = vlarge
      else
        d(i) = 1.0D+00 / dni ** 2
      end if

    end do
!
!  Sort the singular values and singular vectors.
!
    call svsort ( n, d, v )
!
!  Determine the smallest eigenvalue.
!
    dmin = max ( d(n), small )
!
!  The ratio of the smallest to largest eigenvalue determines whether
!  the system is ill conditioned.
!
    if ( dmin < m2 * d(1) ) then
      illc = .true.
    else
      illc = .false.
    end if

    if ( 1 < prin ) then

      if ( 1.0D+00 < scbd ) then
        call r8vec_print ( n, z, '  The scale factors:' )
      end if 

      call r8vec_print ( n, d, '  Principal values of the quadratic form:' )

    end if

    if ( 3 < prin ) then
      call r8mat_print ( n, n, v, '  The principal axes:' )
    end if
!
!  The main loop ends here.
!
  end do

  if ( 0 < prin ) then
    call r8vec_print ( n, x, '  X:' )
  end if

  praxis = fx

  return
end