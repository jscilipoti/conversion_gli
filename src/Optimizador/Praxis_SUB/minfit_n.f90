subroutine minfit_n ( n, tol, a, q )

    !*****************************************************************************80
    !
    !! minfit_n computes the singular value decomposition of an N by N array.
    !
    !  Discussion:
    !
    !    This is an improved version of the EISPACK routine minfit_n
    !    restricted to the case M = N and P = 0.
    !
    !    The singular values of the array A are returned in Q.  A is
    !    overwritten with the orthogonal matrix V such that U * diag(Q) = A * V,
    !    where U is another orthogonal matrix.
    !
    !    Thanks to Andreas Zuend for pointing out a potential for overflow in
    !    the computation z = sqrt ( f_n*f_n + 1 ), 22 March 2012.
    !
    !    Thanks to Martin Horvat for correcting the initial assignment of S,
    !    and for moving the line "L=I" back within the loop, 25 June 2013.
    !
    !    Thanks to Martin Horvat for moving the assignment "X=MAX(X,Y)"
    !    into the loop, 29 October 2013.
    !
    !    Further modifications to this function were made, to remove the use
    !    of statement labels, so that the function looks more like a "normal"
    !    function that relies much less on the peculiarities of FORTRAN,
    !    27 July 2016.
    !
    !    Removed argument M.  A must be dimensions (N,N).
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
    !    James Wilkinson, Christian Reinsch,
    !    Handbook for Automatic Computation,
    !    Volume II, Linear Algebra, Part 2,
    !    Springer Verlag, 1971.
    !
    !    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, Yasuhiko Ikebe, 
    !    Virginia Klema, Cleve Moler,
    !    Matrix Eigensystem Routines, EISPACK Guide,
    !    Lecture Notes in Computer Science, Volume 6,
    !    Springer Verlag, 1976,
    !    ISBN13: 978-3540075462,
    !    LC: QA193.M37.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the matrix A.
    !
    !    Input, real ( kind = 8 ) TOL, a tolerance which determines when a vector
    !    (a column or part of a column of the matrix) may be considered
    !    "essentially" equal to zero.
    !
    !    Input/output, real ( kind = 8 ) A(N,N).  On input, an N by N array whose
    !    singular value decomposition is desired.  On output, the
    !    SVD orthogonal matrix factor V.
    !
    !    Input/output, real ( kind = 8 ) Q(N), the singular values.
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) a(n,n)
      real ( kind = 8 ) c
      real ( kind = 8 ) e(n)
      real ( kind = 8 ) eps
      real ( kind = 8 ) f_n
      real ( kind = 8 ) g
      real ( kind = 8 ) h
      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      integer ( kind = 4 ) k
      integer ( kind = 4 ) kt
      integer ( kind = 4 ), parameter :: kt_max = 30
      integer ( kind = 4 ) l
      integer ( kind = 4 ) l2
      real ( kind = 8 ) q(n)
      real ( kind = 8 ) r8_hypot_n
      real ( kind = 8 ) s
      logical skip
      real ( kind = 8 ) temp
      real ( kind = 8 ) tol
      real ( kind = 8 ) x
      real ( kind = 8 ) y
      real ( kind = 8 ) z
    !
    !  Householder's reduction to bidiagonal form.
    !
      if ( n == 1 ) then
        q(1) = a(1,1)
        a(1,1) = 1.0D+00
        return
      end if
    
      eps = epsilon ( eps )
      g = 0.0D+00
      x = 0.0D+00
    
      do i = 1, n
    
        e(i) = g
        l = i + 1
    
        s = sum ( a(i:n,i) ** 2 )
    
        g = 0.0D+00
    
        if ( tol <= s ) then
    
          f_n = a(i,i)
    
          g = sqrt ( s )
          if ( 0.0D+00 <= f_n ) then
            g = - g
          end if
    
          h = f_n * g - s
          a(i,i) = f_n - g
    
          do j = l, n
    
            f_n = dot_product ( a(i:n,i), a(i:n,j) ) / h
    
            a(i:n,j) = a(i:n,j) + f_n * a(i:n,i)
    
          end do 
    
        end if
    
        q(i) = g
    
        s = sum ( a(i,l:n) ** 2 )
    
        g = 0.0D+00
    
        if ( tol <= s ) then
    
          if ( i /= n ) then
            f_n = a(i,i+1)
          end if
    
          g = sqrt ( s )
          if ( 0.0D+00 <= f_n ) then
            g = - g
          end if
    
          h = f_n * g - s
    
          if ( i /= n ) then
    
            a(i,i+1) = f_n - g
            e(l:n) = a(i,l:n) / h
    
            do j = l, n
    
              s = dot_product ( a(j,l:n), a(i,l:n) )
    
              a(j,l:n) = a(j,l:n) + s * e(l:n)
    
            end do
    
          end if
    
        end if
    
        y = abs ( q(i) ) + abs ( e(i) )
    
        x = max ( x, y )
    
      end do
    !
    !  Accumulation of right-hand transformations.
    !
      a(n,n) = 1.0D+00
      g = e(n)
      l = n
    
      do i = n - 1, 1, -1
    
        if ( g /= 0.0D+00 ) then
    
          h = a(i,i+1) * g
    
          a(l:n,i) = a(i,l:n) / h
    
          do j = l, n
    
            s = dot_product ( a(i,l:n), a(l:n,j) )
    
            a(l:n,j) = a(l:n,j) + s * a(l:n,i)
    
          end do
    
        end if
    
        a(i,l:n) = 0.0D+00
        a(l:n,i) = 0.0D+00
        a(i,i) = 1.0D+00
    
        g = e(i)
    
        l = i
    
      end do
    !
    !  Diagonalization of the bidiagonal form.
    !
      eps = eps * x
    
      do k = n, 1, -1
    
        kt = 0
    
        do
    
          kt = kt + 1
    
          if ( kt_max < kt ) then
            e(k) = 0.0D+00
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'minfit_n - Fatal error!'
            write ( *, '(a)' ) '  The QR algorithm failed to converge.'
            stop 1
          end if
    
          skip = .false.
    
          do l2 = k, 1, -1
    
            l = l2
    
            if ( abs ( e(l) ) <= eps ) then
              skip = .true.
              exit
            end if
    
            if ( l /= 1 ) then
              if ( abs ( q(l-1) ) <= eps ) then
                exit
              end if
            end if
    
          end do
    !
    !  Cancellation of E(L) if 1 < L.
    !
          if ( .not. skip ) then
    
            c = 0.0D+00
            s = 1.0D+00
    
            do i = l, k
       
              f_n = s * e(i)
              e(i) = c * e(i)
              if ( abs ( f_n ) <= eps ) then
                exit
              end if
              g = q(i)
    !
    !  q(i) = h = sqrt(g*g + f_n*f_n).
    !
              h = r8_hypot_n ( f_n, g )
      
              q(i) = h
    
              if ( h == 0.0D+00 ) then
                g = 1.0D+00
                h = 1.0D+00
              end if
    
              c =   g / h
              s = - f_n / h
    
            end do
    
          end if
    !
    !  Test for convergence for this index K.
    !
          z = q(k)
    
          if ( l == k ) then
            if ( z < 0.0D+00 ) then
              q(k) = - z
              a(1:n,k) = - a(1:n,k)
            end if
            exit
          end if
    !
    !  Shift from bottom 2*2 minor.
    !
          x = q(l)
          y = q(k-1)
          g = e(k-1)
          h = e(k)
          f_n = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) &
            / ( 2.0D+00 * h * y )
    
          g = r8_hypot_n ( f_n, 1.0D+00 )
    
          if ( f_n < 0.0D+00 ) then
            temp = f_n - g
          else
            temp = f_n + g
          end if
    
          f_n = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
    !
    !  Next QR transformation.
    !
          c = 1.0D+00
          s = 1.0D+00
    
          do i = l + 1, k
    
            g = e(i)
            y = q(i)
            h = s * g
            g = g * c
    
            z = r8_hypot_n ( f_n, h )
    
            e(i-1) = z
    
            if ( z == 0.0D+00 ) then
              f_n = 1.0D+00
              z = 1.0D+00
            end if
    
            c = f_n / z
            s = h / z
            f_n =   x * c + g * s
            g = - x * s + g * c
            h = y * s
            y = y * c
    
            do j = 1, n
              x = a(j,i-1)
              z = a(j,i)
              a(j,i-1) = x * c + z * s
              a(j,i) = - x * s + z * c
            end do
    
            z = r8_hypot_n ( f_n, h )
    
            q(i-1) = z
    
            if ( z == 0.0D+00 ) then
              f_n = 1.0D+00
              z = 1.0D+00
            end if
    
            c = f_n / z
            s = h / z
            f_n =   c * g + s * y
            x = - s * g + c * y
    
          end do
    
          e(l) = 0.0D+00
          e(k) = f_n
          q(k) = x
    
        end do
    
      end do
    
      return
    end