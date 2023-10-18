function flin ( n, jsearch, l, f, x, nf, v, q0, q1, qd0, qd1, qa, qb, qc )

    !*****************************************************************************80
    !
    !! FLIN is the function of one variable to be minimized by MINNY.
    !
    !  Discussion:
    !
    !    F(X) is a scalar function of a vector argument X.
    !
    !    A minimizer of F(X) is sought along a line or parabola.
    !
    !    This function has been modified, by removing the occurrence of a 
    !    common block, so that it looks more like a "normal" function that does
    !    not rely so strongly on peculiarities of FORTRAN.
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
    !    If J is a legal column index, linear search in direction of V(*,JSEARCH).
    !    Otherwise, then the search is parabolic, based on X, Q0 and Q1.
    !
    !    Input, real ( kind = 8 ) L, is the parameter determining the particular
    !    point at which F is to be evaluated.  
    !    For a linear search, L is the step size.
    !    For a quadratic search, L is a parameter which specifies
    !    a point in the plane of X, Q0 and Q1.
    !
    !    Input, external F, is the name of the function to be minimized.
    !    The function should have the form 
    !      function f(x,n)
    !      integer ( kind = 4 ) n
    !      real ( kind = 8 ) f
    !      real ( kind = 8 ) x(n)
    !    and accepts X and N as input, returning in F the function value.
    !
    !    Input, real ( kind = 8 ) X(N), the base point of the search.
    !
    !    Input/output, integer ( kind = 4 ) NF, the function evaluation counter.
    !
    !    Input, real ( kind = 8 ) V(N,N), a matrix whose columns constitute 
    !    search directions.
    !
    !    Input, real ( kind = 8 ) Q0(N), Q1(N), two auxiliary points used to
    !    determine the plane when a quadratic search is performed.
    !
    !    Input, real ( kind = 8 ) QD0, QD1, values needed to compute the 
    !    coefficients QA, QB, QC.
    !
    !    Output, real ( kind = 8 ) QA, QB, QC, coefficients used to combine
    !    Q0, X, and A1 if a quadratic search is used.
    !
    !    Output, real ( kind = 8 ) FLIN, the value of the function at the 
    !    minimizing point.
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ), external :: f
      real ( kind = 8 ) flin
      integer ( kind = 4 ) jsearch
      real ( kind = 8 ) l
      integer ( kind = 4 ) nf
      real ( kind = 8 ) q0(n)
      real ( kind = 8 ) q1(n)
      real ( kind = 8 ) qa
      real ( kind = 8 ) qb
      real ( kind = 8 ) qc
      real ( kind = 8 ) qd0
      real ( kind = 8 ) qd1
      real ( kind = 8 ) t(n)
      real ( kind = 8 ) v(n,n)
      real ( kind = 8 ) x(n)
    !
    !  The search is linear.
    !
      if ( 1 <= jsearch ) then
    
        t(1:n) = x(1:n) + l * v(1:n,jsearch)
    !
    !  The search is along a parabolic space curve.
    !
      else
    
        qa =                 l * ( l - qd1 ) /       ( qd0 + qd1 ) / qd0
        qb = - ( l + qd0 ) *     ( l - qd1 ) / qd1                 / qd0
        qc =   ( l + qd0 ) * l               / qd1 / ( qd0 + qd1 )
    
        t(1:n) = qa * q0(1:n) + qb * x(1:n) + qc * q1(1:n)
    
      end if
    !
    !  The function evaluation counter NF is incremented.
    !
      nf = nf + 1
    !
    !  Evaluate the function.
    !
      flin = f ( t, n )
    
      return
    end