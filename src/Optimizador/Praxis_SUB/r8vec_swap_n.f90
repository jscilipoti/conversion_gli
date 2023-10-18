subroutine r8vec_swap_n ( n, a1, a2 )

    !*****************************************************************************80
    !
    !! r8vec_swap_n swaps the entries of two R8VECs.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries in the arrays.
    !
    !    Input/output, real ( kind = 8 ) A1(N), A2(N), the vectors to swap.
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) a1(n)
      real ( kind = 8 ) a2(n)
      real ( kind = 8 ) a3(n)
    
      a3(1:n) = a1(1:n)
      a1(1:n) = a2(1:n)
      a2(1:n) = a3(1:n)
    
      return
    end
    subroutine svsort_n ( n, d, v ) 
    
    !*****************************************************************************80
    !
    !! svsort_n descending sorts singular values D and adjusts V.
    !
    !  Discussion:
    !
    !    A simple bubble sort is used on D.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 February 2002
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
    !    Input, integer ( kind = 4 ) N, the length of D, and the order of V.
    !
    !    Input/output, real ( kind = 8 ) D(N), the vector to be sorted.  
    !    On output, the entries of D are in descending order.
    !
    !    Input/output, real ( kind = 8 ) V(N,N), an N by N array to be adjusted 
    !    as D is sorted.  In particular, if the value that was in D(I) on input is
    !    moved to D(J) on output, then the input column V(*,I) is moved to
    !    the output column V(*,J).
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      real ( kind = 8 ) d(n)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) j
      integer ( kind = 4 ) j2
      integer ( kind = 4 ) j3
      real ( kind = 8 ) t
      real ( kind = 8 ) v(n,n)
    
      do j = 1, n - 1
    
        j3 = j;
        do j2 = j + 1, n
          if ( d(j3) < d(j2) ) then
            j3 = j2
          end if
        end do
    
        t     = d(j)
        d(j)  = d(j3)
        d(j3) = t
    
        do i = 1, n
          t       = v(i,j)
          v(i,j)  = v(i,j3)
          v(i,j3) = t
        end do
    
      end do
    
      return
    end