subroutine svsort ( n, d, v ) 

    !*****************************************************************************80
    !
    !! SVSORT descending sorts singular values D and adjusts V.
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