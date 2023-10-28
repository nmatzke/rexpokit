      subroutine dgcoov ( x, y )
      implicit none

c     2023-10-28:
c      double precision x(*), y(*)
      REAL, dimension(:), allocatable :: x, y
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed here to be under the COOrdinates storage format.
*
      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer, dimension(nzmax) :: ia
      integer, dimension(nzmax) :: ja
      REAL, dimension(nzmax) :: a
      common /RMAT/ a, ia, ja, nz, n
      integer i, j
 
      do j = 1,n
         y(j) = 0.0d0
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END
