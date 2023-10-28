* Copyright: See /inst/LAPACK_LICENSE.txt for 
* original FORTRAN code in /src.
*
* The FORTRAN lapack/blas code in rexpokit was 
* originally copied from the EXPOKIT package
* with permission of Roger Sidje (who is
* thus listed as coauthor on rexpokit).
*
* The FORTRAN has since had various minor 
* modifications to satisfy new checks as
* CRAN updates their FORTRAN, OSs, and
* R CMD check function.
* 

* 2023-10-28:
* Fix: 
* Version: 0.26.6.9
* Check: usage of KIND in Fortran files
* Result: WARN
*     Found the following files with non-portable usage of KIND:
*      itscale5.f
*      mataid.f
*      my_expokit.f
* 
* mataid.f
* c      complex(kind=8)
*       complex
*
*
* double precision
* replaced with:
* REAL(kind=selected_real_kind(15)) ::
*
* complex
* replaced with
* complex(kind=selected_real_kind(15)) :: 



*-------------------------------NOTE-----------------------------------*
*     This is an accessory to Expokit and it is not intended to be     *
*     complete. It is supplied primarily to ensure an unconstrained    *
*     distribution and portability of the package. The matrix-vector   *
*     multiplication routines supplied here fit the non symmetric      *
*     storage and for a symmetric matrix, the entire (not half) matrix *
*     is required.  If the sparsity pattern is known a priori, it is   *
*     recommended to use the most advantageous format and to devise    *
*     the most advantageous matrix-vector multiplication routine.      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
 
c      subroutine ixsrt1( nx, ix, xx )
c      subroutine ixsrt1( xx )

*---  IDSRT1: indirect sort -- sort ix and carry xx along
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
*----------------------------------------------------------------------|
c      implicit none
c      integer          nx

c      complex(kind=8) xx(1)

c      integer, dimension(nx) :: ix

c ERROR:
c     REAL(kind=selected_real_kind(15)), dimension(nx) :: xx
c			complex xx
c			complex, dimension(nx) :: xx
c			double complex :: xx(nx)
c			REAL(KIND=selected_real_kind(p=15)), dimension(nx) :: xx

c Compile failure:
c			use iso_fortran_env, only: xx(nx) => real64
c			REAL xx(nx)
c			complex, parameter :: xx(nx) = selected_real_kind(33, 4931)
c			REAL xx(nx)

c			COMPLEX(KIND=selected_complex_kind(15,15)) :: xx(nx)

c     SAME ERROR
c			COMPLEX*16 :: xx(nx)
      
c      USE ISO_FORTRAN_ENV, xx => real64
c      COMPLEX(KIND=8) xx(nx)
			
c 300  CONTINUE
c      RETURN
c      END


      subroutine  zcopy(n,zx,incx,zy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 4/11/78.
c
      complex(kind=8) zx(1),zy(1)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return


