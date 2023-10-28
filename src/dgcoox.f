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
      subroutine ixsrt1( xx )

*---  IDSRT1: indirect sort -- sort ix and carry xx along
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
*----------------------------------------------------------------------|
c      implicit none
      integer          nx

      complex(kind=8) xx(nx)

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
      RETURN
      END

