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
 
      subroutine ixsrt1( nx, ix, xx )

*---  IDSRT1: indirect sort -- sort ix and carry xx along
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
*----------------------------------------------------------------------|

      implicit none
      integer          nx
      integer, dimension(nx) :: ix
c      REAL(kind=selected_real_kind(15)), dimension(nx) :: xx

 
 300  CONTINUE
      RETURN
      END
