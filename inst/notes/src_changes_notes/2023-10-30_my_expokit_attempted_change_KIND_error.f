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


* 2019-06-26 NJM edits:
* 
* fixing by changing 

*  DASUM  to DASUMX
*  DAXPY  to DAXPX (for reasons of lines too long)
*  DCOPY  to DCOPYX
*  DDOT   to DDOTX
*  DGEMM  to DGEXX
*  DGEMX  to DGEMX
*  DNRM2  to DNRM2X
*  DSCAL  to DSCALX
*  ZSWAP  to ZSWAPX
*  ZAXPY  to ZAXPX  (for reasons of lines too long)
*  
* Fixed these errors:
* 
* rexpokit.out:(.text+0x0): multiple definition of `lsame_'
* rexpokit.out:(.text+0x0): multiple definition of `dasum_'
* rexpokit.out:(.text+0x0): multiple definition of `DAXPX_'
* rexpokit.out:(.text+0x0): multiple definition of `DCOPYX_'
* rexpokit.out:(.text+0x0): multiple definition of `DDOTX_'
* rexpokit.out:(.text+0x0): multiple definition of `DGEXX_'
* rexpokit.out:(.text+0x0): multiple definition of `DGEMX_'
* rexpokit.out:(.text+0x0): multiple definition of `DNRM2X_'
* rexpokit.out:(.text+0x0): multiple definition of `DSCALX_'
* rexpokit.out:(.text+0x0): multiple definition of `dswap_'
* rexpokit.out:(.text+0x0): multiple definition of `idamax_'
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
* my_expokit.f




*----------------------------------------------------------------------|
* myDMEXPV:
      subroutine myDMEXPV( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, itrace,iflag,ia,ja,a,nz )

      implicit none
      integer n,nz,m,lwsp,liwsp, itrace,iflag,iwsp(liwsp),ia(nz),ja(nz)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp), a(nz)

*-----Purpose----------------------------------------------------------|
*
*---  DMEXPV computes w = exp(t*A)*v - Customised for MARKOV CHAINS.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*     This is a customised version for Markov Chains. This means that a
*     check is done within this code to ensure that the resulting vector 
*     w is a probability vector, i.e., w must have all its components 
*     in [0,1], with sum equal to 1. This check is done at some expense
*     and the user may try DGEXPV which is cheaper since it ignores 
*     probability constraints.
*
*     IMPORTANT: The check assumes that the transition rate matrix Q
*                satisfies Qe = 0, where e=(1,...,1)'. Don't use DMEXPV
*                if this condition does not hold. Use DGEXPV instead.
*                DMEXPV/DGEXPV require the matrix-vector product 
*                y = A*x = Q'*x, i.e, the TRANSPOSE of Q times a vector.
*                Failure to remember this leads to wrong results.
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested acurracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H     wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*              IMPORTANT: DMEXPV requires the product y = Ax = Q'x, i.e.
*              the TRANSPOSE of the transition rate matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*     Markov chains are usually well-conditioned problems.
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1,
     .                 roundoff, s_round, x_round

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
c     2019-11-04_NJM
      double precision DDOTX, DNRM2X, DASUMX, w1
c      double precision DDOTX, DNRM2X, DASUMX

*---  check restrictions on input parameters ...
      iflag = 0
*      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) stop 'iflag = -1'
*      if ( liwsp.lt.m+2 ) stop 'iflag = -2'
*      if ( m.ge.n .or. m.le.0 ) stop 'iflag = -3'
       if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
       if ( liwsp.lt.m+2 ) iflag = -2
       if ( m.ge.n .or. m.le.0 ) iflag = -3

*      if ( iflag.ne.0 ) stop 'bad sizes input DMEXPV njm2'
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DMEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      
*     Starting point for finding H, the transition matrix, in wsp
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      sgn      = SIGN( 1.0d0,t )
      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      s_round  = 0.0d0
      x_error  = 0.0d0
      x_round  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      call DCOPYX( n, v,1, w,1 )
      beta = DNRM2X( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call mydgcoov( wsp(j1v-n), wsp(j1v) , n , nz, ia, ja, a)
         do i = 1,j
            hij = DDOTX( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPX( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2X( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            ireject = ireject + 0
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCALX( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call mydgcoov( wsp(j1v-n), wsp(j1v) , n , nz, ia, ja, a)
      avnorm = DNRM2X( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0

c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402
      if ( ireject.eq.1 ) then
         goto 402
      endif
      
 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ..
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

 402  continue
* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            ireject = ireject + 0
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            ireject = ireject + 0
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
*	Original:

      call DGEMX( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )

*******************
*  DGEMX  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*******************
* This says the equivalent would be:
*
* http://computer-programming-forum.com/49-fortran/fe18020d2b9fe2d3.htm
* 
* mathematically, the result of (alpha*y'*A) is the transpose of
* (alpha*A'*y) but as the results are usually stored in 1-D arrays in
* fortran, there is no practical difference (you *cannot* tell the
* difference between a 1-d row vector and a 1-d column vector in fortran
* - they're just 1-d arrays). so blas can do what you want - however as
* your Q is symmetric this is unnecessary in your case anyway...
* 
*		Maybe we can reverse the inputs to switch alpha*A*y --> alpha*y'*A'
*		or something...nah, none of this works...
*
*	Nope:
*      call DGEMX( 'T', n,mx,beta,wsp(iexph),n,wsp(iv),1,0.0d0,w,1 )
*      call DGEMX( 'T', 1,n,beta,wsp(iexph),n,wsp(iv),1,0.0d0,w,1 )
*      call DGEMX( 'T', 1,n,beta,wsp(iexph),n,wsp(iv),n,0.0d0,w,1 )
*      call DGEMX( 'T', n,1,beta,wsp(iexph),n,wsp(iv),n,0.0d0,w,1 )
*      call DGEMX( 'T', n,1,beta,wsp(iexph),n,wsp(iv),1,0.0d0,w,1 )
*            call DGEMX( 'n', n,1,beta,wsp(iexph),n,wsp(iv),n,0.0d0,w,1 )
*            call DGEMX( 'n', n,1,beta,wsp(iexph),n,wsp(iv),1,0.0d0,w,1 )
*           call DGEMX( 'T', 1,n,beta,wsp(iexph),1,wsp(iv),1,0.0d0,w,1 )
*           call DGEMX( 'T', 1,n,beta,wsp(iexph),1,wsp(iv),n,0.0d0,w,1 )

*	Nope:
*      call DGEMX( 'T', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
*      call DGEMX( 'c', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )

*      call DGEMX( 'n', n,mx,wsp(iexph),wsp(iv),n,beta,1,0.0d0,w,1 )
*      call DGEMX( 'T', n,mx,wsp(iexph),wsp(iv),n,beta,1,0.0d0,w,1 )
*		...and others...
      beta = DNRM2X( n, w,1 )
      hump = MAX( hump, beta )
*
*---  Markov model constraints ...
*
      j = 0
      do i = 1,n
         if ( w(i).lt.0.0d0 ) then
            w(i) = 0.0d0
            j = j + 1
         endif
      enddo
      p1 = DASUMX( n, w,1 )

c 2019-10-08
c      if ( j.gt.0 ) call DSCALX( n, 1.0d0/p1, w,1 )

c 2019-11-04_NJM
      w1 = w(1)
      if ( j.gt.0 ) call DSCALX( n, 1.0d0/p1, w1,1 )
c      if ( j.gt.0 ) call DSCALX( n, 1.0d0/p1, w,1 )
      roundoff = DABS( 1.0d0-p1 ) / DBLE( n )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc, roundoff )
      err_loc = MAX( err_loc, rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         ireject = ireject + 0
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      s_round = s_round + roundoff
      x_error = MAX( x_error, err_loc )
      x_round = MAX( x_round, roundoff )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = x_round
      wsp(4)  = s_round
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
      END


*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

      implicit none
      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a general matrix in
*     full, using the irreducible rational Pade approximation to the 
*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degree of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*                  
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr), 
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
*     and set scale = t/2^ns ...
*
* set matrix to 0
      do i = 1,m
         wsp(i) = 0.0d0
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      enddo
      hnorm = ABS( t*hnorm )
*      hnorm = t
* njm1
*      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of DGPADM.'
*      if ( hnorm.eq.0.0d0 ) hnorm=t/2

* This error may happen with DMEXPV, with equal starting probabilities in v
* and a non-symmetrical matrix
*      if ( hnorm.eq.0.0d0 ) stop 'NJMerr1-nullH DMEXPVmbe=inprobs'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
*
*---  H2 = scale2*H*H ...
*
      call DGEXX( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEXX( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEXX( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEXX( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPX( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
*      if ( iflag.ne.0 ) stop 'Problem in DGESV (within DGPADM)'
      
      call DSCALX( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCALX( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEXX( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

      implicit none
      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a symmetric matrix
*     in full, using the irreducible rational Pade approximation to the 
*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix (both lower and upper parts).
*
*     t         : (input) time-scale (can be < 0).
*                  
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr), 
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DSPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2; 
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = 0.0d0
      enddo
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      enddo
      hnorm = ABS( t*hnorm )
*      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of DSPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
*
*---  H2 = scale2*H*H ...
*
      call DGEXX( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEXX( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      enddo
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEXX( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEXX( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPX( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DSYSV( 'U',m,m,wsp(iq),m,ipiv,wsp(ip),m,wsp(ih2),mm,iflag )
*      if ( iflag.ne.0 ) stop 'Problem in DSYSV (within DSPADM)'
      call DSCALX( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      enddo
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCALX( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEXX( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      enddo
 200  continue
      iexph = iput
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DGCHBV( m, t, H,ldh, y, wsp, iwsp, iflag )

      implicit none
      integer          m, ldh, iflag, iwsp(m)
      double precision t, H(ldh,m), y(m)
*      complex(kind=8)       wsp(m*(m+2))
      complex(8)       wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DGCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is a General matrix.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) argument matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double precision alpha0
*     complex(kind=8) alpha(ndeg), theta(ndeg)
      complex(8) alpha(ndeg), theta(ndeg)

      intrinsic DBLE
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         enddo
         call ZGESV( M, 1, WSP(iH),M, IWSP, WSP(iY),M, IFLAG )
*         if ( IFLAG.ne.0 ) stop 'Error in DGCHBV'
*---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSCHBV( m, t, H,ldh, y, wsp, iwsp, iflag )

      implicit none
      integer          m, ldh, iflag, iwsp(m)
      double precision t, H(ldh,m), y(m)
*      complex(kind=8)       wsp(m*(m+2))
      complex(8)       wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DSCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be symmetric.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) symmetric matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double precision alpha0
*      complex(kind=8) alpha(ndeg), theta(ndeg), w
      complex(8) alpha(ndeg), theta(ndeg), w

      intrinsic ABS,CMPLX,DBLE,MIN
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         enddo
         call ZSYSV('U', M, 1, WSP(iH),M, IWSP, WSP(iY),M, W,1, IFLAG )
*         if ( IFLAG.ne.0 ) stop 'Error in DSCHBV'
*---     Accumulate the partial result in y ...     
         do i = 1,m
            y(i) = y(i) + DBLE( alpha(ip)*wsp(iy+i-1) )
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DNCHBV( m, t, H,ldh, y, wsp )

      implicit none
      integer          m, ldh
c     2019-07-02_NJM:
c     double precision t, H(ldh,m), y(m), wsp(m*(m+2))
      double precision t, H(ldh,m), y(m), wsp(m*(m+2))
*      complex(kind=8) wspc,wspd,wspe,wspf
      complex(8) wspc,wspd,wspe,wspf

*-----Purpose----------------------------------------------------------|
*
*---  DNCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise. 
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
*      complex(kind=8) ZERO
      complex(8) ZERO
c     2019-07-02_NJM:
c     integer ndeg, i, j, k, ip, ih, iy, iz, tempn
      integer ndeg, i, j, k, ip, ih, iy, iz, tempn
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
*      complex(kind=8) alpha(ndeg), theta(ndeg), tmpc, wspee, wspff
      complex(8) alpha(ndeg), theta(ndeg), tmpc, wspee, wspff

      intrinsic ABS,DBLE,MIN
      
*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*     
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      enddo
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            enddo
c            2018-09-30_NJM:
c            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
       wsp(ih+(j-1)*m+j-1)=REAL(wsp(ih+(j-1)*m+j-1)-theta(ip),KIND=8)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            enddo
         enddo
         do i = 1,m-1
*---     Get pivot and exchange rows ...
c        2019-07-01_NJM: putting 1st ZSWAP on 1 row
         if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
c         2019-07-02_NJM:
c         call ZSWAPX(m-i+1,wsp(ih+(i-1)*m+i-1),m,wsp(ih+(i-1)*m+i),m)
          tempn = m-i+1
c /usr/local/gcc10/bin/gfortran -fno-optimize-sibling-calls
c -fpic  -g -O2 -mtune=native -Wall -fallow-argument-mismatch
c -c my_expokit.f -o my_expokit.o
c my_expokit.f:1157:26:
c 1156 |           call ZSWAPX(tempn,wspc,m,wspd,m)
c      |                            2
c 1157 |           call ZSWAPX( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
c      |                          1
c Warning: Type mismatch between actual argument at (1) and 
c actual argument at (2) (REAL(8)/COMPLEX(8)).

* 2023-06-21 error:
* wspc = complex(wsp(ih+(i-1)*m+i-1),0)
* where complex is a GNU extension not supported by the flang compiler:
* the F77 standard function is CMPLX.
* error: Semantic errors in my_expokit.f
* ./my_expokit.f:1174:18: error: No explicit type declared for 'complex'
*             wspc = complex(wsp(ih+(i-1)*m+i-1),0)

*          wspc = complex(wsp(ih+(i-1)*m+i-1),0)
*          wspd = complex(wsp(ih+(i-1)*m+i),0)

          wspc = CMPLX(wsp(ih+(i-1)*m+i-1),0,KIND=8)
          wspd = CMPLX(wsp(ih+(i-1)*m+i),0,KIND=8)
          call ZSWAPX(tempn,wspc,m,wspd,m)
c         call ZSWAPX( 1, wsp(iy+i-1),1, wsp(iy+i),1 )

* 2023-06-21
*          wspee = complex(wsp(iy+i-1),0)
*          wspff = complex(wsp(iy+i),0) 
          wspee = CMPLX(wsp(iy+i-1),0,KIND=8)
          wspff = CMPLX(wsp(iy+i),0,KIND=8) 
          call ZSWAPX( 1, wspee,1, wspff,1 )
          
         endif
*---     Forward eliminiation ... 
c        2019-07-02_NJM:
c        tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
         tmpc = -1*(wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1))
c        call ZAXPX(m-i,-tmpc,wsp(ih+i*m+i-1),m,wsp(ih+i*m+i),m )
c        2019-07-02_NJM:

* 2023-06-21:
*         wspe = complex(wsp(ih+i*m+i-1),0)
*         wspf = complex(wsp(ih+i*m+i),0)
         wspe = CMPLX(wsp(ih+i*m+i-1),0,KIND=8)
         wspf = CMPLX(wsp(ih+i*m+i),0,KIND=8)
         call ZAXPX(m-i,tmpc,wspe,m,wspf,m )
c         2018-09-30_NJM:
c         wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
          wsp(iy+i) = REAL( wsp(iy+i) - tmpc*wsp(iy+i-1), KIND=8 )
         enddo
*---     Backward substitution ...    
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            enddo
c           2018-09-30_NJM:
c            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
            wsp(iy+i-1) = REAL( tmpc / wsp(ih+(i-1)*m+i-1), KIND=8 )
         enddo
*---     Accumulate the partial result in y ...     
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         enddo
      enddo
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine mydgcoov ( x, y , n , nz, ia, ja, a)
      implicit none
      double precision x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed here to be under the COOrdinates storage format.
*
      integer n, nz
      integer ia(nz), ja(nz)
      double precision a(nz)
      integer i, j
 
      do j = 1,n
         y(j) = 0.0d0
      enddo
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      enddo
      END










* myDGEXPV:      
      subroutine myDGEXPV( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp,iwsp,liwsp,itrace,iflag,ia,ja,a,nz )

      implicit none
      integer n,nz,m,lwsp,liwsp,itrace,iflag,iwsp(liwsp),ia(nz),ja(nz)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp), a(nz)

*-----Purpose----------------------------------------------------------|
*
*---  DGEXPV computes w = exp(t*A)*v - for a General matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*                      
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in wsp and iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 1000,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOTX, DNRM2X

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call DCOPYX( n, v,1, w,1 )
      beta = DNRM2X( n, w,1 )
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call mydgcoov( wsp(j1v-n), wsp(j1v), n , nz, ia, ja, a)
         do i = 1,j
            hij = DDOTX( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPX( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2X( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCALX( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call mydgcoov( wsp(j1v-n), wsp(j1v), n , nz, ia, ja, a )
      avnorm = DNRM2X( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402
      if ( ireject.eq.1 ) then
         goto 402
      endif

 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

 402  continue
* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            ireject = ireject + 0
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            ireject = ireject + 0
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMX( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2X( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         ireject = ireject + 0
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSEXPV( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

      implicit none
      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp)
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DSEXPV computes w = exp(t*A)*v - for a Symmetric matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*                      
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, nbr of Tridiagonal matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, nbr of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hjj, hump, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOTX, DNRM2X

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DSEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call DCOPYX( n, v,1, w,1 )
      beta = DNRM2X( n, w,1 )
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      enddo
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         if ( j.gt.1 )
     .     call DAXPX(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = DDOTX( n, wsp(j1v-n),1, wsp(j1v),1 )
         call DAXPX( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DNRM2X( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            ireject = ireject + 0
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         wsp(ih+j*mh+j-1) = hj1j
         call DSCALX( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2X( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m-1) = 0.0d0
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402
      if ( ireject.eq.1 ) then
         goto 402
      endif

      
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         enddo
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue
* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            ireject = ireject + 0
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            ireject = ireject + 0
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMX( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2X( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         ireject = ireject + 0
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DGPHIV( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag ) 

      implicit none
      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
c     2019-11-04_NJM
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp), u1
c      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp)
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a General matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*   
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 1000, 
     .           mxreject = 0,
     .           ideg     = 6, 
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hij, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOTX, DNRM2X

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm
 
      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call DCOPYX( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv) )
c     2019-10-08
c     call DAXPX( n, 1.0d0, u,1, wsp(iv),1 )
c     2019-11-04_NJM
c      call DAXPX( n, 1.0d0, u(:),1, wsp(iv),1 )
      u1 = u(1)
      call DAXPX( n, 1.0d0, u1,1, wsp(iv),1 )
      beta = DNRM2X( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call DSCALX( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = DDOTX( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPX( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = DNRM2X( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            ireject = ireject + 0
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCALX( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2X( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = 1.0d0
      wsp(ih+(m-1)*mh+m) = 0.0d0
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = 1.0d0
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402
      if ( ireject.eq.1 ) then
         goto 402
      endif

 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
*---  irreducible rational Pade approximation ...
      mx = mbrkdwn + MAX( 1,k1 )
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .             wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      
c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402 (which uses 401)
      hj1j = 0.0d0
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)
 
 402  continue
*---  error estimate ...
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif

*---  reject the step-size if the error is not acceptable ...   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and. 
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            ireject = ireject + 0
         endif 
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            ireject = ireject + 0
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMX( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1 

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step 
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         ireject = ireject + 0
      endif
 
      step_min = MIN( step_min, t_step ) 
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )
 
      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1
 
 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      END
*----------------------------------------------------------------------|
*----------------------------------------------------------------------|
      subroutine DSPHIV( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag ) 

      implicit none
      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
c     2019-11-04_NJM
c      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp)
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp)
      double precision u1
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DSPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution 
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a Symmetric matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*                      
*     m      : (input) maximum size for the Krylov basis.
*                      
*     t      : (input) time at wich the solution is needed (can be < 0).
*   
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*  
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u 
* 
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500, 
     .           mxreject = 0,
     .           ideg     = 6, 
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hjj, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOTX, DNRM2X

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
*      if ( iflag.ne.0 ) stop 'bad sizes (in input of DSPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1 
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm
 
      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call DCOPYX( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv) )

c     2019-11-04_NJM
c      call DAXPX( n, 1.0d0, u,1, wsp(iv),1 )
      u1 = u(1)
      call DAXPX( n, 1.0d0, u1,1, wsp(iv),1 )
      beta = DNRM2X( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call DSCALX( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      enddo

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         if ( j.gt.1 )
     .     call DAXPX(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = DDOTX( n, wsp(j1v-n),1, wsp(j1v),1 )
         call DAXPX( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DNRM2X( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            ireject = ireject + 0
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         wsp(ih+j*mh+j-1) = hj1j
         call DSCALX( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2X( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = 1.0d0
      wsp(ih+m*mh+m-1)   = 0.0d0
      wsp(ih+(m-1)*mh+m) = 0.0d0
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = 1.0d0
      enddo
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402
      if ( ireject.eq.1 ) then
         goto 402
      endif

 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .              wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns

c     2018-09-30_NJM: This will never happen, as ireject=0
c     But, satifies need to use 402 (which uses 401)
      hj1j = 0.0d0

      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)
 
 402  continue
* 
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm 
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and. 
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            ireject = ireject + 0
         endif 
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            ireject = ireject + 0
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMX( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1 

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step 
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         ireject = ireject + 0
      endif

      step_min = MIN( step_min, t_step ) 
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )
 
      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1
 
 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      END
*----------------------------------------------------------------------|
