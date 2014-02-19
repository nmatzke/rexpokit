setMethod("expm", signature(x="matrix"),
  function(x, t=2.1)
  {
    ret <- expokit_dgpadm_Qmat(Qmat=x, t=t, transpose_needed=TRUE)
    
    return( ret )
  }
)


