### Checks for memory leaks

set.seed(1234)


shutup_and_time <- function(s)
{
  capture.output(suppressWarnings(t <- system.time(s)))
  
  return( t[3] )
}


test <- function(fun)
{
  x <- NULL
  
  suppressPackageStartupMessages(library(rexpokit))
  pkg_fun <- eval(parse(text=fun))
  shutup_and_time(A <- pkg_fun(Qmat=x))
  
  invisible()
}


tests <- c("expokit_dgpadm_Qmat", 
           "expokit_dmexpv_Qmat",
           "expokit_wrapalldmexpv_tvals",
           "expokit_dgexpv_Qmat",
           "expokit_wrapalldgexpv_tvals")


ignore <- lapply(tests, test)
