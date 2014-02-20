set.seed(1234)

shutup_and_time <- function(s)
{
  capture.output(suppressWarnings(t <- system.time(s)))
  
  return( t[3] )
}

test <- function(mat, fun, bench=FALSE)
{
  if (bench)
  {
    n <- 50
    x <- matrix(rnorm(n*n), n, n)
  }
  else
    x <- NULL
  
  suppressPackageStartupMessages(library(rexpokit))
  pkg_fun <- eval(parse(text=fun))
  t_new <- shutup_and_time(A <- pkg_fun(Qmat=x))
  detach("package:rexpokit", unload=TRUE)
  library.dynam.unload("rexpokit", system.file(package = "rexpokit"))
  
  
  suppressPackageStartupMessages(library(rexpokit, lib.loc="~/tmp/tmp"))
  pkg_fun <- eval(parse(text=fun))
  t_old <- shutup_and_time(B <- pkg_fun(Qmat=x))
  
  
  cat(paste(fun, ":  ", all.equal(A, B), "\n", sep=""))
  
  if (bench)
    cat(sprintf("Old: %.3f   New: %.3f\n", t_old, t_new))
  
  invisible()
}

