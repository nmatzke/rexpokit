set.seed(1234)

randsparse <- function(n, sparsity=.95)
{
  x <- double(n*n)
  for (i in 1L:(n*n))
    x[i] <- sample(c(0.0, rnorm(1)), size=1, prob=c(sparsity, 1-sparsity))
  
  dim(x) <- c(n, n)
  
  x
}

shutup_and_time <- function(s)
{
  capture.output(suppressWarnings(t <- system.time(s)))
  
  return( t[3] )
}

test <- function(mat, fun, bench=FALSE)
{
  if (bench)
  {
    n <- 100
    x <- randsparse(n)
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

