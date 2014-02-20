set.seed(1234)
n <- 5

x <- matrix(rnorm(n*n), n, n)


test <- function(mat, fun)
{
  suppressPackageStartupMessages(library(rexpokit))
  pkg_fun <- eval(parse(text=fun))
  A <- pkg_fun(Qmat=x)
  detach("package:rexpokit", unload=TRUE)
  library.dynam.unload("rexpokit", system.file(package = "rexpokit"))
  
  suppressPackageStartupMessages(library(rexpokit, lib.loc="~/tmp/tmp"))
  pkg_fun <- eval(parse(text=fun))
  B <- pkg_fun(Qmat=x)
  
  
  cat(paste(fun, ":  ", all.equal(A, B), "\n", sep=""))
}

