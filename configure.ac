## autoconf configure.ac && autoreconf -if

AC_INIT(DESCRIPTION)

# Obtain from R.
R_SCMD="${R_HOME}/bin/R CMD config"
FLIBS=`${R_SCMD} FLIBS`
BLAS_LIBS=`${R_SCMD} BLAS_LIBS`
LAPACK_LIBS=`${R_SCMD} LAPACK_LIBS`


# Check tools.
AC_CHECK_PROG([ls_ok], [ls], [yes], [no], ,)
AC_CHECK_PROG([echo_ok], [echo], [yes], [no], ,)
AC_CHECK_PROG([sed_ok], [sed], [yes], [no], ,)
AC_CHECK_PROG([grep_ok], [grep], [yes], [no], ,)


# Check Rlapack.
if test "X$echo_ok" = "Xyes" -a "X$grep_ok" = "Xyes"; then
  Rlapack_ok=`echo ${LAPACK_LIBS} | grep "\-lRlapack"`
  if test "X$Rlapack_ok" = "X"; then
    dnl Rlapack is not found.
    Rlapack_ok="NA"
  else
    echo "Found Rlpack from ${Rlapack_ok}"
  fi
else
  dnl echo and grep are not found.
  Rlapack_ok="NA"
fi

if test "X$Rlapack_ok" = "XNA"; then
  # An external LAPACK is used in run time.
  REXPOKIT_OBJS=""
else
  dnl Rlapack is used.
  
  if test "X$ls_ok" = "Xyes" -a "X$sed_ok" = "Xyes" -a "X$grep_ok" = "Xyes"; then
    ### Check R_ext/Lapack.h
    LAPACK_HEADER="${R_HOME}/include/R_ext/Lapack.h"
  
    funcs=`ls src/expokit/src/*.f | sed -e "s/.*\/\(.*\)\.f/\\1/"`
    REXPOKIT_OBJS=""
    for func in ${funcs}
    do
      func_ok=`grep "F77_NAME($func)" ${LAPACK_HEADER}`
      echo "Found $func from ${func_ok}"
      if test "X${func_ok}" = "X"; then
        REXPOKIT_OBJS="${REXPOKIT_OBJS} expokit/$func.o"
      fi
    done
  else
    ### Use default if all above fail.
    REXPOKIT_OBJS="\$(REXPOKIT_OBJS)"
  fi
fi


AC_SUBST(REXPOKIT_OBJS)
AC_OUTPUT(src/Makevars)

