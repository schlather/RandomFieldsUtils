
AC_INIT([RandomFieldsUtils], 0.5)
CXX=`"${R_HOME}/bin/R" CMD config CXX`
AC_PROG_CXX
AC_LANG(C++)



CPU_FLAGS=""
AC_CHECK_PROG(LSCPU, [lscpu], yes, no)

if test "x${LSCPU}" == xyes; then
  INFO_CPU=`lscpu | grep Flags | tr "[[:upper:]]" "[[:lower:]]"`
else
  AC_CHECK_PROG(SYSCTL, [sysctl], yes, no)
  if test  "x${SYSCTL}" == xyes; then
    INFO_CPU=`sysctl -a 2> /dev/null | grep machdep.cpu.features | tr "[[:upper:]]" "[[:lower:]]"`
  fi
fi


AC_MSG_CHECKING([which SIMD flags are available])
FX=""
## sysctl -a | grep "cpu.features:" # for OS X
for SET in sse2 ssse3 avx2; do 
  CPPFLAGS="-m$SET"
  info=`echo "${INFO_CPU}" | grep " $SET "`
  if test "${info}" != ""; then
     AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]])], [ FX="$FX $CPPFLAGS" ], [ ])
  fi
done

CPPFLAGS="$FX"

##CPPFLAGS="" ## to be deleted if the above does not cause any problems
if test "$CPPFLAGS" != "" ; then
  AC_MSG_RESULT([$CPPFLAGS])
else
  AC_MSG_RESULT([none])
fi

## install.packages("miraculix_0.9.8.tar.gz", configure.args="CXX_FLAGS=-march=native", repos=NULL)
RANDOMFIELDSUTILS_CXXFLAGS="$CXX_FLAGS $CPPFLAGS"
AC_SUBST(RANDOMFIELDSUTILS_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
