AC_INIT(jackFitter,1.0,shultz@jlab.org)
AC_CONFIG_AUX_DIR(config)

AC_LANG_CPLUSPLUS
AC_PROG_CXX(g++ cxx CC)
AC_PROG_CC(gcc cc)
AC_PROG_RANLIB

dnl Initialise automake
AM_INIT_AUTOMAKE
dnl
dnl
dnl --with-xxxx and --enable-xxxx switches 
dnl

dnl --with-adat=DIR
AC_ARG_WITH(adat,
  AC_HELP_STRING(
    [--with-adat=DIR],
    [Building on top of ADAT where ADAT is installed in DIR]
  ),
  [ADAT_HOME="$with_adat"]
)


if test "X${ADAT_HOME}X" = "XX" ; then
  AC_PATH_PROG(ADAT_CONFIG, [adat-config], [])
else
  AC_PATH_PROG(ADAT_CONFIG, [adat-config], [], [${ADAT_HOME}/bin:${PATH}])
fi

if test "X${ADAT_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([ADAT configuration program adat-config not found.])
fi

AC_MSG_NOTICE([Found ADAT configuration program ${ADAT_CONFIG}])
AC_SUBST(CXX, "`${ADAT_CONFIG} --cxx`")
AC_MSG_NOTICE([ADAT compiler: ${CXX}])
AC_SUBST(ADAT_CXXFLAGS, "`${ADAT_CONFIG} --cxxflags`")
AC_MSG_NOTICE([ADAT compile flags: ${ADAT_CXXFLAGS}])
AC_SUBST(ADAT_LDFLAGS,  "`${ADAT_CONFIG} --ldflags`")
AC_MSG_NOTICE([ADAT linking flags: ${ADAT_LDFLAGS}])
AC_SUBST(ADAT_LIBS,     "`${ADAT_CONFIG} --libs`")
AC_MSG_NOTICE([ADAT libraries flags: ${ADAT_LIBS}])

dnl include the m4 macro PAC_ADAT_LINK_CXX_FUNC
m4_include([m4/pac_adat_link_cxx_func.m4])
dnl This function is defined in m4/pac_adat_link_cxx_func.m4
dnl it tries to try to compile the program 
AC_MSG_CHECKING([if we can compile/link a simple adat program])
PAC_ADAT_LINK_CXX_FUNC(
  ${ADAT_CXXFLAGS},
  ${ADAT_LDFLAGS},
  ${ADAT_LIBS},
  ,
  ,
  [adat_link_ok=yes],
  [adat_link_ok=no]
)

if test "X${adat_link_ok}X" = "XyesX" ; then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with adat.
    Use --with-adat=<dir> to select a working version.])
fi 



dnl itpp support
AC_ARG_WITH(itpp, 
  AC_HELP_STRING(
    [--with-itpp=DIR],
    [Use ITPP in DIR]
    ),
    [ITPP_HOME="$with_itpp"],
)

if test "X${ITPP_HOME}X" = "XX"; then
  AC_PATH_PROG(ITPP_CONFIG, [itpp-config], [])
else
  AC_PATH_PROG(ITPP_CONFIG, [itpp-config], [], [${ITPP_HOME}/bin:$PATH] )
fi

if test "X${ITPP_CONFIG}X" = "XX"; then 
  AC_MSG_ERROR([Program itpp-config not found])
fi

AC_MSG_NOTICE([Found ITPP configuration program ${LIBXML2XX_CONFIG}])
AC_SUBST(ITPP_CXXFLAGS, "`${ITPP_CONFIG} --cflags`")
AC_MSG_NOTICE([ITPP compile flags: ${ITPP_CXXFLAGS}])
AC_SUBST(ITPP_LIBS,     "`${ITPP_CONFIG} --static --libs`")
AC_MSG_NOTICE([ITPP libraries flags: ${ITPP_LIBS}])


AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(jackFitter-config)
AC_CONFIG_FILES(lib/Makefile)
AC_CONFIG_FILES(src/Makefile)

AC_OUTPUT