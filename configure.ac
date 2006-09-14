dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([classify.cc])
AM_INIT_AUTOMAKE(classify, 1.0.04)
AM_CONFIG_HEADER(config.h)

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

dnl Build only static libs by default
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

smr_WITH_BUILD_PATH

# Check for required C libraries
mni_REQUIRE_VOLUMEIO
mni_REQUIRE_EBTKS

# Find the version of perl in use
AC_PATH_PROGS(PERL, perl)

AC_CONFIG_FILES([Makefile ann/Makefile bayes/Makefile fcm/Makefile
hcm/Makefile knn/Makefile mindist/Makefile scripts/Makefile models/Makefile
cleantag/Makefile])
AC_OUTPUT(
epm-header
scripts/classify_clean
scripts/insect
)
