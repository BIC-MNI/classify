SUBDIRS = ann bayes fcm hcm knn mindist scripts models cleantag

EXTRA_DIST = class_globals.h \
	classify.h \
	class_protos.h \
	unsuper_globals.h \
	m4/smr_WITH_BUILD_PATH.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	epm-header.in

AUTOMAKE_OPTIONS = check-news

bin_PROGRAMS = classify

classify_SOURCES = time_stamp.c classify.cc

classify_LDADD = \
	$(top_builddir)/ann/libann.la \
	$(top_builddir)/bayes/libbayes.la \
	$(top_builddir)/fcm/libfcm.la \
	$(top_builddir)/hcm/libhcm.la \
	$(top_builddir)/knn/libknn.la \
	$(top_builddir)/mindist/libmindist.la
