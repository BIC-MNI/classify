SUBDIRS = ann bayes fcm hcm knn mindist scripts models cleantag

EXTRA_DIST = class_globals.h \
	classify.h \
	class_protos.h \
	unsuper_globals.h

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
