bin_SCRIPTS = insect classify_clean
EXTRA_DIST = insect.pl classify_clean.pl

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

CLEANFILES = insect classify_clean
