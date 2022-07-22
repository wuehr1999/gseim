include mkgseim
include mkfilter_so

all: gsmain libfilter.so

.PHONY: clean install_libfilter
install_libfilter: libfilter.so
	cp libfilter.so ../../gseim/exec/libfilter.so

clean: clean_gsmain clean_libfilter_so
	$(RM) ../../gseim/exec/libfilter.so
