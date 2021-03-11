EXTRA_DIST = README.md AUTHORS _config.yml build/init.C build/exec.C build/plan.C build/templ.C build/wrap.C

include_HEADERS = include/p3dfft.h include/p3dfft3config.h include/Fwrap.h include/Cwrap.h sample/FORTRAN/p3dfft_plus_plus.mod sample/FORTRAN/wrap.o sample/FORTRAN/wrap.f90 sample/FORTRAN/p3dfft++.o sample/FORTRAN/p3dfft++.f90

SUBDIRS = build sample

library: p3dfft3config.h
	$(MAKE) -C build

sample: library
	$(MAKE) -C sample

all-local:
	@echo $(SUCCESS)
clean-local:
	-[ -z "include/p3dfft.mod" ] || rm -f include/p3dfft.mod
