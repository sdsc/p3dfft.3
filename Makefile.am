EXTRA_DIST = README.md AUTHORS _config.yml build/init.C build/exec.C build/plan.C build/templ.C build/wrap.C

include_HEADERS = include/p3dfft.h

SUBDIRS = build# sample

all-local:
	@echo $(SUCCESS)
clean-local:
	-[ -z "include/p3dfft.mod" ] || rm -f include/p3dfft.mod
