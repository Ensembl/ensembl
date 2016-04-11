DEP_FILE_PATHS := $(wildcard deps/*)
DEP_FILES := $(notdir $(DEP_FILE_PATHS))

all:
	@echo "Not implemented, use target install-dep to install dependencies to PREFIX"

install-dep: perldeps
	for dep in $(DEP_FILES) ; do \
		$(MAKE) -C deps -f $$dep ; \
	done

perldeps:
	cpanm -v --installdeps --notest .
	cpanm -n DBD::SQLite
	cp travisci/MultiTestDB.conf.travisci.mysql  modules/t/MultiTestDB.conf.mysql
	cp travisci/MultiTestDB.conf.travisci.SQLite modules/t/MultiTestDB.conf.SQLite

install:
	@echo "Not implemented, local install only of dependencies to PREFIX"
