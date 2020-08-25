#!/bin/bash

ENSDIR="${ENSDIR:-$PWD}"

export PERL5LIB=$ENSDIR/bioperl-live:$ENSDIR/ensembl-test/modules:$PWD/modules:$ENSDIR/ensembl-io/modules:$ENSDIR/ensembl-variation/modules:$ENSDIR/ensembl-compara/modules:$PWD/misc-scripts/xref_mapping
export TEST_AUTHOR=$USER

if [ "$DB" = 'mysql' ]; then
    (cd modules/t && ln -sf MultiTestDB.conf.mysql MultiTestDB.conf)
    ln -sf testdb.conf.mysql testdb.conf
elif [ "$DB" = 'sqlite' ]; then
    (cd modules/t && ln -sf MultiTestDB.conf.SQLite MultiTestDB.conf)
    ln -sf testdb.conf.SQLite testdb.conf
    SKIP_TESTS="--skip dbConnection.t,schema.t,schemaPatches.t,strainSlice.t,sliceVariation.t,mappedSliceContainer.t"
else
    echo "Don't know about DB '$DB'"
    exit 1;
fi
ln -sf ../../../modules/t/MultiTestDB.conf misc-scripts/xref_mapping/t/

echo "Running test suite"
rt=0
if [ "$COVERALLS" = 'true' ]; then
  PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl-variation,ensembl-compara' perl $ENSDIR/ensembl-test/scripts/runtests.pl -verbose modules/t $SKIP_TESTS
  rt=$?
  if [ "$DB" = 'mysql' ]; then
    PERL5OPT='-MDevel::Cover=+ignore,bioperl,+ignore,ensembl-test,+ignore,ensembl-variation,ensembl-compara' perl $ENSDIR/ensembl-test/scripts/runtests.pl -verbose misc-scripts/xref_mapping/t
    rt=$(($rt+$?))
  fi
else
  perl $ENSDIR/ensembl-test/scripts/runtests.pl modules/t $SKIP_TESTS
  rt=$?
  if [ "$DB" = 'mysql' ]; then
    perl $ENSDIR/ensembl-test/scripts/runtests.pl misc-scripts/xref_mapping/t
    rt=$(($rt+$?))
  fi
fi

if [ $rt -eq 0 ]; then
  if [ "$COVERALLS" = 'true' ]; then
    echo "Running Devel::Cover coveralls report"
    cover --nosummary -report coveralls
  fi
  exit $?
else
  exit 255
fi
