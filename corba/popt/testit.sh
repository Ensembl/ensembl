#!/bin/sh

testname=$1
case "$testname" in
*/*) ;;
*) testname=./$testname ;;
esac

run() {
    prog=$1; shift
    name=$1; shift
    answer=$1; shift

    echo Running test $name.

    result=`$prog $*`
    if [ "$answer" != "$result" ]; then
	echo "Test \"$*\" failed with: $result"
	exit 2
    fi
}

run ${testname} "test1 - 1" "arg1: 1 arg2: (none)" --arg1
run ${testname} "test1 - 2" "arg1: 0 arg2: foo" --arg2 foo
run ${testname} "test1 - 3" "arg1: 1 arg2: something" --arg1 --arg2 something
run ${testname} "test1 - 4" "arg1: 0 arg2: another" --simple another
run ${testname} "test1 - 5" "arg1: 1 arg2: alias" --two
run ${testname} "test1 - 6" "arg1: 1 arg2: (none) rest: --arg2" --arg1 -- --arg2 
run ${testname} "test1 - 7" "arg1: 0 arg2: abcd rest: --arg1" --simple abcd -- --arg1 
run ${testname} "test1 - 8" "arg1: 1 arg2: (none) rest: --arg2" --arg1 --takerest --arg2 
run ${testname} "test1 - 9" "arg1: 0 arg2: foo" -2 foo
run ${testname} "test1 - 10" "arg1: 0 arg2: (none) arg3: 50" -3 50
run ${testname} "test1 - 11" "arg1: 0 arg2: bar" -T bar
run ${testname} "test1 - 12" "arg1: 1 arg2: (none)" -O 
run ${testname} "test1 - 13" "arg1: 1 arg2: foo" -OT foo
run ${testname} "test1 - 14" "arg1: 0 arg2: (none) inc: 1" --inc
run ${testname} "test1 - 15" "arg1: 0 arg2: foo inc: 1" -i --arg2 foo
POSIX_ME_HARDER=1 run ${testname} "test1 - 16" "arg1: 1 arg2: (none) rest: foo --arg2 something" --arg1 foo --arg2 something
POSIXLY_CORRECT=1 run ${testname} "test1 - 17" "arg1: 1 arg2: (none) rest: foo --arg2 something" --arg1 foo --arg2 something
run ${testname} "test1 - 18" "callback: c sampledata bar arg1: 1 arg2: (none)" --arg1 --cb bar
run ${testname} "test1 - 19" "${testname} ;" --echo-args
run ${testname} "test1 - 20" "${testname} ; --arg1" --echo-args --arg1
run ${testname} "test1 - 21" "${testname} ; --arg2 something" -T something -e
run ${testname} "test1 - 22" "${testname} ; --arg2 something -- more args" -T something -a more args
run ${testname} "test1 - 23" "${testname} ; --echo-args -a" --echo-args -e -a
run ${testname} "test1 - 24" "arg1: 0 arg2: (none) short: 1" -shortoption
run ${testname} "test1 - 25" "arg1: 0 arg2: (none) short: 1" --shortoption
run ${testname} "test1 - 26" "callback: c arg for cb2 foo arg1: 0 arg2: (none)" --cb2 foo
run ${testname} "test1 - 27" "arg1: 0 arg2: (none) -" -
run ${testname} "test1 - 28" "arg1: 0 arg2: foo -" - -2 foo

echo ""
echo "Passed."

