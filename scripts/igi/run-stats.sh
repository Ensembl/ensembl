#!/bin/sh -x

# wrapper to rename the sources:
./stats-from-merge-files.pl "$@" | ./remap-sources.sed
