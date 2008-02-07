#!/bin/sh

gunzip -c *.xml.gz | grep '<entry_count>' | perl -e 'while (<>) {$_ =~ /\D+(\d+)\D+/; $t+=$1;} print "$t\n";'