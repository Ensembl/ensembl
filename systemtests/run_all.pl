#!/usr/local/bin/perl

use strict;

system ("perl dna_test.pl -dbname ens100");
system ("perl every_atleast.pl -dbname ens100");
system ("perl stop_codons.pl -dbname ens100");
system ("perl transcript_strand.pl -dbname ens100");
system ("perl exon_duplicates.pl -dbname ens100");
