#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

my $old;
my $new;
while (<>) {
    
    if (/grep (\S+)/) {
	$old=$1;
	next;
    }

    if (/(ENSE\S+)\s(ENSE\S+)/) {
	if ($2 eq $old) {
	    $new=$1;
	    if ($old ne $1) {
		print "update translation set start_exon\=\"$1\" where start_exon\=\"$old\"\;\n";
		print "update translation set end_exon\=\"$1\" where end_exon\=\"$old\"\;\n";
		print "update supporting_feature set exon\=\"$1\" where exon\=\"$old\"\;\n";
		print "update exon_transcript set exon\=\"$1\" where exon\=\"$old\"\;\n";
	    }
	}
	elsif ($1 eq $old) {
	    if ($new ne $2) {
		print "update translation set start_exon\=\"$new\" where start_exon\=\"$2\"\;\n";
		print "update translation set end_exon\=\"$new\" where end_exon\=\"$2\"\;\n";
		print "update supporting_feature set exon\=\"$new\" where exon\=\"$2\"\;\n";
		print "update exon_transcript set exon\=\"$new\" where exon\=\"$2\"\;\n";
	    }
	    $old = $2;
	}
    }
}
