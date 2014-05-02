#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use File::Basename;

use Bio::EnsEMBL::Utils::IO::FileFaidx;

my $current_time = time();

my $file = $ARGV[0];
my $faidx = Bio::EnsEMBL::Utils::IO::FileFaidx->new($file, 'write_index');

my $index = $faidx->index_path($file);
if(-f $index) {
  print "Removing the existing index at $index\n";
  unlink $index;
}

print "Generating the index\n";
$faidx->lookup();
my $finished_time = time();
my $elapsed = $finished_time - $current_time;
my $time = gmtime($finished_time);
print "Finished @ $time. Took $elapsed second(s)\n";