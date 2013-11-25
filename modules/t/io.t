# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

################################################ 
#                                              #
# io.t                                         #
#                                              #
# A set of tests to verify various subroutines #
# in the Bio::EnsEMBL::Utils::IO module        #
#                                              #
################################################

use strict;
use warnings;

use Test::More;
use File::Temp;
use Bio::EnsEMBL::Utils::IO qw /:all/;

ok(1, 'module compiles');

#
# test filtering the content of a directory
#
my $tmpdir = File::Temp->newdir();
my $dirname = $tmpdir->dirname;

my $perl_tmp_file1 = File::Temp->new(DIR => $dirname,
				     SUFFIX => '.pl');
my $perl_tmp_file2 = File::Temp->new(DIR => $dirname,
				     SUFFIX => '.pl');
my $perl_tmp_file3 = File::Temp->new(DIR => $dirname,
				     SUFFIX => '.pl');

my $other_tmp_file1 = File::Temp->new(DIR => $dirname,
				      SUFFIX => '.dat');
my $other_tmp_file2 = File::Temp->new(DIR => $dirname,
				      SUFFIX => '.dat');

is(scalar @{Bio::EnsEMBL::Utils::IO::filter_dir($dirname, sub {
						  my $file = shift;
						  return $file if $file =~ /\.pl$/;
						})}, 3, "filter_dir: number of entries in dir");

done_testing();
