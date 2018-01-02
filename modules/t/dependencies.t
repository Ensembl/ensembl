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

use Test::More;
use Test::Warnings;
use File::Spec;

my ($volume, $directory, $file) = File::Spec->splitpath(__FILE__);
$directory = File::Spec->rel2abs($directory);
my $modules_dir = File::Spec->catdir($directory, File::Spec->updir(), qw/Bio EnsEMBL/);

#test for dependencies on Variation, Compara and Funcgen APIs
my @result = `egrep -r "^use Bio::EnsEMBL::(Variation|Compara|Funcgen){1}" $modules_dir`;

my %result = map{$_ => 1} @result;

my @exceptions = ('/Bio/EnsEMBL/Utils/TranscriptAlleles.pm', 
   	       	 '/Bio/EnsEMBL/Utils/ensembl_init.example');

my $exceptions = join("|",@exceptions);
$exceptions =~ s/\//\\\//g;


foreach my $key (keys %result) {
	if ( $key =~ /($exceptions)/ ) {
	   delete($result{$key});
	}
}

ok(!%result);

if (%result) { 
   warn "Dependencies found in the following files:\n";
   warn keys %result;
}

done_testing();
