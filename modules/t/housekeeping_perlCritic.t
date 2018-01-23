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

use Cwd;
use File::Spec;
use File::Basename qw/dirname/;
use Test::More;
use Test::Warnings;

if ( not $ENV{TEST_AUTHOR} ) {
  my $msg = 'Author test. Set $ENV{TEST_AUTHOR} to a true value to run.';
  plan( skip_all => $msg );
}

eval {
  require Test::Perl::Critic;
  require Perl::Critic::Utils;
};
if($@) {
  plan( skip_all => 'Test::Perl::Critic required.' );
  note $@;
}

#chdir into the file's target & request cwd() which should be fully resolved now.
#then go back
my $file_dir = dirname(__FILE__);
my $original_dir = cwd();
chdir($file_dir);
my $cur_dir = cwd();
chdir($original_dir);
my $root = File::Spec->catdir($cur_dir, File::Spec->updir(),File::Spec->updir());

# Configure critic
Test::Perl::Critic->import(-profile => File::Spec->catfile($root, 'perlcriticrc'), -severity => 5, -verbose => 8);

#Find all files & run
my @perl_files = Perl::Critic::Utils::all_perl_files(
  File::Spec->catdir($root, 'modules')
);
foreach my $perl (@perl_files) {
  critic_ok($perl);
}

done_testing();
