# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
use File::Temp;
use Test::More;
use Test::Warnings;

#chdir into the file's target & request cwd() which should be fully resolved now.
#then go back
my $file_dir = dirname(__FILE__);
my $original_dir = cwd();
chdir($file_dir);
my $cur_dir = cwd();
chdir($original_dir);
my $root = File::Spec->catdir($cur_dir, File::Spec->updir(),File::Spec->updir());

my $script = File::Spec->catfile($root, "misc-scripts/xref_mapping/xref_config2sql.pl");
my $config = File::Spec->catfile($root, "misc-scripts/xref_mapping/xref_config.ini");

my $tmpdir = File::Temp->newdir();
my $dirname = $tmpdir->dirname;
my $tmpfile = File::Temp->new(DIR => $dirname, SUFFIX => ".txt");
my $cmd = "perl $script $config > $tmpfile";
system($cmd);

my $content = do { local(@ARGV, $/) = $tmpfile; <> };

# Test the file was written in entirety
like($content, qr/FINISHED SUCCESSFULLY/, "Config was parsed successfully");


done_testing();
