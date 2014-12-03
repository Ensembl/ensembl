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
use Bio::EnsEMBL::Utils::IO qw/:all/;

ok(1, 'Module compiles');

my $tmpdir = File::Temp->newdir();
my $dirname = $tmpdir->dirname;

#
# test working with bzip2 and zip files
#
# create a dummy file
my $tmpfile = File::Temp->new(DIR => $dirname, SUFFIX => '.txt');
my $tmpfilename = $tmpfile->filename;
print $tmpfile "test data\n";
print $tmpfile "some more data.";
$tmpfile->close;

my $BZIP2_OK = 0;
eval {
  require IO::Compress::Bzip2;
  require IO::Uncompress::Bunzip2;
  $BZIP2_OK = 1;
};

SKIP: {
  skip "Cannot run Bzip2/Bunzip2 tests, install related IO::[Un]Compress modules first",
    2 unless $BZIP2_OK;

  # send the content of the tmpfile to another
  # bzip2 compressed file
  my $file_content = slurp($tmpfilename);
  my $bz2tmpfile = File::Temp->new(DIR => $dirname, SUFFIX => '.bz2');

  bz_work_with_file($bz2tmpfile->filename, 'w', sub {
    my ($fh) = @_;
    print $fh $file_content;
    return;
  });

  # check the content of the compressed file
  my $content = bz_slurp($bz2tmpfile->filename);
  like($content, qr/test data/, "Bzip2: correct content");
  like($content, qr/more data/, "Bzip2: more correct content");

  # test bz_slurp_to_array
  my $content_array = bz_slurp_to_array($bz2tmpfile->filename, 1);
  ok($content_array->[0] eq 'test data', "Bzip2 slurped file first element");
  ok($content_array->[1] eq 'some more data.', "Bzip2 slurped file second element");

}

my $ZIP_OK = 0;
eval {
  require IO::Compress::Zip;
  require IO::Uncompress::Unzip;
  $ZIP_OK = 1;
};

SKIP: {
  skip "Cannot run Zip/Unzip tests, install related IO::[Un]Compress modules first",
    2 unless $ZIP_OK;

  # send the content of the tmpfile to another
  # bzip2 compressed file
  my $file_content = slurp($tmpfilename);
  my $ziptmpfile = File::Temp->new(DIR => $dirname, SUFFIX => '.zip');

  zip_work_with_file($ziptmpfile->filename, 'w', sub {
    my ($fh) = @_;
    print $fh $file_content;
    return;
  });

  # check the content of the compressed file
  my $content = zip_slurp($ziptmpfile->filename);
  like($content, qr/test data/, "Zip: correct content");
  like($content, qr/more data/, "Zip: more correct content");

  # test zip_slurp_to_array
  my $content_array = zip_slurp_to_array($ziptmpfile->filename, 1);
  ok($content_array->[0] eq 'test data', "Zip slurped file first element");
  ok($content_array->[1] eq 'some more data.', "Zip slurped file second element");

}


#
# test filtering the content of a directory
#
my $perl_tmp_file1 = File::Temp->new(DIR => $dirname, SUFFIX => '.pl');
my $perl_tmp_file2 = File::Temp->new(DIR => $dirname, SUFFIX => '.pl');
my $perl_tmp_file3 = File::Temp->new(DIR => $dirname, SUFFIX => '.pl');
my $other_tmp_file1 = File::Temp->new(DIR => $dirname, SUFFIX => '.dat');
my $other_tmp_file2 = File::Temp->new(DIR => $dirname, SUFFIX => '.dat');

is(scalar @{Bio::EnsEMBL::Utils::IO::filter_dir($dirname, sub {
						  my $file = shift;
						  return $file if $file =~ /\.pl$/;
						})}, 3, "Number of entries in filtered and sorted dir");

done_testing();
