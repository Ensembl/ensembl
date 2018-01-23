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
use Test::Exception;
use File::Temp;
use Bio::EnsEMBL::Utils::IO qw/:all/;
use IO::String;

my ($tmp_fh, $file) = File::Temp::tempfile();
close($tmp_fh);
unlink $file;


my $contents = <<'EOF';
>X
AAAAGGGTTCCC
TTGGCCAAAAAA
ATTC
EOF

my $expected_array = [qw/>X AAAAGGGTTCCC TTGGCCAAAAAA ATTC/];

{
  throws_ok { slurp($file) } qr/No such file/, 'File does not currently exist so die';
  
  work_with_file($file, 'w', sub {
    my ($fh) = @_;
    print $fh $contents;
    return;
  });
  
  my $written_contents = slurp($file);
  is($contents, $written_contents, 'Contents should be the same');
  
  my $written_contents_ref = slurp($file, 1);
  is('SCALAR', ref($written_contents_ref), 'Asked for a ref so expect one back');
  is($contents, $$written_contents_ref, 'Contents should be the same');
  
  work_with_file($file, 'r', sub {
    my ($fh) = @_;
    my $line = <$fh>;
    chomp($line);
    is($line, '>X', 'First line expected to be FASTA header'); 
  });
  
  my $chomp = 1;
  is_deeply(slurp_to_array($file, $chomp), $expected_array, 'Checking slurp to array with chomp');
  $chomp = 0;
  is_deeply(slurp_to_array($file, $chomp), [ map { $_."\n" } @{$expected_array}], 'Checking slurp to array with chomp');

  my $iterator_counter = 0;  
  iterate_file($file, sub {
    my ($line) = @_;
    chomp($line);
    is($line, $expected_array->[$iterator_counter++], sprintf('Checking line %d is ok', $iterator_counter+1));
    return;
  });
  
  unlink $file;
  
  dies_ok { slurp($file) } 'File no longer exists so die';

  spurt($file, $contents);
  my $rewritten_contents = slurp($file);
  is($contents, $rewritten_contents, 'Contents should still be the same');

  spurt($file, $contents, 'append');
  $rewritten_contents = slurp($file);
  is($contents.$contents, $rewritten_contents, 'Contents should be doubled');
  unlink $file;

}

{
  my $content = 'ABCDE1198473' x 33012;
  my $src = IO::String->new($content);
  
  {
    my $trg = IO::String->new();
    move_data($src, $trg);
    is(${$trg->string_ref()}, $content, 'Checking copied data is as expected');
  }
  
  {
    $src->setpos(0);
    my $trg = IO::String->new();
    move_data($src, $trg, (8*1024*1024)); #8MB
    is(${$trg->string_ref()}, $content, 'Checking large buffer copied data is as expected');
  }
}

{
  gz_work_with_file($file, 'w', sub {
    my ($fh) = @_;
    print $fh $contents;
    return;
  });
  
  my $written_contents = gz_slurp($file);
  is($contents, $written_contents, 'Gzipped Contents should be the same');
  my $non_gz_written_contents = slurp($file);
  isnt($contents, $non_gz_written_contents, 'Reading normally should not return the same contents');
  
  my $chomp = 1;
  is_deeply(gz_slurp_to_array($file, $chomp), $expected_array, 'Checking slurp to array with chomp');
  $chomp = 0;
  is_deeply(gz_slurp_to_array($file, $chomp), [ map { $_."\n" } @{$expected_array}], 'Checking slurp to array with chomp');
  
  unlink $file;
  
  dies_ok { slurp($file) } 'File no longer exists so die';
}

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
    4 unless $BZIP2_OK;

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
    4 unless $ZIP_OK;

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
