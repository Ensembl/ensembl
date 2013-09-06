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
