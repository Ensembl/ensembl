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


use warnings;
use strict;

sub BEGIN {
  use FindBin qw/$Bin/;
  use lib $Bin;
  use EnsEMBL::PerlFilter;
}

use Getopt::Long;
use Pod::Usage;
my $opts = {
  help => 0,
  verbose => 0
};
GetOptions( $opts, qw/verbose help/) or pod2usage(-msg => 'Error during command line parsing', -exitlevel => 1, -verbose => 1);
pod2usage(-exitlevel => 0, -verbose => 2) if $opts->{help};
my $file = $ARGV[0];
pod2usage(-msg => 'No file given', -exitlevel => 1, -verbose => 1) if ! $file;
pod2usage(-msg => 'File '.$file.' does not exist', -exitlevel => 1, -verbose => 1) if ! -f $file;
my ($ext) = lc($ARGV[0]) =~ /\.([a-z]+)$/;

if( $ext eq 'pl' || $ext eq 'pm' || $ext eq 'perl' ) {
  print STDERR "Working with a Perl file\n" if $opts->{verbose};
  my $filter = EnsEMBL::PerlFilter->new(\*STDOUT);
  $filter->filter($file);
}
else {
  print STDERR "Passing file through\n" if $opts->{verbose};
  print <>;
}

exit 0;
__END__
=pod

=head1 NAME

ensembldoxygenfilter.pl

=head1 SYNOPSIS

  ./ensembldoxygenfilter.pl -v modules/Bio/EnsEMBL/Registry.pm

  ./ensembldoxygenfilter.pl -help

=head1 DESCRIPTION

Generate Doxygen compatible filtered files from Ensembl POD. Bring this onto
your PATH and edit your Doxygen configuration file and specify 

INPUT_FILTER=ensembldoxygenfilter.pl

The code has been flagged as executable in CVS and so should be self
executing once you check it out.

=head1 OPTIONS
 
=over 8 

=item B<--verbose>

Prints messages to STDERR

=item B<--help>

Prints this message

=back
