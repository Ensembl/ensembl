#!/usr/bin/env perl

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

The code has been flagged as executable in the 

=head1 OPTIONS
 
=over 8 

=item B<--verbose>

Prints messages to STDERR

=item B<--help>

Prints this message

=back