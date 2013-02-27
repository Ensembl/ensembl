#!/usr/bin/env perl

package Script;

use strict;
use warnings;

use Bio::DB::Sam;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->check();
  $self->setup();
  $self->process();
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {
    port => 3306
  };
  GetOptions(
    $opts, qw/
      file=s
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub opts {
  my ($self) = @_;
  return $self->{'opts'};
}

sub check {
  my ($self) = @_;
  my $o = $self->opts();

  my @required_params = qw/file/;

  foreach my $r (@required_params) {
    if (!$o->{$r}) {
      pod2usage(
        -message => "-${r} has not been given at the command line but is a required parameter",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  return;
}

sub process {
  my ($self) = @_;
  my $file = $self->opts->{file};
  die "Cannot find file" unless -f $file;
  my $bam = Bio::DB::Bam->open($file);
  my $header = $bam->header();
  my $a = $bam->read1();
  my $query_start  = $a->query->start;
  my $query_end    = $a->query->end;
  my $query_dna    = $a->query->dna;
  printf('%s:%d-%d', $query_dna, $query_start, $query_end);
  print "\n";
  return;
}


Script->run();

1;
__END__

=pod

=head1 NAME

first_bam_alignment.pl

=head1 SYNOPSIS

  #BASIC
  ./first_bam_alignment.pl -file FILE [-help | -man]
  
  #EXAMPLE
  ./first_bam_alignment.pl -file my.bam
  
=head1 DESCRIPTION

A script which returns the first region with BAM data in it. Used as a way of sanity checking that
a BAM has some kind of alignment available.

=head1 OPTIONS

=over 8

=item B<--file>

The location of the BAM

=item B<--help>

Help message

=item B<--man>

Man page

=back

=head1 REQUIREMENTS

=over 8

=item Perl 5.8+

=item Bio::DB::Bam

=back

=end

