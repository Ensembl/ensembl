use strict;
use warnings;

package ErrCode;

use Exporter;

use vars qw(@EXPORT_OK @ISA);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&push_err &pop_err &ec2str);


use constant OK => 0;

use constant PART_EXON_CDS_DELETE_TOO_LONG             => 1;
use constant PART_EXON_CDS_FRAMESHIFT_DELETE_TOO_LONG  => 2;
use constant PART_EXON_CDS_DELETE_ENTIRE               => 3;
use constant PART_EXON_CDS_INSERT_TOO_LONG             => 4;
use constant PART_EXON_CDS_FRAMESHIFT_INSERT_TOO_LONG  => 5;
use constant PART_EXON_CONFUSED                        => 6;

use constant NO_CDS_LEFT => 7;

use constant EXON_DELETE_CODING => 10;
use constant EXON_STRAND_FLIP   => 11;
use constant EXON_INVERT        => 12;
use constant EXON_SCAFFOLD_SPAN => 13;

use constant TRANSCRIPT_STRAND_FLIP   => 20;
use constant TRANSCRIPT_INVERT        => 21;
use constant TRANSCRIPT_SCAFFOLD_SPAN => 22;

use constant TRANSCRIPT_DOES_NOT_TRANSLATE => 30;

use constant PART_EXON_CDS_DELETE_FIVE_PRIME_TOO_LONG  => 40;
use constant PART_EXON_CDS_DELETE_THREE_PRIME_TOO_LONG => 41;
use constant PART_EXON_CDS_DELETE_MIDDLE_TOO_LONG      => 42;


my @ERR_STACK = ();

sub push_err {
  my $err      = shift || return;
  my $err_desc = shift || '';

  push @ERR_STACK, [$err, $err_desc];
}

sub pop_err {
  my ($err, $err_desc);

  return () if(!@ERR_STACK);

  my $arref = pop(@ERR_STACK);
  return @$arref;
}


#converts an error code to a string
sub ec2str {
  my $ec = shift;

  my $str = '';

  if($ec == OK) {
    $str = 'OK';
  }
  elsif($ec == PART_EXON_CDS_DELETE_TOO_LONG) {
    $str =  'Exon CDS deletion too long.';
  }
  elsif($ec == PART_EXON_CDS_FRAMESHIFT_DELETE_TOO_LONG) {
    $str = 'Exon frameshifting deletion too long.';
  }
  elsif($ec == PART_EXON_CDS_DELETE_ENTIRE) {
    $str = 'Entire CDS deleteion.';
  }
  elsif($ec == PART_EXON_CDS_INSERT_TOO_LONG) {
    $str = 'Exon CDS insertion too long.';
  }
  elsif($ec == PART_EXON_CDS_FRAMESHIFT_INSERT_TOO_LONG) {
    $str = 'Exon frameshifting insertion too long.';
  }
  elsif($ec == NO_CDS_LEFT) {
    $str = 'No CDS left following CDS adjustments due to inserts/deletions';
  }
  elsif($ec == PART_EXON_CONFUSED) {
    $str = 'Confused due to insert following short match consumed ' .
      'by introduced frameshift intron.';
  }
  elsif($ec == EXON_DELETE_CODING) {
    $str = 'Entire coding exon deletion.';
  }
  elsif($ec == EXON_STRAND_FLIP) {
    $str = 'Exon flips strands.';
  }
  elsif($ec == EXON_INVERT) {
    $str = 'Exon inversion.';
  }
  elsif($ec == EXON_SCAFFOLD_SPAN) {
    $str = 'Exon spans multiple scaffolds.';
  }
  elsif($ec == TRANSCRIPT_STRAND_FLIP) {
    $str = 'Transcript flips strands.';
  }
  elsif($ec == TRANSCRIPT_INVERT) {
    $str = 'Transcript inversion.';
  }
  elsif($ec == TRANSCRIPT_SCAFFOLD_SPAN) {
    $str = 'Transcript spans multiple scaffolds.';
  }
  elsif($ec == PART_EXON_CDS_DELETE_FIVE_PRIME_TOO_LONG) {
    $str = "CDS deletion at 5' end of exon too long.";
  }
  elsif($ec == PART_EXON_CDS_DELETE_THREE_PRIME_TOO_LONG) {
    $str = "CDS deletion at 3' end of exon too long.";
  }
  elsif($ec == PART_EXON_CDS_DELETE_MIDDLE_TOO_LONG) {
    $str = "CDS deletion in middle of exon too long.";
  }

  return $str;
}

1;
