use strict;
use warnings;

package StatMsg;

use Exporter;

#
# Status codes
# The status codes are bitfields that can be combined with biwise or
#
# Examples :
#
# $code = StatMsg::EXON  | StatMsg::DELETE | StatMsg::CDS |
#         StatMsg::SHORT | StatMsg::FRAMESHIFT;
#

use constant DELETE                 => 0x00000001;
use constant INSERT                 => 0x00000002;
use constant CDS                    => 0x00000004;
use constant UTR                    => 0x00000008;

use constant ENTIRE                 => 0x00000010;
use constant LONG                   => 0x00000020;
use constant MEDIUM                 => 0x00000040;
use constant SHORT                  => 0x00000080;

use constant FRAMESHIFT             => 0x00000100;

use constant NO_CDS_LEFT            => 0x00000200;
use constant STRAND_FLIP            => 0x00000400;
use constant INVERT                 => 0x00000800;
use constant SCAFFOLD_SPAN          => 0x00001000;
use constant DOESNT_TRANSLATE       => 0x00002000;

use constant TRANSCRIPT             => 0x00004000;
use constant EXON                   => 0x00008000;

use constant FIVE_PRIME             => 0x00010000;
use constant THREE_PRIME            => 0x00020000;
use constant MIDDLE                 => 0x00040000;

use constant CONFUSED               => 0x00080000;
use constant ALL_INTRON             => 0x00100000;
use constant TRANSLATES             => 0x00200000;
use constant SPLIT                  => 0x00400000;
use constant NO_SEQUENCE_LEFT       => 0x00800000;
use constant PARTIAL                => 0x01000000;

use vars qw(@EXPORT_OK @ISA);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&push_err &pop_err &ec2str);


our $CUR_ID = 0;
our $LOGGER = undef;


sub new {
  my $class = shift;
  my $code  = shift;

  if(!defined($code)) {
    die("Status Code Argument is required.\n");
  }

  my $sm = bless {'code' => $code, 'id' => $CUR_ID++}, $class;

  if($LOGGER) {
    $LOGGER->add_StatMsg($sm);
  }

  return $sm;
}


#
# Sets a logger for stat msg logging
# If set, every created StatMsg will be logged to a file
# if undef, nothing will be logged.
#
sub set_logger {
  my $logger = shift;
  $LOGGER = $logger;
}


#
# Unique id is given to each statmsg so that they can be distinguished.
# Sometimes the same statmsg will be given to multiple exons and it is handy
# to be able to tell the messages apart, for example.
#
sub id {
  my $self = shift;
  return $self->{'id'};
}

#
# The code is a bitvector - a combination of bitwise or'd constants.
#
sub code {
  my $self = shift;
  return $self->{'code'};
}

#
# Returns a human readible version of what the code of this statmsg means.
#
sub code_str {
  my $self = shift;
  return code2str($self->{'code'});
}

#
# converts a code to a string
#
sub code2str {
  my $code = shift;

  my $str = "StatCode $code:";

  if(!$code) {
    $str .= '  OK';
  }
  if($code & ENTIRE) {
    $str .= ' entire';
  }
  if($code & PARTIAL) {
    $str .= ' partial';
  }
  if($code & LONG) {
    $str .= ' long';
  }
  if($code & MEDIUM) {
    $str .= ' medium';
  }
  if($code & SHORT) {
    $str .= ' short';
  }
  if($code & FRAMESHIFT) {
    $str .= ' frameshifting';
  }
  if($code & EXON) {
    $str .=  ' exon';
  }
  if($code & TRANSCRIPT) {
    $str .=  ' transcript';
  }
  if($code & DELETE) {
    $str .= ' deletion';
  }
  if($code & INSERT) {
    $str .= ' insertion';
  }
  if($code & CDS) {
    $str .= ' in CDS';
  }
  if($code & UTR) {
    $str .= ' in UTR';
  }
  if($code & FIVE_PRIME) {
    $str .= " at 5' end";
  }
  if($code & THREE_PRIME) {
    $str .= " at 3' end";
  }
  if($code & MIDDLE) {
    $str .= " in middle";
  }
  if($code & DOESNT_TRANSLATE) {
    $str .= " does not translate";
  }
  if($code & SCAFFOLD_SPAN) {
    $str .= " spans multiple scaffolds";
  }
  if($code & INVERT) {
    $str .= " inversion";
  }
  if($code & STRAND_FLIP) {
    $str .= " flips strands";
  }
  if($code & NO_CDS_LEFT) {
    $str .= " all CDS deleted";
  }
  if($code & CONFUSED) {
    $str .= " confused";
  }
  if($code & ALL_INTRON) {
    $str .= " consumed by frameshift intron";
  }
  if($code & TRANSLATES) {
    $str .= " translates";
  }
  if($code & SPLIT) {
    $str .= " split";
  }
  if($code & NO_SEQUENCE_LEFT) {
    $str .= " no sequence left";
  }


  return "$str\n";
}




1;
