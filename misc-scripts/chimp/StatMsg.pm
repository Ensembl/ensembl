use strict;
use warnings;



package StatMsg;
use Exporter;

use vars qw(@EXPORT_OK @ISA);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&push_err &pop_err &ec2str);

sub new {
  my $class = shift;
  my $code  = shift;

  if(!defined($code)) {
    die("Status Code Argument is required.\n");
  }

  return bless {'code' => $code}, $class;
}

sub code {
  return $self->{'code'};
}

sub code_str {
  return code2str($self->{'code'});
}


#converts a code to a string
sub code2str {
  my $code = shift;

  my $str = '';

  my print "StatCode $code:";

  if(!$ec) {
    $str .= '  OK';
  }
  if($ec & ENTIRE) {
    $str .= ' entire';
  }
  if($ec & LONG) {
    $str .= ' long';
  }
  if($ec & MEDIUM) {
    $str .= ' medium';
  }
  if($ec & SHORT) {
    $str .= ' short';
  }
  if($ec & FRAMESHIFT) {
    $str .= ' frameshifting';
  }
  if($ec & EXON) {
    $str .=  ' exon';
  }
  if($ec & TRANSCRIPT) {
    $str .=  ' transcript';
  }
  if($ec & DELETE) {
    $str .= ' deletion';
  }
  if($ec & INSERT) {
    $str .= ' insertion';
  }


  


  return "$str\n";
}

#
# Status codes.
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


1;
