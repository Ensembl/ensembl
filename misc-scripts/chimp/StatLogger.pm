use strict;
use warnings;

use FileHandle;
use IO::File;

use Bio::EnsEMBL::Utils::Exception qw(throw);

#
# This is a class for logging StatMsgs to a file
# all it does is write stat codes to a file and
# ensure that duplicates are not written twice.
#

package StatLogger;

sub new {
  my $class = shift;
  my $filename = shift;

  my $fh;

  if($filename) {
    $fh = IO::File->new();
    $fh->open(">$filename") or throw("Could not open file $filename.");
    $fh->autoflush();
  } else {
    $fh = \*STDOUT;
  }

  return bless {'fh' => $fh,
                'seen_msgs' => {}}, $class;
}


sub add_StatMsg {
  my $self = shift;
  my $msg  = shift;

  my $fh = $self->{'fh'};

  if(!$self->{'seen_msgs'}->{$msg->id}) {
    print $fh $msg->code(), "\n";
    $self->{'seen_msgs'}->{$msg->id()} = 1;
  }
}


1;
