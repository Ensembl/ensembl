
#
# BioPerl module for Bio::EnsEMBL::Translation
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Translation - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Translation;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;



@ISA = qw(Bio::Root::Object);


# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: return or assign the value of start, which is a position within
           the exon given by start_exon_id.
 Returns : value of start
 Args    : newvalue (optional)


=cut

sub start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}

=head2 start_exon_id

 Title   : start_exon_id
 Usage   : $obj->start_exon_id($newval)
 Function: return or assign the value of start_exon_id, which denotes the
           exon at which translation starts (and within this exon, at the
           position indicated by start, see above).
 Returns : value of start_exon_id
 Args    : newvalue (optional)


=cut

sub start_exon_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start_exon_id'} = $value;
    }
    return $obj->{'start_exon_id'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: return or assign the value of end, which is a position within
           the exon given by end_exon_id. 
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}

=head2 end_exon_id

 Title   : end_exon_id
 Usage   : $obj->end_exon_id($newval)
 Function: return or assign the value of end_exon_id, which denotes the
           last exon of a translation (and within this exon, at the
           position indicated by end, see above).
 Returns : value of end_exon_id
 Args    : newvalue (optional)


=cut

sub end_exon_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end_exon_id'} = $value;
    }
    return $obj->{'end_exon_id'};

}
