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

use Bio::EnsEMBL::Root;



@ISA = qw(Bio::EnsEMBL::Root);


# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;



# set stuff in self from @args
 return $self; # success - we hope!
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


=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: return or assign the value of end, which is a position within
           the exon given by end_exon.
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}


=head2 start_Exon

 Title   : start_exon
 Usage   : $obj->start_Exon($newval)
 Function: return or assign the value of start_exon, which denotes the
           exon at which translation starts (and within this exon, at the
           position indicated by start, see above).
 Returns : value of start_exon (Exon object)
 Args    : newvalue (optional)


=cut

sub start_Exon {
   my $self = shift;

   if( @_ ) {
      my $value = shift;
      if( !ref $value || !$value->isa('Bio::EnsEMBL::Exon') ) {
         $self->throw("Got to have an Exon object, not a $value");
      }
      $self->{'start_exon'} = $value;
    }
   return $self->{'start_exon'};
}




=head2 end_Exon

 Title   : end_exon
 Usage   : $obj->end_Exon($newval)
 Function: return or assign the value of end_exon, which denotes the
           exon at which translation ends (and within this exon, at the
           position indicated by end, see above).
 Returns : value of end_exon (Exon object)
 Args    : newvalue (optional)


=cut

sub end_Exon {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if( !ref $value || !$value->isa('Bio::EnsEMBL::Exon') ) {
         $self->throw("Got to have an Exon object, not a $value");
      }
      $self->{'end_exon'} = $value;
    } 

    return $self->{'end_exon'};
}


=head2 version

 Title   : version
 Usage   : $obj->version()
 Function: 
 Returns : value of version
 Args    : 

=cut

sub version{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  modified dates are loaded. Ignoring set value $value");
      return;
    }

    if( exists $self->{'_version'} ) {
      return $self->{'_version'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_version'};

}



=head2 stable_id

 Title   : stable_id
 Usage   : $obj->stable_id
 Function: 
 Returns : value of stable_id
 Args    : 

=cut

sub stable_id{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_stable_id'} = $value;
      return;
    }

    if( exists $self->{'_stable_id'} ) {
      return $self->{'_stable_id'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_stable_id'};

}


sub _get_stable_entry_info {
   my $self = shift;

   if( !defined $self->adaptor ) {
     return undef;
   }

   $self->adaptor->get_stable_entry_info($self);

}



=head2 temporary_id

 Title   : temporary_id
 Usage   : $obj->temporary_id($newval)
 Function: 
 Returns : value of temporary_id
 Args    : newvalue (optional)

=cut

sub temporary_id {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'tempID'} = $value;
    }
    return $self->{'tempID'};

}



=head2 dbID

 Title   : dbID
 Usage   : $obj->dbID($newval)
 Function: 
 Returns : value of dbID
 Args    : newvalue (optional)

=cut

sub dbID {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'dbID'} = $value;
    }
    return $self->{'dbID'};

}



=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: 
 Returns : value of adaptor
 Args    : newvalue (optional)

=cut

sub adaptor {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}



=head2 transform

  Arg  1    : hashref $old_new_exon_map
              a hash that maps old to new exons for a whole gene
  Function  : maps start end end exon according to mapping table
              if an exon is not mapped, just keep the old one
  Returntype: none
  Exceptions: none
  Caller    : Transcript->transform() 

=cut

sub transform {
  my $self = shift;
  my $href_exons = shift;

  my $start_exon = $self->start_Exon();
  my $end_exon = $self->end_Exon();

  if ( exists $$href_exons{$start_exon} ) {
    $self->start_exon($$href_exons{$start_exon});
  } else {
    # do nothing, the start exon wasnt mapped
  }

  if ( exists $$href_exons{$end_exon} ) {
    $self->end_exon($$href_exons{$end_exon});
  } else { 
    # do nothing, the end exon wasnt mapped
  }
}



=head2 id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use stable_id or dbID instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub id{
   my $self = shift;
   my $value = shift;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l id deprecated. Please choose from stable_id or dbID");


  if( defined $value ) {
    $self->warn("$f:$l stable ids are loaded separately and dbIDs are " .
		"generated on writing. Ignoring set value $value");
    return;
  }

   if( defined $self->stable_id ) {
     return $self->stable_id();
   } else {
     return $self->dbID;
   }
}



=head2 end_exon_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use end_exon to set the end exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub end_exon_id{
   my $self = shift;
   my $value = shift;

   if ( defined $value ) {
     $self->throw( "Please use end_exon to set the end exon in Translation" .
		   "object.Translation objects changed." );
   }

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l end_exon_id is deprecated. Please use end_exon to get " .
	       "out the object. We have returned the dbID of the Exon object. "
             . "This might not be what you want!");

   return $self->end_exon->dbID;
}



=head2 start_exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use start_Exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub start_exon {
  my ($self, @args) = @_;

  $self->warn("start_exon has been renamed start_Exon " . caller);

  $self->start_Exon(@args);
}


=head2 start_exon_id

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use start_exon to set the start exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub start_exon_id{
   my $self = shift;
   my $value = shift;

   if( defined $value ) {
     $self->throw( "Please use start_exon to set the start exon in " .
		   "Translation object.Translation objects changed." );
   }
   my ($p,$f,$l) = caller;
   $self->warn("$f:$l start_exon_id is deprecated. Please use start_exon " .
	       "to get out the object. We have returned the dbID of the " .
	       "Exon object. This might not be what you want!");

   return $self->start_exon->dbID;
}






1;
