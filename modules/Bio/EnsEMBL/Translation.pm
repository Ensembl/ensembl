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


package Bio::EnsEMBL::Translation;
use vars qw($AUTOLOAD @ISA);
use strict;


use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  return bless {},$class;
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
      $self->{'_version'} = $value;    
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
      $self->{'adaptor'} = shift;
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

  if ( exists $href_exons->{$start_exon} ) {
    $self->start_Exon($href_exons->{$start_exon});
  } else {
    # do nothing, the start exon wasnt mapped
  }

  if ( exists $href_exons->{$end_exon} ) {
    $self->end_Exon($href_exons->{$end_exon});
  } else { 
    # do nothing, the end exon wasnt mapped
  }
}


=head2 get_all_DBEntries

  Arg [1]    : none
  Example    : @dbentries = @{$gene->get_all_DBEntries()};
  Description: Retrieves DBEntries (xrefs) for this translation.  

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the translation (i.e. they have not already been added or 
               loaded).
  Returntype : list reference to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, TranslationAdaptor::store

=cut

sub get_all_DBEntries {
  my $self = shift;

  #if not cached, retrieve all of the xrefs for this gene
  if(!defined $self->{'dbentries'} && $self->adaptor()) {
    $self->{'dbentries'} = 
      $self->adaptor->db->get_DBEntryAdaptor->fetch_all_by_Translation($self);
  }

  $self->{'dbentries'} ||= [];

  return $self->{'dbentries'};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : @dbentries = @{$gene->get_all_DBEntries()};
  Description: Associates a DBEntry with this gene. Note that adding DBEntries
               will prevent future lazy-loading of DBEntries for this gene
               (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general

=cut

sub add_DBEntry {
  my $self = shift;
  my $dbe = shift;

  unless($dbe && ref($dbe) && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    $self->throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}


=head2 get_all_DBLinks

  Arg [1]    : see get_all_DBEntries
  Example    : see get_all_DBEntries
  Description: This is here for consistancy with the Transcript and Gene 
               classes.  It is a synonym for the get_all_DBEntries method.
  Returntype : see get_all_DBEntries
  Exceptions : none
  Caller     : general

=cut

sub get_all_DBLinks {
  my $self = shift;

  return $self->get_all_DBEntries(@_);
}


1;
