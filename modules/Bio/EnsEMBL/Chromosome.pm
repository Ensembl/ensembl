# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 07.04.2000
# Last modified : 16.02.2003 by James Smith
#
# Copyright EMBL-EBI/WTSI 2000
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Chromosome

=head1 SYNOPSIS

=head1 DESCRIPTION

Contains very basic information of a chromosome and access methods
for global features of a chromosome. It does not have the sequence or
more detailed information - check out SliceAdaptor for that (you will
want to make a slice of the chromosome)
    
=head1 CONTACT 

Post questions to the EnsEMBL developer mailing list: <ensembl-dev@ebi.ac.uk>

=cut

# Let the code begin...

package Bio::EnsEMBL::Chromosome;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;

@ISA = qw( Bio::EnsEMBL::Root );

## Changes to new by James Smith - no longer takes statistics
## on creation.. these get added later.

=head2 new

  Args [...] : List of named arguments 
  Example    : $chr = new Chromosome(-chr_name      => $name,
                                     -dbID          => $dbID,
                                     -adaptor       => $adaptor,
                                     -length        => $length);
  Description: Creates a new chromosome object
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : thrown if the adaptor or chr_name argument is not supplied
  Caller     : Bio::EnsEMBL::DBSQL::ChromosomeAdaptor

=cut

sub new {
  my ($class,@args) = @_;
  my $self = {};
  bless($self, $class);
   
  my ( $chr_name, $chromosome_id, $adaptor, $length ) =
    $self->_rearrange([qw(CHR_NAME DBID ADAPTOR LENGTH)], @args);

  $self->throw("Badly formed chromosome")
    unless defined $chr_name && defined $adaptor;

  $self->adaptor(  $adaptor       );
  $self->chr_name( $chr_name      );
  $self->dbID(     $chromosome_id );
  $self->length(   $length        );
  $self->{'stats'} ={};
  return $self;
}

=head2 chr_name

  Arg [1]    : string $chr_name
  Example    : none
  Description: get/set for attribute chromosome name
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub chr_name{
  my( $self,$value ) = @_;
         $self->{'chr_name'} = $value if  defined $value;
  return $self->{'chr_name'};
}

=head2 stats

  Arg [1]    : 
  Example    : $obj->stats()
  Description: returns hashref of additional statistics stored in Chromosome table
  Returntype : hashref
  Exceptions : none
  Caller     : mapview

=cut

sub stats {
   my $self = shift;
   return %{$self->{'stats'}};
}

=head2 stat

  Arg [1]    : Name of attribute to set
  Arg [2]    : Optional value to set to...
  Example    : $obj->stat( $value)
  Description: get/set for a statistic...
  Returntype : number
  Exceptions : none
  Caller     : mapview

=cut

sub stat {
  my( $self, $key, $value ) = @_;
         $self->{'stats'}{$key} = $value if defined $value;
  return $self->{'stats'}{$key};
}

=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::ChromosomeAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::ChromsomeAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub adaptor {
  my( $self,$value ) = @_;
         $self->{'adaptor'} = $value if defined $value;
  return $self->{'adaptor'};
}

=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: get/set for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub dbID {
  my ($self, $value) = @_;
         $self->{'_dbID'} = $value if defined $value;
  return $self->{'_dbID'};
}



=head2 length

  Arg [1]    : int $length
  Example    : none
  Description: get/set for the attribute length, the Chromosomes length in 
               basepairs
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub length {
  my ($self, $length) = @_;
         $self->{'length'} = $length if defined $length;
  return $self->{'length'};
}

## Deprecated calls - these should now use "stat"
sub xref_genes    { return $_[0]->stat('xref_genes'); }
sub known_genes   { return $_[0]->stat('known_genes'); }
sub unknown_genes { return $_[0]->stat('unknown_genes'); }
sub snps          { return $_[0]->stat('snps'); }

1;

