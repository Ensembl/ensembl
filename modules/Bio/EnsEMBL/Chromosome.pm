# BioPerl module for Bio::EnsEMBL::Chromosome
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 07.04.2000
# Last modified : 09.04.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
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

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _
    
=cut


# Let the code begin...

package Bio::EnsEMBL::Chromosome;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Root;

@ISA = qw( Bio::EnsEMBL::Root );


=head2 new

  Args [...] : List of named arguments 
  Example    : $chr = new Chromosome(-chr_name      => $name,
                                     -dbID          => $dbID,
                                     -adaptor       => $adaptor,
                                     -length        => $length,
                                     -known_genes   => $known_genes,
                                     -xref_genes    => $xref_genes,
                                     -unknown_genes => $unknown_genes,
                                     -snps          => $snps);
  Description: Creates a new chromosome object
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : thrown if the adaptor or chr_name argument is not supplied
  Caller     : Bio::EnsEMBL::DBSQL::ChromosomeAdaptor

=cut

sub new {
    my ($class,@args) = @_;
    
    my $self = {};

    bless($self, $class);
	 
    my ( $chr_name, $chromosome_id, $adaptor, $length, 
	 $known_genes, $xref_genes, $unknown_genes, $snps) = 
	 $self->_rearrange([qw(CHR_NAME
			       DBID
			       ADAPTOR 
			       LENGTH 
			       KNOWN_GENES
			       XREF_GENES
			       UNKNOWN_GENES
			       SNPS)], 
			   @args);

    if( !defined $chr_name || !defined $adaptor ) {
      $self->throw("Badly formed chromosome");
    }

    $self->adaptor($adaptor);
    $self->chr_name($chr_name);
    $self->dbID($chromosome_id);
    $self->unknown_genes($unknown_genes);
    $self->length($length);
    $self->xref_genes($xref_genes);
    $self->known_genes($known_genes);
    $self->snps($snps);

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
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'chr_name'} = $value;
    }
    return $obj->{'chr_name'};

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
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};
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

  if(defined $value) {
    $self->{'_dbID'} = $value;
  }

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

  if(defined $length) {
    $self->{'length'} = $length;
  }

  return $self->{'length'};
}



=head2 xref_genes

  Arg [1]    : int $number_of_xref_genes
  Example    : none
  Description: get/set for the attribute xref_genes, the number of xref genes
               on this chromosome
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub xref_genes {
  my ($self, $xref_genes) = @_;

  if(defined $xref_genes) {
    $self->{'xref_genes'} = $xref_genes;
  }

  return $self->{'xref_genes'};
}

=head2 known_genes

  Arg [1]    : int $number_of_known_genes
  Example    : none
  Description: get/set for the attribute known_genes, the number of known genes
               on this chromosome
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub known_genes {
  my ($self, $known_genes) = @_;

  if(defined $known_genes) {
    $self->{'known_genes'} = $known_genes;
  }

  return $self->{'known_genes'};
}



=head2 unknown_genes

  Arg [1]    : int $number_of_unknown_genes
  Example    : none
  Description: get/set for the attribute unknown_genes, the number of unknown 
               genes on this chromosome
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub unknown_genes {
  my ($self, $unknown_genes) = @_;
  
  if(defined $unknown_genes) {
    $self->{'unknown_genes'} = $unknown_genes;
  }

  return $self->{'unknown_genes'};
}

			    

=head2 snps

  Arg [1]    : int $number_of_snps
  Example    : none
  Description: get/set for the attribute snps. The SNP count 
               on this chromosome
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub snps {
  my($self, $snps) = @_;

  if(defined $snps) {
    $self->{'snps'} = $snps;
  }

  return $self->{'snps'}
}



=head2 chromosome_id

  Args       : none
  Example    : none
  Description: deprecated, use dbID or chr_name
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub chromosome_id {
  my ($self, $id ) = @_;

  my ($package, $filename, $line) = caller();

  $self->warn("Chromosome::chromosome_id is deprecated, use Chromosome::dbID 
              instead\n line:$line package:$package filename:$filename");

  return $self->dbID($id);
}



=head2 get_landmark_MarkerFeatures

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED 
               use Bio::EnsEMBL::Slice::get_all_landmark_MarkerFeatures instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_landmark_MarkerFeatures{
   my ($self,@args) = @_;

   $self->warn("Chromosome::get_landmark_MarkerFeatures is deprecated. \n" .
	       "Use Slice::get_landmark_MarkerFeatures instead\n");

   return $self->adaptor->get_landmark_MarkerFeatures($self->chr_name);
}


1;




