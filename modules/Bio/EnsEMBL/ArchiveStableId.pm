package Bio::EnsEMBL::ArchiveStableId;

=head1 NAME

Bio::EnsEMBL::ArchiveStableId

=head1 SYNOPSIS


=head1 DESCRIPTION

ArchiveStableId objects are the main workunit for retrieving stable id archived
information from EnsEMBL core database.

Attributes:
  type: Gene, Transcript, Translation, Exon, other, undef
  stable_id: eg. ENSG00000000001
  version: e.g. 1
  db_name: eg. homo_sapiens_core_12_31
  release: e.g. 35
  assembly: e.g. NCBI35
  successors: listref of Bio::EnsEMBL::ArchiveStableIds
  adaptor: Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor

Status: At Risk. This module is in development.
 
=head1 METHODS

    new
    new_fast
    get_all_predecessors
    get_all_successors
    get_peptide
    get_all_transcript_archive_ids
    get_all_translation_archive_ids

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Ensembl core API team
Currently maintained by Patrick Meidl <meidl@ebi.ac.uk>

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut
  

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Root;
our @ISA = qw(Bio::EnsEMBL::Root);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate);


=head2 new

  Arg [STABLE_ID] : String $stable_id 
  Arg [VERSION]   : Int $version 
  Arg [DB_NAME]   : String $db_name 
  Arg [RELEASE]   : String $release
  Arg [ASSEMBLY_NAME] : String $assembly
  Arg [TYPE]      : String $type - "Gene", "Transcript", "Translation", "Exon"
  Arg [ADAPTOR]   : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor $adaptor 
  Example         : none
  Description     : standard constructor with named arguments to create
                    ArchiveStableId
  Returntype      : Bio::EnsEMBL::ArchiveStableId
  Exceptions      : none
  Caller          : general, Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
  Status          : At Risk
                  : under development

=cut

sub new {
  my $class = shift;
  $class = ref( $class ) || $class;

  my $self = bless {}, $class;

  my ($stable_id, $version, $db_name, $release, $assembly, $type, $adaptor) =
    rearrange([qw( STABLE_ID VERSION DB_NAME RELEASE ASSEMBLY TYPE ADAPTOR)],
    @_ );

  $self->{'stable_id'} = $stable_id;
  $self->{'version'} = $version;
  $self->{'db_name'} = $db_name;
  $self->{'release'} = $release;
  $self->{'assembly'} = $assembly;
  $self->{'type'} = $type;
  $self->{'adaptor'} = $adaptor;

  return $self;
}


=head2 new_fast

  Arg [1]     : String $stable_id 
  Arg [2]     : Int $version 
  Arg [3]     : String $db_name 
  Arg [4]     : String $release
  Arg [5]     : String $assembly
  Arg [6]     : String $type - "Gene", "Transcript", "Translation", "Exon"
  Arg [7]     : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor $adaptor 
  Example     : none
  Description : faster version of above constructor
  Returntype  : Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general, Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
  Status      : At Risk
              : under development

=cut

sub new_fast {
  my $class = shift;
  
  $class = ref ($class) || $class;

  my $self = bless {
    'stable_id' => $_[0],
      'version' => $_[1],
      'db_name' => $_[2],
      'release' => $_[3],
      'assembly' => $_[4],
      'type' => $_[5],
      'adaptor' => $_[6]
  }, $class;

  return $self;
}


=head2 get_all_predecessors

  Args        : none
  Example     : none
  Description : Retrieve a list of ArchiveStableIds that were mapped to this
                one.
  Returntype  : listref of Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_predecessors {
  my $self = shift;
  
  my $predecessors = $self->adaptor->fetch_predecessors_by_archive_id($self);
  
  foreach my $pre (@$predecessors) {
    $pre->successors($self);
  }

  return $predecessors;
}


=head2 get_all_successors

  Args        : none
  Example     : none
  Description : Retrieve a list of ArchiveStableIds that this one was mapped to.
  Returntype  : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_successors {
  my $self = shift;

  if ($self->{'successors'}) {
    return $self->{'successors'};
  } else {
    my $successors = $self->adaptor->fetch_successors_by_archive_id($self);
    return $self->successors(@$successors);
  }
}


=head2 get_peptide

  Example     : none
  Description : Retrieves the peptide string for this ArchiveStableId.
  Returntype  : String, or undef if this is not a Translation or cant be found
                in the database.
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_peptide {
  my $self = shift;

  if( $self->type() eq "Translation" ) {
    return $self->adaptor->get_peptide( $self );
  } else { 
    return undef;
  }
}


=head2 get_all_gene_archive_ids

  Example     : my @archived_genes = @{ $arch_id->get_all_gene_archive_ids };
  Description : Returns gene ArchiveStableIds associated with this
                ArchiveStableId. If this is a gene, it returns itself.
  Returntype  : listref of Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_gene_archive_ids {
  my $self = shift;

  if ($self->type eq "Gene") {
    return [$self];
  } else {
    return $self->adaptor->fetch_all_by_archive_id($self, 'Gene');
  }
}


=head2 get_all_transcript_archive_ids

  Example     : none
  Example     : my @archived_transcripts =
                  @{ $arch_id->get_all_transcript_archive_ids };
  Description : Returns transcript ArchiveStableIds associated with this
                ArchiveStableId. If this is a gene, it returns itself.
  Returntype  : listref of Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_transcript_archive_ids {
  my $self = shift;

  if ($self->type eq "Transcript") {
    return [$self];
  } else {
    return $self->adaptor->fetch_all_by_archive_id($self, 'Transcript');
  }
}


=head2 get_all_translation_archive_ids

  Example     : my @archived_peptides =
                  @{ $arch_id->get_all_translation_archive_ids };
  Description : Returns translation ArchiveStableIds associated with this
                ArchiveStableId. If this is a gene, it returns itself.
  Returntype  : listref of Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_translation_archive_ids {
  my $self = shift;

  if ($self->type eq "Translation") {
    return [$self];
  } else {
    return $self->adaptor->fetch_all_by_archive_id($self, 'Translation');
  }
}


# getter/setters for attributes

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if (@_);
  return $self->{'stable_id'};
}

sub db_name {
  my $self = shift;
  $self->{'db_name'} = shift if (@_);
  return $self->{'db_name'};
}

sub release {
  my $self = shift;
  $self->{'release'} = shift if (@_);
  return $self->{'release'};
}

sub assembly {
  my $self = shift;
  $self->{'assembly'} = shift if (@_);
  return $self->{'assembly'};
}

sub adaptor {
  my $self = shift;
  $self->{'adaptor'} = shift if (@_);
  return $self->{'adaptor'};
}

sub type {
  my $self = shift;
  $self->{'type'} = shift if (@_);
  return $self->{'type'};
}

sub successors {
  my $self = shift;
  $self->{'successors'} = \@_;
  return $self->{'successors'};
}


# lazy loading 

sub version {
  my $self = shift;
  if( @_ ) {
    $self->{'version'} = shift;

  } else {
    if( ! defined $self->{'version'} ) {
      if( defined $self->{'db_name'} && defined $self->{'adaptor'} ) {
	# lazy loading
	$self->{'adaptor'}->_lookup_version( $self );
      }
    }       
  }
  return $self->{'version'};
}


# deprecated methods (changed to more descriptive names)

sub get_translation_archive_id {
    my $self = shift;
    
    deprecate("Use get_all_translation_archive_ids() instead");
    
    return $self->get_all_translation_archive_ids;
}


1;
