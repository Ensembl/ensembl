=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::ArchiveStableId

=head1 DESCRIPTION

ArchiveStableId objects are the main workunit for retrieving stable id
archived information from EnsEMBL core database.

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
 
=head1 SEE ALSO

Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
Bio::EnsEMBL::StableIdEvent
Bio::EnsEMBL::StableIdHistoryTree

=cut

package Bio::EnsEMBL::ArchiveStableId;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

=head2 new

  Arg [STABLE_ID]      : String $stable_id 
  Arg [VERSION]        : Int $version 
  Arg [CURRENT_VERSION]: Int $current_version 
  Arg [DB_NAME]        : String $db_name 
  Arg [RELEASE]        : String $release
  Arg [ASSEMBLY_NAME]  : String $assembly
  Arg [TYPE]           : String $type - "Gene", "Transcript", "Translation", "Exon"
  Arg [ADAPTOR]        : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor $adaptor 
  Description          : standard constructor with named arguments to create
                         ArchiveStableId
  Returntype           : Bio::EnsEMBL::ArchiveStableId
  Exceptions           : none
  Caller               : general, Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor
  Status               : At Risk
                       : under development

=cut

sub new {
  my $class = shift;
  $class = ref( $class ) || $class;

  my $self = bless {}, $class;

  my ($stable_id, $version, $current_version, $db_name, $release, $assembly,
    $type, $adaptor) = rearrange([qw( STABLE_ID VERSION CURRENT_VERSION DB_NAME
    RELEASE ASSEMBLY TYPE ADAPTOR)], @_ );

  $self->{'stable_id'} = $stable_id;
  $self->{'version'} = $version;
  $self->{'current_version'} = $current_version;
  $self->{'db_name'} = $db_name;
  $self->{'release'} = $release;
  $self->{'assembly'} = $assembly;
  $self->{'type'} = $type;
  $self->adaptor($adaptor);

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
  Arg [8]     : Int $current_version 
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
      'adaptor' => $_[6],
      'current_version' => $_[7],
  }, $class;

  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );

  return $self;
}


=head2 get_history_tree

  Arg[1]      : (optional) Int $num_high_scorers
                number of mappings per stable ID allowed when filtering
  Arg[2]      : (optional) Int $max_rows
                maximum number of stable IDs in history tree (used for
                filtering)
  Example     : my $history_tree = $archive_id->get_history_tree;
  Description : Returns the history tree of this ArchiveStableId
  Return type : Bio::EnsEMBL::StableIdHistoryTree
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_history_tree {
  my ($self, $num_high_scorers, $max_rows) = @_;
  
  unless ($self->{'history'}) {
    $self->{'history'} = $self->adaptor->fetch_history_tree_by_stable_id(
      $self->stable_id, $num_high_scorers, $max_rows);
  }

  return $self->{'history'};
}


=head2 get_event

  Args        : stable_id
  Description : Retrieve a specific event for this archive and a given stable id
  Returntype  : listref of Bio::EnsEMBL::StableIdEvent
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_event {
  my ($self, $stable_id) = @_;

  my $event = $self->adaptor->fetch_stable_id_event($self, $stable_id);

  return $event;
}


=head2 get_all_predecessors

  Args        : none
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

  if ( lc( $self->type() ) eq 'translation' ) {
    return $self->adaptor->get_peptide($self);
  } else {
    return undef;
  }
}


=head2 get_all_associated_archived

  Example     : my ($arch_gene, $arch_tr, $arch_tl, $pep_seq) =
                  @{ $arch_id->get_all_associated_archived };
  Description : Fetches associated archived stable IDs from the db for this
                ArchiveStableId (version is taken into account).
  Return type : Listref of
                  ArchiveStableId archived gene
                  ArchiveStableId archived transcript
                  (optional) ArchiveStableId archived translation
                  (optional) peptide sequence
  Caller      : webcode, general
  Status      : At Risk
              : under development

=cut

sub get_all_associated_archived {
  my $self = shift;
  return $self->adaptor->fetch_associated_archived($self);
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

  Example     : my @archived_transcripts =
                  @{ $arch_id->get_all_transcript_archive_ids };
  Description : Returns transcript ArchiveStableIds associated with this
                ArchiveStableId. If this is a transcript, it returns itself.
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
                ArchiveStableId. If this is a translation, it returns itself.
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


=head2 current_version

  Example     : if (my $v = $arch_id->current_version) {
                  print "Current version of this stable ID ", $v, "\n";
                } else {
                  print "This stable ID is not in the current db.\n";
                }
  Description : Lazy-loads the current version of stable ID
  Return type : Boolean (TRUE is current version found, else FALSE)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub current_version {
  my $self = shift;

  if (@_) {
    $self->{'current_version'} = shift;
  } elsif (! defined $self->{'current_version'}) {
    if (defined $self->adaptor()) {
      # lazy load
      $self->adaptor()->lookup_current($self);
    }       
  }

  return $self->{'current_version'};
}


=head2 is_current

  Example     : if ($arch_id->is_current) {
                  print $arch_id->version, " is the current version of this
                    stable ID.\n";
                }
  Description : Determines if the version of this object is the current version
                of this stable ID. Note that this method doesn't lazy-load the
                current version of an ArchiveStableId; if you want to be sure,
                use current_version() instead.
  Return type : Boolean (TRUE if it is current, else FALSE)
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_current {
  my $self = shift;
  return ($self->{'version'} == $self->{'current_version'});
}


=head2 get_latest_incarnation

  Example     : my $latest = $arch_id->get_latest_incarnation;
                print "Latest version of ".$arch_id->stable_id." is ".
                  $latest->version."\n";
  Description : Returns the ArchiveStableId representing the latest version
                of this stable ID. Returns itself if this already is the latest
                version, otherwise fetches it from the db.
  Return type : Bio::EnsEMBL::ArchiveStableId
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_latest_incarnation {
  my $self = shift;

  return $self if ($self->is_latest);

  my $latest = $self->adaptor->fetch_by_stable_id($self->stable_id);
  return $latest;
}


=head2 is_latest

  Arg[1]      : (optional) Boolean $is_latest
  Example     : if ($arch_id->is_latest) {
                  print "Version ".$arch_id->version." is the latest version 
                    of ".$arch_id->stable_id."\n";
                }
  Description : Indicates whether this is the latest version of this stable ID.
                Can also be used as a setter if we know this is the latest
                version.
  Return type : Boolean (TRUE if yes, FALSE if no)
  Exceptions  : none
  Caller      : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor->fetch_by_stable_id, general
  Status      : At Risk
              : under development

=cut

sub is_latest {
  my $self = shift;
  $self->{'is_latest'} = shift if (@_);
  return ($self->{'is_latest'} || $self->is_current);
}


#
# getter/setters for attributes
#

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if (@_);
  return $self->{'stable_id'};
}

sub version {
  my $self = shift;
  $self->{'version'} = shift if (@_);
  return $self->{'version'};
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

sub type {
  my $self = shift;
  $self->{'type'} = shift if (@_);
  return $self->{'type'};
}

sub adaptor {
  my $self = shift;
  weaken($self->{'adaptor'} = shift) if (@_);
  return $self->{'adaptor'};
}

sub successors {
  my $self = shift;
  $self->{'successors'} = \@_;
  return $self->{'successors'};
}

1;

