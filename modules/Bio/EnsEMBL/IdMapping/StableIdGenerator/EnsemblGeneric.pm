=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric - default Ensembl
StableIdGenerator implementation

=head1 SYNOPSIS

  # inject the confiured StableIdGenerator plugin
  my $stable_id_generator = $conf->param('plugin_stable_id_generator');
  inject($stable_id_generator);

  # create a new StableIdGenerator object
  my $generator_instance = $stable_id_generator->new(
    -LOGGER => $self->logger,
    -CONF   => $self->conf,
    -CACHE  => $self->cache
  );

  # determine starting stable ID for new assignments
  my $new_stable_id = $generator_instance->initial_stable_id('gene');

  # loop over genes
  foreach my $target_gene (@all_target_genes) {

    # if the stable Id for this gene was mapped, assign it
    if ( $mapping_exists{ $target_gene->id } ) {
      my $source_gene = $mappings{ $target_gene->id };
      $target_gene->stable_id( $source_gene->stable_id );

      # calculate and set version
      my $version =
        $generator_instance->calculate_version( $source_gene,
        $target_gene );
      $target_gene->version($version);

      # no mapping exists, assign a new stable Id
    } else {
      $target_gene->stable_id($new_stable_id);
      $target_gene->version('1');

    # increment the stable Id (to be assigned to the next unmapped gene)
      $new_stable_id =
        $generator_instance->increment_stable_id($new_stable_id);
    }
  }

=head1 DESCRIPTION

This is the default implementation for a StableIdGenerator, which
is used by Bio::EnsEMBL::IdMapping::StableIdMapper to generate new
stable Ids and increment versions on mapped stable Ids.  Refer to the
documentation in this module if you would like to implement your own
StableIdGenerator.

The stable Id mapping application allows you to plugin your own
implementation by specifying it with the --plugin_stable_id_generator
configuration parameter.

Requirements for a StableIdGenerator plugin:

  - inherit from Bio::EnsEMBL::IdMapping::BaseObject
  - implement all methods listed in METHODS below (see method POD for
    signatures)

=head1 METHODS

  initial_stable_id
  increment_stable_id
  calculate_version

=cut

package Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 initial_stable_id

  Arg[1]      : String $type - an entity type (gene|transcript|translation|exon)
  Example     : my $new_stable_id = $generator->initial_stable_id('gene');
  Description : Determine the initial stable Id to use for new assignments. This
                method is called once at the beginning of stable Id mapping.
  Return type : String - a stable Id of appropriate type
  Exceptions  : none
  Caller      : Bio::EnsEMBL::IdMapping::StableIdMapper::map_stable_ids()
  Status      : At Risk
              : under development

=cut

sub initial_stable_id {
  my $self = shift;
  my $type = shift;

  my $init_stable_id;

  # use stable ID from configuration if set
  if ($init_stable_id = $self->conf->param("starting_${type}_stable_id")) {
    $self->logger->debug("Using pre-configured $init_stable_id as base for new $type stable IDs.\n");
    return $init_stable_id;
  }

  my $s_dba = $self->cache->get_DBAdaptor('source');
  my $s_dbh = $s_dba->dbc->db_handle;

  # look in the ${type}_stable_id table first
  my $sql = qq(SELECT MAX(stable_id) FROM ${type}_stable_id);
  $init_stable_id = $self->fetch_value_from_db($s_dbh, $sql);

  # also look in gene_archive to make sure there are no larger Ids there
  unless ($type eq 'exon') {
    $sql = qq(SELECT MAX(${type}_stable_id) FROM gene_archive);
    my $archived_stable_id = $self->fetch_value_from_db($s_dbh, $sql);
    if ($archived_stable_id and $self->is_valid($archived_stable_id) and
        ($archived_stable_id gt $init_stable_id)) {
      $init_stable_id = $archived_stable_id;
    }
  }

  if ($init_stable_id) {
    # since $init_stable_id now is the highest existing stable Id for this
    # object type, we need to increment it to find the first one we want to use
    # for new assignments
    $init_stable_id = $self->increment_stable_id($init_stable_id);
    
    $self->logger->debug("Using $init_stable_id as base for new $type stable IDs.\n");

  } else {
    $self->logger->warning("Can't find highest ${type}_stable_id in source db.\n");
  }

  return $init_stable_id;
}


=head2 increment_stable_id

  Arg[1]      : String $stable_id - the stable Id to increment
  Example     : $next_stable_id = $generator->increment_stable_id(
                  $current_stable_id);
  Description : Increments the stable Id used for new assigments. This method is
                called after each new stable Id assigment to generate the next
                stable Id to be used.
  Return type : String - the next new stable Id
  Exceptions  : thrown on missing or malformed argument
  Caller      : Bio::EnsEMBL::IdMapping::StableIdMapper::map_stable_ids()
  Status      : At Risk
              : under development

=cut

sub increment_stable_id {
  my $self = shift;
  my $stable_id = shift;

  unless ($self->is_valid($stable_id)) {
    throw("Unknown or missing stable ID: $stable_id.");
  }
  
  $stable_id =~ /ENS([A-Z]{1,4})(\d{11})/;
  my $number = $2;
  my $new_stable_id = 'ENS'.$1.(++$number);

  return $new_stable_id;
}


=head2 is_valid

  Arg[1]      : String $stable_id - the stable Id to check
  Example     : unless ($generator->is_valid($stable_id)) {
                  die "Invalid stable Id: $stable_id.\n";
                }
  Description : Tests a stable Id to be valid (according to the Ensembl stable
                Id format definition).
  Return type : Boolean - TRUE if valid, FALSE otherwise
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_valid {
  my ($self, $stable_id) = @_;
  return ($stable_id and ($stable_id =~ /ENS([A-Z]{1,4})(\d{11})/));
}


=head2 calculate_version

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyFeature $s_obj - source object
  Arg[2]      : Bio::EnsEMBL::IdMapping::TinyFeature $t_obj - target object
  Example     : my $version = $generator->calculate_version($source_gene,
                  $target_gene);
                $target_gene->version($version);
  Description : Determines the version for a mapped stable Id. For Ensembl
                genes, the rules for incrementing the version number are:
                    - exons: if exon sequence changed
                    - transcript: if spliced exon sequence changed
                    - translation: if transcript changed
                    - gene: if any of its transcript changed
  Return type : String - the version to be used
  Exceptions  : thrown on wrong argument
  Caller      : Bio::EnsEMBL::IdMapping::StableIdMapper::map_stable_ids()
  Status      : At Risk
              : under development

=cut

sub calculate_version {
  my $self = shift;
  my $s_obj = shift;
  my $t_obj = shift;

  my $version = $s_obj->version;

  if ($s_obj->isa('Bio::EnsEMBL::IdMapping::TinyExon')) {
    
    # increment version if sequence changed
    $version++ unless ($s_obj->seq eq $t_obj->seq);
  
  } elsif ($s_obj->isa('Bio::EnsEMBL::IdMapping::TinyTranscript')) {
  
    # increment version if spliced exon sequence changed
    $version++ unless ($s_obj->seq_md5_sum eq $t_obj->seq_md5_sum);

  } elsif ($s_obj->isa('Bio::EnsEMBL::IdMapping::TinyTranslation')) {

    # increment version if transcript changed
    my $s_tr = $self->cache->get_by_key('transcripts_by_id', 'source',
      $s_obj->transcript_id);
    my $t_tr = $self->cache->get_by_key('transcripts_by_id', 'target',
      $t_obj->transcript_id);

    $version++ unless ($s_tr->seq_md5_sum eq $t_tr->seq_md5_sum);
    
  } elsif ($s_obj->isa('Bio::EnsEMBL::IdMapping::TinyGene')) {
    
    # increment version if any transcript changed
    my $s_tr_ident = join(":", map { $_->stable_id.'.'.$_->version }
      sort { $a->stable_id cmp $b->stable_id }
        @{ $s_obj->get_all_Transcripts });
    my $t_tr_ident = join(":", map { $_->stable_id.'.'.$_->version }
      sort { $a->stable_id cmp $b->stable_id }
        @{ $t_obj->get_all_Transcripts });

    $version++ unless ($s_tr_ident eq $t_tr_ident);
    
  } else {
    throw("Unknown object type: ".ref($s_obj));
  }

  return $version;
}


1;

