package Bio::EnsEMBL::IdMapping::StableIdGenerator::EnsemblGeneric;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http:#www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub initial_stable_id {
  my $self = shift;
  my $type = shift;

  my $init_stable_id;

  # use stable ID from configuration if set
  if ($init_stable_id = $self->conf->param("starting_${type}_stable_id")) {
    $self->logger->debug("Using pre-configured $init_stable_id as base for new $type stable IDs\n");
    return $init_stable_id;
  }

  my $s_dba = $self->cache->get_DBAdaptor('source');
  my $s_dbh = $s_dba->dbc->db_handle;
  my $sql = qq(SELECT MAX(stable_id) FROM ${type}_stable_id);
  $init_stable_id = $self->fetch_value_from_db($s_dbh, $sql);

  if ($init_stable_id) {
    $self->logger->debug("Using $init_stable_id as base for new $type stable IDs\n");
  } else {
    $self->logger->warning("Can't find highest ${type}_stable_id in source db.\n");
  }

  return $init_stable_id;
}


sub increment_stable_id {
  my $self = shift;
  my $stable_id = shift;

  unless ($stable_id and ($stable_id =~ /ENS([A-Z]{1,4})(\d{11})/)) {
    throw("Unknown or missing stable ID: $stable_id.");
  }

  my $number = $2;
  my $new_stable_id = 'ENS'.$1.(++$number);

  return $new_stable_id;
}


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

