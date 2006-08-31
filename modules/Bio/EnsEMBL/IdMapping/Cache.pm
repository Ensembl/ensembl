package Bio::EnsEMBL::IdMapping::Cache;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes);
use Bio::EnsEMBL::IdMapping::TinyFeature;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Storable qw(nfreeze thaw nstore retrieve);


# define available cache names here
my @cache_names = qw(
    exons_by_id
);


=head2 new

  Arg[1]      : 
  Example     : 
  Description : constructor
  Return type : Bio::EnsEMBL::IdMapping::Cache object
  Exceptions  : 
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logger, $conf) = rearrange(['LOGGER', 'CONF'], @_);

  unless ($logger->isa('Bio::EnsEMBL::Utils::Logger')) {
    throw("You must provide a Bio::EnsEMBL::Utils::Logger for logging.");
  }
  
  unless ($conf->isa('Bio::EnsEMBL::Utils::ConfParser')) {
    throw("You must provide configuration as a Bio::EnsEMBL::Utils::ConfParser object.");
  }
  
  my $self = {};
  bless ($self, $class);

  # initialise
  $self->logger($logger);
  $self->conf($conf);
  
  return $self;
}


sub build_cache {
  my $self = shift;
  my $dbtype = shift;
  my $slice_name = shift;
  
  my $dba = $self->get_DBAdaptor($dbtype);
  my $sa = $dba->get_SliceAdaptor;
  my $slice = $sa->fetch_by_name($slice_name);
  my $genes = $slice->get_all_Genes(undef, undef, 1);

  # biotype filter
  $genes = $self->filter_biotype($genes) if ($self->conf->param('biotypes'));
  my $i = scalar(@$genes);

  # find common coord_system
  $self->find_common_coord_systems;

  # find out whether native and common coord_system are identical
  my $need_project = 1;
  my $csid = join(':', $slice->coord_system_name, $slice->coord_system_version);
  $need_project = 0 if ($self->is_common_cs($csid));
  
  # build cache
  my $type = "$dbtype.$slice_name";
  $self->build_cache_from_genes($type, $genes, $need_project);
  undef $genes;

  # write cache to file, then flush cache to reclaim memory
  my $size = $self->write_all_to_file($type);

  return $i, $size;
}


sub filter_biotypes {
  my $self = shift;
  my $genes = shift;

  my $filtered = [];

  foreach my $biotype ($self->conf->param('biotypes')) {
    push @$filtered, grep { $_->biotype eq $biotype } @$genes;
  }

  return $filtered;
}


sub build_cache_from_genes {
  my $self = shift;
  my $type = shift;
  my $genes = shift;
  my $need_project = shift;
  
  throw("You must provide a type.") unless $type;
  throw("You must provide a listref of genes.") unless (ref($genes) eq 'ARRAY');

  #my $i = 0;
  #my $num_genes = scalar(@$genes);
  #$self->logger->init_progressbar('index_genes', $num_genes);

  # loop over genes sorted by gene location.
  # the sort will hopefully improve assembly mapper cache performance and
  # therefore speed up exon sequence retrieval
  foreach my $gene (sort { $a->start <=> $b->start } @$genes) {
    #$self->logger->log_progressbar('index_genes', ++$i, 2);
    #$self->logger->log_progress($num_genes, ++$i, 20, 3, 1);

    # create lightweigt gene
    my $lgene = Bio::EnsEMBL::IdMapping::TinyGene->new_fast([
        $gene->dbID,
        $gene->stable_id,
        $gene->start,
        $gene->end,
        $gene->strand,
    ]);

    # build gene caches
    #$self->add($type, 'genes_by_id', $gene->dbID, $lgene);
    #$self->add($type, 'genes_by_stable_id', $gene->stable_id, $lgene);
    
    # transcripts
    foreach my $tr (@{ $gene->get_all_Transcripts }) {
      my $ltr = Bio::EnsEMBL::IdMapping::TinyTranscript->new_fast([
          $tr->dbID,
          $tr->stable_id,
          $tr->start,
          $tr->end,
          $tr->strand,
      ]);

      #$lgene->add_Transcript($ltr);

      # build transcript caches
      #$self->add($type, 'transcripts_by_id', $tr->dbID, $ltr);
      #$self->add($type, 'transcripts_by_stable_id', $tr->stable_id, $ltr);
      #$self->add($type, 'genes_by_transcript_id', $tr->dbID, $lgene);

      # translation (if there is one)
      if (my $tl = $tr->translation) {
        my $ltl = Bio::EnsEMBL::IdMapping::TinyTranslation->new_fast([
            $tl->dbID,
            $tl->stable_id,
        ]);

        #$ltr->add_Translation($ltl);

        #$self->add($type, 'translations_by_id', $tl->dbID, $ltl);
        #$self->add($type, 'translations_by_stable_id', $tl->stable_id, $ltl);
        #$self->add($type, 'translations_by_transcript_id', $tr->dbID, $ltl);

        undef $tl;
      }

      # exons
      foreach my $exon (@{ $tr->get_all_Exons }) {
        my $lexon = Bio::EnsEMBL::IdMapping::TinyExon->new_fast([
            $exon->dbID,
            $exon->stable_id,
            $exon->start,
            $exon->end,
            $exon->strand,
            $exon->slice->seq_region_name,
            $exon->slice->coord_system_name,
            $exon->slice->coord_system->version,
            #$exon->slice->subseq($exon->start, $exon->end, $exon->strand),
            $need_project,
        ]);

        # get coordinates in common coordinate system if needed
        if ($need_project) {
          my @seg = @{ $exon->project($self->highest_common_cs,
                                      $self->highest_common_cs_version) };

          if (scalar(@seg) == 1) {
            my $sl = $seg[0]->to_Slice;
            $lexon->common_start($sl->start);
            $lexon->common_end($sl->end);
            $lexon->common_strand($sl->strand);
            $lexon->common_sr_name($sl->seq_region_name);
          }
        }
        
        $ltr->add_Exon($lexon);

        $self->add('exons_by_id', $type, $exon->dbID, $lexon);
        #$self->add($type, 'genes_by_exon_id', $exon->dbID, $lgene);
        #$self->add_list($type, 'transcripts_by_exon_id', $exon->dbID, $ltr);

        undef $exon;
      }

      undef $tr;
    }

    undef $gene;
  }

  #use Data::Dumper;
  #warn Data::Dumper::Dumper($genes);

}


sub add {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;
  my $val = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  $self->{'cache'}->{$name}->{$type}->{$key} = $val;

  return $self->{'cache'}->{$name}->{$type}->{$key};
}

sub add_list {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;
  my @vals = @_;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  push @{ $self->{'cache'}->{$name}->{$type}->{$key} }, @vals;

  return $self->{'cache'}->{$name}->{$type}->{$key};
}

sub get_by_key {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  return $self->{'cache'}->{$name}->{$type}->{$key};
}

sub get_by_name {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  
  return $self->{'cache'}->{$name}->{$type} || {};
}


sub merge {
  my $self = shift;
  my $name = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;

  foreach my $type (keys %{ $self->{'cache'}->{$name} || {} }) {
    (my $merged_type = $type) =~ s/^(\w+)\..+/$1/;
    
    foreach my $key (keys %{ $self->{'cache'}->{$name}->{$type} || {} }) {
      if (defined $self->{'cache'}->{$name}->{$merged_type}->{$key}) {
        warning("Duplicate key in cache: $name|$merged_type|$key. Skipping.\n");
      } else {
        $self->{'cache'}->{$name}->{$merged_type}->{$key} =
          $self->{'cache'}->{$name}->{$type}->{$key};
      }

      delete $self->{'cache'}->{$name}->{$type}->{$key};
    }
    
    delete $self->{'cache'}->{$name}->{$type};

  }
}


sub find_common_coord_systems {
  my $self = shift;

  # get adaptors for source db
  my $s_dba = $self->dba('source');
  my $s_csa = $s_dba->get_CoordSystemAdaptor;
  my $s_sa = $s_dba->get_SliceAdaptor;

  # get adaptors for target db
  my $t_dba = $self->dba('target');
  my $t_csa = $t_dba->get_CoordSystemAdaptor;
  my $t_sa = $t_dba->get_SliceAdaptor;

  # find common coord_systems
  my @s_coord_systems = $s_csa->fetch_all;
  my @t_coord_systems = $t_csa->fetch_all;
  my $found_highest = 0;

  SOURCE:
  foreach my $s_cs (@s_coord_systems) {
    
    TARGET:
    foreach my $t_cs (@t_coord_systems) {
      if ($s_cs->name eq $t_cs->name) {

        # test for identical coord_system version
        if ($s_cs->version and ($s_cs->version ne $t_cs->version)) {
          next TARGET;
        }

        # test for at least 50% identical seq_regions
        if ($self->seq_regions_compatible($s_cs, $s_sa, $t_sa)) {
          $self->add_common_cs($s_cs);
          
          unless ($found_higest) {
            $self->highest_common_cs = $s_cs->name;
            $self->highest_common_cs_version = $s_cs->version;
          }

          $found_highest = 1;

          next SOURCE;
        }
      }
    }
  }
  
}


sub seq_regions_compatible {
  my $self = shift;
  my $cs = shift;
  my $s_sa = shift;
  my $t_sa = shift;

  unless ($cs and $cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('You must provide a CoordSystem');
  }

  unless ($s_sa and $t_sa and $s_sa->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor')
          and $t_sa->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor')) {
    throw('You must provide a source and target SliceAdaptor');
  }

  my %sr_match;
  my $equal = 0;

  my $s_seq_regions = $s_sa->fetch_all($cs->name, $cs->version);
  my $t_seq_regions = $t_sa->fetch_all($cs->name, $cs->version);
  
  # sanity check to prevent divison by zero
  my $s_count = scalar(@$s_seq_regions);
  my $t_count = scalar(@$t_seq_regions)
  return(0) if ($s_count == 0 or $t_count == 0);
  
  foreach my $s_sr (@$s_seq_regions) {
    $sr_match{$s_sr->seq_region_name} = $s_sr->length;
  }

  foreach my $t_sr (@$t_seq_regions) {
    if (exists($sr_match{$t_sr->seq_region_name})) {
      $equal++;

      # return false if we have a region with same name but different length
      return(0) unless ($sr_match{$t_sr->seq_region_name} == $t_sr->length);
    }
  }

  if ($equal/$s_count > 0.5 and $equal/$t_count > 0.5) {
    return(1);
  } else {
    $self->logger->log("Only $equal seq_regions identical for ".$cs->name." ".$cs->version."\n");
    return(0);
  }
  
}


sub get_DBAdaptor {
  my $self = shift;
  my $prefix = shift;

  unless ($self->{'_dba'}->{$prefix}) {
    # connect to database
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
            -host   => $self->conf->param("${prefix}host"),
            -port   => $self->conf->param("${prefix}port"),
            -user   => $self->conf->param("${prefix}user"),
            -pass   => $self->conf->param("${prefix}pass"),
            -dbname => $self->conf->param("${prefix}dbname"),
            -group  => $prefix,
    );
    
    # explicitely set the dnadb to itself - by default the Registry assumes
    # a group 'core' for this now
    $dba->dnadb($dba);

    $self->{'_dba'}->{$prefix} = $dba;
  }

  return $self->{'_dba'}->{$prefix};
}


sub all_cache_files_exist {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  my $i = 1;
  foreach my $name (@{ $self->cache_names }) {
    $i *= $self->cache_file_exists($name, $type);
  }

  return $i;
}


sub cache_names {
  return \@cache_names;
}


sub cache_file_exists {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;

  my $cache_file = $self->cache_file($name, $type);

  if (-s $cache_file) {
    $self->logger->log("Cache file found. Will read from $cache_file.\n", 3);
    return 1;
  } else {
    $self->logger->log("No cache file found. Will build cache from db.\n", 3);
    return 0;
  }
}


sub cache_file {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;

  my $cache_file = ($self->conf->param('dumppath') || '.').
    "/$name.$type.object_cache.ser";

  return $cache_file;
}


sub instance_file {
  my $self = shift;

  my $instance_file = ($self->conf->param('dumppath') || '.').
    "/cache_instance.ser";

  return $instance_file;
}


sub write_all_to_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  foreach my $name (@{ $self->cache_names }) {
    $self->write_to_file($name, $type);
  }

  $self->write_instance_to_file;
}


sub write_to_file {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;

  # create dump directory if it doesn't exist
  if (my $dump_path = $self->conf->param('dumppath')) {
    unless (-d $dump_path) {
      system("mkdir -p $dump_path") == 0 or
        throw("Unable to create directory $dump_path.\n");
    }
  }
  
  my $cache_file = $self->cache_file($name, $type);

  eval { nstore($self->{'cache'}->{$name}->{$type}, $cache_file) };
  if ($@) {
    throw("Unable to store $cache_file: $@\n");
  }

  my $size = -s $cache_file;
  return parse_bytes($size);
}


sub write_instance_to_file {
  my $self = shift;

  # create dump directory if it doesn't exist
  if (my $dump_path = $self->conf->param('dumppath')) {
    unless (-d $dump_path) {
      system("mkdir -p $dump_path") == 0 or
        throw("Unable to create directory $dump_path.\n");
    }
  }
  
  my $instance_file = $self->instance_file;

  eval { nstore($self->{'instance'}, $instance_file) };
  if ($@) {
    throw("Unable to store $instance_file: $@\n");
  }

  my $size = -s $instance_file;
  return parse_bytes($size);
}


sub read_from_file {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;

  my $cache_file = $self->cache_file($name, $type);

  unless (-s $cache_file) {
    throw("No valid cache file found at $cache_file.");
  }

  $self->logger->log_stamped("Reading cache from file...\n");
  $self->logger->log("Cache file $cache_file.\n", 1);
  eval { $self->{'cache'}->{$name}->{$type} = retrieve($cache_file); };
  if ($@) {
    throw("Unable to retrieve cache: $@");
  }
  $self->logger->log_stamped("Done.\n");

  return $self->{'cache'}->{$name}->{$type};
}


sub read_instance_from_file {
  my $self = shift;

  my $instance_file = $self->instance_file;

  unless (-s $instance_file) {
    throw("No valid cache instance file found at $instance_file.");
  }

  eval { $self->{'instance'} = retrieve($instance_file); };
  if ($@) {
    throw("Unable to retrieve cache instance: $@");
  }

  return $self->{'instance'};
}


sub slice_names {
  my $self = shift;
  my $dbtype = shift;

  throw("You must provide a db type (source|target).") unless $dbtype;

  my $dba = $self->get_DBAdaptor($dbtype);
  my $sa = $dba->get_SliceAdaptor;

  my @slice_names = ();

  if ($self->conf->param('chromosomes')) {
    # filter by chromosome
    foreach my $chr ($self->conf->param('chromosomes')) {
      my $slice = $sa->fetch_by_region('chromosome', $chr);
      push @slice_names, $slice->name;
    }
    
  } elsif ($self->conf->param('region')) {
    # filter by region (specific slice)
    my $slice = $sa->fetch_by_name($self->conf->param('region'));
    push @slice_names, $slice->name;

  } else {
    # fetch all genes, but do in junks to save memory
    my $ga = $dba->get_GeneAdaptor;
    foreach my $srid (@{ $ga->list_seq_region_ids }) {
      my $slice = $sa->fetch_by_seq_region_id($srid);
      push @slice_names, $slice->name;
    }
  }

  return \@slice_names;
}

#
# getters/setters
#

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub logger {
  my $self = shift;
  $self->{'_logger'} = shift if (@_);
  return $self->{'_logger'};
}


=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub conf {
  my $self = shift;
  $self->{'_conf'} = shift if (@_);
  return $self->{'_conf'};
}


sub highest_common_cs {
  my $self = shift;
  $self->{'instance'}->{'hccs'} = shift if (@_);
  return $self->{'instance'}->{'hccs'};
}


sub highest_common_cs_version {
  my $self = shift;
  $self->{'instance'}->{'hccsv'} = shift if (@_);
  return $self->{'instance'}->{'hccsv'};
}


sub add_common_cs {
  my $self = shift;
  my $cs = shift;

  unless ($cs and $cs->isa('Bio::EnsEMBL::CoordSystem')) {
    throw('You must provide a CoordSystem');
  }

  my $csid = join(':', $slice->coord_system_name, $slice->coord_system_version);

  $self->{'instance'}->{'ccs'}->{$csid} = 1;
}


sub is_common_cs {
  my $self = shift;
  my $csid = shift;

  return $self->{'instance'}->{'ccs'}->{$csid};
}


1;

