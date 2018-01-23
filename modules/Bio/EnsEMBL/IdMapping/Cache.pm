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

Bio::EnsEMBL::IdMapping::Cache - a cache to hold data objects used by the 
IdMapping application

=head1 DESCRIPTION

=head1 METHODS

=cut


package Bio::EnsEMBL::IdMapping::Cache;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes path_append);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::IdMapping::TinyGene;
use Bio::EnsEMBL::IdMapping::TinyTranscript;
use Bio::EnsEMBL::IdMapping::TinyTranslation;
use Bio::EnsEMBL::IdMapping::TinyExon;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Storable qw(nstore retrieve);
use Digest::MD5 qw(md5_hex);

# define available cache names here
my @cache_names = qw(
    exons_by_id
    transcripts_by_id
    transcripts_by_exon_id
    translations_by_id
    genes_by_id
    genes_by_transcript_id
);


=head2 new

  Arg [LOGGER]: Bio::EnsEMBL::Utils::Logger $logger - a logger object
  Arg [CONF]  : Bio::EnsEMBL::Utils::ConfParser $conf - a configuration object
  Example     : my $cache = Bio::EnsEMBL::IdMapping::Cache->new(
                  -LOGGER => $logger,
                  -CONF   => $conf,
                );
  Description : constructor
  Return type : Bio::EnsEMBL::IdMapping::Cache object
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logger, $conf, $load_instance) =
    rearrange(['LOGGER', 'CONF', 'LOAD_INSTANCE'], @_);

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

  if ($load_instance) {
    $self->read_instance_from_file;
  }
  
  return $self;
}


=head2 build_cache_by_slice

  Arg[1]      : String $dbtype - db type (source|target)
  Arg[2]      : String $slice_name - the name of a slice (format as returned by
                Bio::EnsEMBL::Slice->name)
  Example     : my ($num_genes, $filesize) = $cache->build_cache_by_slice(
                  'source', 'chromosome:NCBI36:X:1:1000000:-1');
  Description : Builds a cache of genes, transcripts, translations and exons
                needed by the IdMapping application and serialises the resulting
                cache object to a file, one slice at a time.
  Return type : list of the number of genes processed and the size of the
                serialised cache file
  Exceptions  : thrown on invalid slice name
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub build_cache_by_slice {
  my $self       = shift;
  my $dbtype     = shift;
  my $slice_name = shift;

  # set cache method (required for loading cache later)
  $self->cache_method('BY_SEQ_REGION');

  my $dba = $self->get_DBAdaptor($dbtype);
  my $sa  = $dba->get_SliceAdaptor;

  my $slice = $sa->fetch_by_name($slice_name);
  unless ($slice) {
    throw("Could not retrieve slice $slice_name.");
  }

  my $genes = $slice->get_all_Genes( undef, undef, 1 );

  # find common coord_system
  my $common_cs_found = $self->find_common_coord_systems;

  # find out whether native coord_system is a common coord_system.
  # if so, you don't need to project.
  # also don't project if no common coord_system present
  my $need_project = 1;

  my $csid = join( ':',
                   $slice->coord_system_name,
                   $slice->coord_system->version );

  if ( $self->is_common_cs($csid) or !$self->highest_common_cs ) {
    $need_project = 0;
  }

  # build cache
  my $type = "$dbtype.$slice_name";
  my $num_genes =
    $self->build_cache_from_genes( $type, $genes, $need_project );
  undef $genes;

  # write cache to file, then flush cache to reclaim memory
  my $size = $self->write_all_to_file($type);

  return $num_genes, $size;
} ## end sub build_cache_by_slice


=head2 build_cache_all

  Arg[1]      : String $dbtype - db type (source|target)
  Example     : my ($num_genes, $filesize) = $cache->build_cache_all('source');
  Description : Builds a cache of genes, transcripts, translations and exons
                needed by the IdMapping application and serialises the
                resulting cache object to a file. All genes across the genome
                are processed in one go. This method should be used when
                build_cache_by_seq_region can't be used due to a large number
                of toplevel seq_regions (e.g. 2x genomes).
  Return type : list of the number of genes processed and the size of the
                serialised cache file
  Exceptions  : thrown on invalid slice name
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub build_cache_all {
  my $self = shift;
  my $dbtype = shift;
  
  # set cache method (required for loading cache later)
  $self->cache_method('ALL');

  my $dba = $self->get_DBAdaptor($dbtype);
  my $ga = $dba->get_GeneAdaptor;
  
  my $genes = $ga->fetch_all;

  # find common coord_system
  my $common_cs_found = $self->find_common_coord_systems;

  # Build cache. Setting $need_project to 'CHECK' will cause
  # build_cache_from_genes() to check the coordinate system for each
  # gene.
  my $type         = "$dbtype.ALL";
  my $need_project = 'CHECK';
  my $num_genes =
    $self->build_cache_from_genes( $type, $genes, $need_project );

  undef $genes;

  # write cache to file, then flush cache to reclaim memory
  my $size = $self->write_all_to_file($type);

  return $num_genes, $size;
}


=head2 build_cache_from_genes 

  Arg[1]      : String $type - cache type
  Arg[2]      : Listref of Bio::EnsEMBL::Genes $genes - genes to build cache
                from
  Arg[3]      : Boolean $need_project - indicate if we need to project exons to
                common coordinate system
  Example     : $cache->build_cache_from_genes(
                  'source.chromosome:NCBI36:X:1:100000:1', \@genes);
  Description : Builds the cache by fetching transcripts, translations and exons
                for a list of genes from the database, and creating lightweight
                Bio::EnsEMBL::IdMapping::TinyFeature objects containing only the
                data needed by the IdMapping application. These objects are
                attached to a name cache in this cache object. Exons only need
                to be projected to a commond coordinate system if their native
                coordinate system isn't common to source and target assembly
                itself.
  Return type : int - number of genes after filtering
  Exceptions  : thrown on wrong or missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub build_cache_from_genes {
  my $self         = shift;
  my $type         = shift;
  my $genes        = shift;
  my $need_project = shift;

  throw("You must provide a type.") unless $type;
  throw("You must provide a listref of genes.")
    unless ( ref($genes) eq 'ARRAY' );

  # biotype filter
  if ( $self->conf()->param('biotypes') ||
       $self->conf()->param('biotypes_include') ||
       $self->conf()->param('biotypes_exclude') )
  {
    $genes = $self->filter_biotypes($genes);
  }
  my $num_genes = scalar(@$genes);

  # initialise cache for the given type.
  $self->{'cache'}->{$type} = {};

  #my $i = 0;
  #my $num_genes = scalar(@$genes);
  #my $progress_id = $self->logger->init_progress($num_genes);

 # loop over genes sorted by gene location.
 # the sort will hopefully improve assembly mapper cache performance and
 # therefore speed up exon sequence retrieval
  foreach my $gene ( sort { $a->start <=> $b->start } @$genes ) {
    #$self->logger->log_progressbar($progress_id, ++$i, 2);
    #$self->logger->log_progress($num_genes, ++$i, 20, 3, 1);

    if ( $need_project eq 'CHECK' ) {
      # find out whether native coord_system is a common coord_system.
      # if so, you don't need to project.
      # also don't project if no common coord_system present
      if ( $self->highest_common_cs ) {
        my $csid = join( ':',
                         $gene->slice->coord_system_name,
                         $gene->slice->coord_system->version );
        if ( $self->is_common_cs($csid) ) {
          $need_project = 0;
        }
      }
      else {
        $need_project = 0;
      }
    }

    # create lightweigt gene
    my $lgene =
      Bio::EnsEMBL::IdMapping::TinyGene->new_fast( [
                          $gene->dbID,          $gene->stable_id,
                          $gene->version,       $gene->created_date,
                          $gene->modified_date, $gene->start,
                          $gene->end,           $gene->strand,
                          $gene->slice->seq_region_name, $gene->biotype,
                          $gene->analysis->logic_name,
                          ] );

    # build gene caches
    $self->add( 'genes_by_id', $type, $gene->dbID, $lgene );

    # transcripts
    foreach my $tr ( @{ $gene->get_all_Transcripts } ) {
      my $ltr =
        Bio::EnsEMBL::IdMapping::TinyTranscript->new_fast( [
                               $tr->dbID,          $tr->stable_id,
                               $tr->version,       $tr->created_date,
                               $tr->modified_date, $tr->start,
                               $tr->end,           $tr->strand,
                               $tr->length, md5_hex( $tr->spliced_seq ),
                                ] );

      $ltr->biotype( $tr->biotype() );
      $ltr->seq_region_name( $tr->slice->seq_region_name() );
      $lgene->add_Transcript($ltr);

      # build transcript caches
      $self->add( 'transcripts_by_id',      $type, $tr->dbID, $ltr );
      $self->add( 'genes_by_transcript_id', $type, $tr->dbID, $lgene );

      # translation (if there is one)
      if ( my $tl = $tr->translation ) {
        my $ltl =
          Bio::EnsEMBL::IdMapping::TinyTranslation->new_fast( [
                         $tl->dbID,          $tl->stable_id,
                         $tl->version,       $tl->created_date,
                         $tl->modified_date, $tr->dbID,
                         $tr->translate->seq,
                       ] );

        $ltr->add_Translation($ltl);

        $self->add( 'translations_by_id', $type, $tl->dbID, $ltl );

        undef $tl;
      }

      # exons
      foreach my $exon ( @{ $tr->get_all_Exons } ) {
        my $lexon =
          Bio::EnsEMBL::IdMapping::TinyExon->new_fast( [
                         $exon->dbID,
                         $exon->stable_id,
                         $exon->version,
                         $exon->created_date,
                         $exon->modified_date,
                         $exon->start,
                         $exon->end,
                         $exon->strand,
                         $exon->slice->seq_region_name,
                         $exon->slice->coord_system_name,
                         $exon->slice->coord_system->version,
                         $exon->slice->subseq( $exon->start, $exon->end,
                                               $exon->strand ),
                         $exon->phase,
                         $need_project, ] );

        # get coordinates in common coordinate system if needed
        if ($need_project) {
          my @seg = @{
            $exon->project( $self->highest_common_cs,
                            $self->highest_common_cs_version ) };

          if ( scalar(@seg) == 1 ) {
            my $sl = $seg[0]->to_Slice;
            $lexon->common_start( $sl->start );
            $lexon->common_end( $sl->end );
            $lexon->common_strand( $sl->strand );
            $lexon->common_sr_name( $sl->seq_region_name );
          }
        }

        $ltr->add_Exon($lexon);

        $self->add( 'exons_by_id', $type, $exon->dbID, $lexon );
        $self->add_list( 'transcripts_by_exon_id',
                         $type, $exon->dbID, $ltr );

        undef $exon;
      } ## end foreach my $exon ( @{ $tr->get_all_Exons...})

      undef $tr;
    } ## end foreach my $tr ( @{ $gene->get_all_Transcripts...})

    undef $gene;
  } ## end foreach my $gene ( sort { $a...})

  return $num_genes;
} ## end sub build_cache_from_genes


=head2 filter_biotypes

  Arg[1]      : Listref of Bio::EnsEMBL::Genes $genes - the genes to filter
  Example     : my @filtered = @{ $cache->filter_biotypes(\@genes) };

  Description : Filters a list of genes by biotype.  Biotypes are
                taken from the IdMapping configuration parameter
                'biotypes_include' or 'biotypes_exclude'.

                If the configuration parameter 'biotypes_exclude' is
                defined, then rather than returning the genes whose
                biotype is listed in the configuration parameter
                'biotypes_include' the method will return the genes
                whose biotype is *not* listed in the 'biotypes_exclude'
                configuration parameter.

                It is an error to define both these configuration
                parameters.

                The old parameter 'biotypes' is equivalent to
                'biotypes_include'.

  Return type : Listref of Bio::EnsEMBL::Genes (or empty list)
  Exceptions  : none
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub filter_biotypes {
  my ( $self, $genes ) = @_;

  my @filtered;
  my @biotypes;
  my $opt_reverse;

  if ( defined( $self->conf()->param('biotypes_include') ) ||
       defined( $self->conf()->param('biotypes') ) )
  {
    if ( defined( $self->conf()->param('biotypes_exclude') ) ) {
      $self->logger()
        ->error( "You may not use both " .
                 "'biotypes_include' and 'biotypes_exclude' " .
                 "in the configuration" );
    }

    if ( defined( $self->conf()->param('biotypes_include') ) ) {
      @biotypes = $self->conf()->param('biotypes_include');
    }
    else {
      @biotypes = $self->conf()->param('biotypes');
    }
    $opt_reverse = 0;
  }
  else {
    @biotypes    = $self->conf()->param('biotypes_exclude');
    $opt_reverse = 1;
  }

  foreach my $gene ( @{$genes} ) {
    my $keep_gene;

    foreach my $biotype (@biotypes) {
      if ( $gene->biotype() eq $biotype ) {
        if   ($opt_reverse) { $keep_gene = 0 }
        else                { $keep_gene = 1 }
        last;
      }
    }

    if ( defined($keep_gene) ) {
      if ($keep_gene) {
        push( @filtered, $gene );
      }
    }
    elsif ($opt_reverse) {
      push( @filtered, $gene );
    }
  }

  return \@filtered;
} ## end sub filter_biotypes


=head2 add

  Arg[1]      : String $name - a cache name (e.g. 'genes_by_id')
  Arg[2]      : String type - a cache type (e.g. "source.$slice_name")
  Arg[3]      : String $key - key of this entry (e.g. a gene dbID)
  Arg[4]      : Bio::EnsEMBL::IdMappping::TinyFeature $val - value to cache
  Example     : $cache->add('genes_by_id',
                  'source.chromosome:NCBI36:X:1:1000000:1', '1234', $tiny_gene);
  Description : Adds a TinyFeature object to a named cache.
  Return type : Bio::EnsEMBL::IdMapping::TinyFeature
  Exceptions  : thrown on wrong or missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub add {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;
  my $val = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  $self->{'cache'}->{$type}->{$name}->{$key} = $val;

  return $self->{'cache'}->{$type}->{$name}->{$key};
}

=head2 add_list

  Arg[1]      : String $name - a cache name (e.g. 'genes_by_id')
  Arg[2]      : String type - a cache type (e.g. "source.$slice_name")
  Arg[3]      : String $key - key of this entry (e.g. a gene dbID)
  Arg[4]      : List of Bio::EnsEMBL::IdMappping::TinyFeature @val - values
                to cache
  Example     : $cache->add_list('transcripts_by_exon_id',
                  'source.chromosome:NCBI36:X:1:1000000:1', '1234',
                  $tiny_transcript1, $tiny_transcript2);
  Description : Adds a list of TinyFeature objects to a named cache.
  Return type : Listref of Bio::EnsEMBL::IdMapping::TinyFeature objects
  Exceptions  : thrown on wrong or missing arguments
  Caller      : internal
  Status      : At Risk
              : under development

=cut

sub add_list {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;
  my @vals = @_;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  push @{ $self->{'cache'}->{$type}->{$name}->{$key} }, @vals;

  return $self->{'cache'}->{$type}->{$name}->{$key};
}

sub get_by_key {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  # transparently load cache from file unless already loaded
  unless ($self->{'instance'}->{'loaded'}->{"$type"}) {
    $self->read_and_merge($type);
  }

  return $self->{'cache'}->{$type}->{$name}->{$key};
}

sub get_by_name {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  
  # transparently load cache from file unless already loaded
  unless ($self->{'instance'}->{'loaded'}->{$type}) {
    $self->read_and_merge($type);
  }

  return $self->{'cache'}->{$type}->{$name} || {};
}


sub get_count_by_name {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  
  # transparently load cache from file unless already loaded
  unless ($self->{'instance'}->{'loaded'}->{$type}) {
    $self->read_and_merge($type);
  }

  return scalar(keys %{ $self->get_by_name($name, $type) });
}


sub find_common_coord_systems {
  my $self = shift;

  # get adaptors for source db
  my $s_dba = $self->get_DBAdaptor('source');
  my $s_csa = $s_dba->get_CoordSystemAdaptor;
  my $s_sa  = $s_dba->get_SliceAdaptor;

  # get adaptors for target db
  my $t_dba = $self->get_DBAdaptor('target');
  my $t_csa = $t_dba->get_CoordSystemAdaptor;
  my $t_sa  = $t_dba->get_SliceAdaptor;

  # find common coord_systems
  my @s_coord_systems = @{ $s_csa->fetch_all };
  my @t_coord_systems = @{ $t_csa->fetch_all };
  my $found_highest   = 0;

SOURCE:
  foreach my $s_cs (@s_coord_systems) {
    if ( !$s_cs->is_default() ) { next SOURCE }

  TARGET:
    foreach my $t_cs (@t_coord_systems) {
      if ( !$t_cs->is_default() ) { next TARGET }

      if ( $s_cs->name eq $t_cs->name ) {

        # test for identical coord_system version
        if ( $s_cs->version and ( $s_cs->version ne $t_cs->version ) ) {
          next TARGET;
        }

        # test for at least 50% identical seq_regions
        if ( $self->seq_regions_compatible( $s_cs, $s_sa, $t_sa ) ) {
          $self->add_common_cs($s_cs);

          unless ($found_highest) {
            $self->highest_common_cs( $s_cs->name );
            $self->highest_common_cs_version( $s_cs->version );
          }

          $found_highest = 1;

          next SOURCE;
        }
      }
    } ## end foreach my $t_cs (@t_coord_systems)
  } ## end foreach my $s_cs (@s_coord_systems)

  return $found_highest;
} ## end sub find_common_coord_systems


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
  my $t_count = scalar(@$t_seq_regions);
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
    $self->logger->info("Only $equal seq_regions identical for ".$cs->name." ".$cs->version."\n");
    return(0);
  }
  
}


sub check_db_connection {
  my $self = shift;
  my $dbtype = shift;
  
  my $err = 0;
  
  eval {
    my $dba = $self->get_DBAdaptor($dbtype);
    $dba->dbc->connect;
  };
  
  if ($@) {
    $self->logger->warning("Can't connect to $dbtype db: $@\n");
    $err++;
  } else {
    $self->logger->debug("Connection to $dbtype db ok.\n");
    $self->{'_db_conn_ok'}->{$dbtype} = 1;
  }

  return $err;
}

  
sub check_db_read_permissions {
  my $self = shift;
  my $dbtype = shift;

  # skip this check if db connection failed (this prevents re-throwing
  # exceptions).
  return 1 unless ($self->{'_db_conn_ok'}->{$dbtype});
  
  my $err = 0;
  my %privs = %{ $self->get_db_privs($dbtype) };
  
  unless ($privs{'SELECT'} or $privs{'ALL PRIVILEGES'}) {
    $self->logger->warning("User doesn't have read permission on $dbtype db.\n");
    $err++;
  } else {
    $self->logger->debug("Read permission on $dbtype db ok.\n");
  }

  return $err;
}

  
sub check_db_write_permissions {
  my $self = shift;
  my $dbtype = shift;
  
  # skip this check if db connection failed (this prevents re-throwing
  # exceptions).
  return 1 unless ($self->{'_db_conn_ok'}->{$dbtype});
  
  my $err = 0;

  unless ($self->do_upload) {
    $self->logger->debug("No uploads, so write permission on $dbtype db not required.\n");
    return $err;
  }

  my %privs = %{ $self->get_db_privs($dbtype) };

  unless ($privs{'INSERT'} or $privs{'ALL PRIVILEGES'}) {
    $self->logger->warning("User doesn't have write permission on $dbtype db.\n");
    $err++;
  } else {
    $self->logger->debug("Write permission on $dbtype db ok.\n");
  }

  return $err;
}


sub do_upload {
  my $self = shift;

  if ($self->conf->param('dry_run') or
    ! ($self->conf->param('upload_events') or
       $self->conf->param('upload_stable_ids') or
       $self->conf->param('upload_archive'))) {
    return 0;
  } else {
    return 1;
  }
}   


sub get_db_privs {
  my ( $self, $dbtype ) = @_;

  my %privs = ();
  my $rs;

  # get privileges from mysql db
  eval {
    my $dbc = $self->get_DBAdaptor($dbtype)->dbc();
    my $sql = qq(SHOW GRANTS FOR ) . $dbc->username();
    my $sth = $dbc->prepare($sql);
    $sth->execute();
    $rs = $sth->fetchall_arrayref();
    #$sth->finish();
  };

  if ($@) {
    $self->logger->warning(
      "Error obtaining privileges from $dbtype db: $@\n");
    return {};
  }

  # parse the output
  foreach my $r ( map { $_->[0] } @{$rs} ) {
    $r =~ s/GRANT (.*) ON .*/$1/i;
    foreach my $p ( split( ',', $r ) ) {
      # trim leading and trailing whitespace
      $p =~ s/^\s+//;
      $p =~ s/\s+$//;
      $privs{ uc($p) } = 1;
    }
  }

  return \%privs;
} ## end sub get_db_privs


sub check_empty_tables {
  my $self = shift;
  my $dbtype = shift;
  
  # skip this check if db connection failed (this prevents re-throwing
  # exceptions).
  return 1 unless ($self->{'_db_conn_ok'}->{$dbtype});
  
  my $err = 0;
  my $c = 0;

  if ($self->conf->param('no_check_empty_tables') or !$self->do_upload) {
    $self->logger->debug("Won't check for empty stable ID and archive tables in $dbtype db.\n");
    return $err;
  }

  eval {
    my @tables =
      qw(
      gene_stable_id
      transcript_stable_id
      translation_stable_id
      exon_stable_id
      stable_id_event
      mapping_session
      gene_archive
      peptide_archive
    );

    my $dba = $self->get_DBAdaptor($dbtype);
    foreach my $table (@tables) {
      if ( $table =~ /^([^_]+)_stable_id/ ) {
        $table = $1;
        if ( $c =
             $self->fetch_value_from_db(
               $dba,
               "SELECT COUNT(*) FROM $table WHERE stable_id IS NOT NULL"
             ) )
        {
          $self->logger->warning(
                        "$table table in $dbtype db has $c stable IDs.\n");
          $err++;
        }
      }
      else {
        if ( $c =
             $self->fetch_value_from_db(
                                     $dba, "SELECT COUNT(*) FROM $table"
             ) )
        {
          $self->logger->warning(
                        "$table table in $dbtype db has $c entries.\n");
          $err++;
        }
      }
    } ## end foreach my $table (@tables)
  };

  if ($@) {
    $self->logger->warning(
"Error retrieving stable ID and archive table row counts from $dbtype db: $@\n"
    );
    $err++;
  }
  elsif ( !$err ) {
    $self->logger->debug(
         "All stable ID and archive tables in $dbtype db are empty.\n");
  }
  return $err;
}


sub check_sequence {
  my ( $self, $dbtype ) = @_;

  # skip this check if db connection failed (this prevents re-throwing
  # exceptions).
  return 1 unless ( $self->{'_db_conn_ok'}->{$dbtype} );

  my $err = 0;
  my $c   = 0;

  eval {
    my $dba = $self->get_DBAdaptor($dbtype);
    unless ( $c =
             $self->fetch_value_from_db(
                               $dba->dnadb(), "SELECT COUNT(*) FROM dna"
             ) )
    {
      $err++;
    }
  };

  if ($@) {
    $self->logger->warning(   "Error retrieving dna table row count "
                            . "from $dbtype database: $@\n" );
    $err++;
  } elsif ($err) {
    $self->logger->warning("No sequence found in $dbtype database.\n");
  } else {
    $self->logger->debug(
                ucfirst($dbtype) . " db has sequence ($c entries).\n" );
  }

  return $err;
} ## end sub check_sequence


sub check_meta_entries {
  my $self = shift;
  my $dbtype = shift;
  
  # skip this check if db connection failed (this prevents re-throwing
  # exceptions).
  return 1 unless ($self->{'_db_conn_ok'}->{$dbtype});
  
  my $err = 0;
  my $assembly_default;
  my $schema_version;
  
  eval {
    my $dba = $self->get_DBAdaptor($dbtype);
    $assembly_default = $self->fetch_value_from_db($dba,
      qq(SELECT meta_value FROM meta WHERE meta_key = 'assembly.default'));
    $schema_version = $self->fetch_value_from_db($dba,
      qq(SELECT meta_value FROM meta WHERE meta_key = 'schema_version'));
  };
  
  if ($@) {
    $self->logger->warning("Error retrieving dna table row count from $dbtype db: $@\n");
    return ++$err;
  }
  
  unless ($assembly_default) {
    $self->logger->warning("No meta.assembly.default value found in $dbtype db.\n");
    $err++;
  } else {
    $self->logger->debug("meta.assembly.default value found ($assembly_default).\n");
  }

  unless ($schema_version) {
    $self->logger->warning("No meta.schema_version value found in $dbtype db.\n");
    $err++;
  } else {
    $self->logger->debug("meta.schema_version value found ($schema_version).\n");
  }

  return $err;
}


sub fetch_value_from_db {
  my ( $self, $dba, $sql ) = @_;

  assert_ref( $dba, 'Bio::EnsEMBL::DBSQL::DBAdaptor' );

  if ( !defined($sql) ) {
    throw("Need an SQL statement to execute.\n");
  }

  my $sth = $dba->dbc->prepare($sql);
  $sth->execute();

  my ($c) = $sth->fetchrow_array;
  return $c;
}

sub get_DBAdaptor {
  my ( $self, $prefix ) = @_;

  unless ( $self->{'_dba'}->{$prefix} ) {
    # connect to database
    my $dba =
      new Bio::EnsEMBL::DBSQL::DBAdaptor(
                       -host   => $self->conf->param("${prefix}host"),
                       -port   => $self->conf->param("${prefix}port"),
                       -user   => $self->conf->param("${prefix}user"),
                       -pass   => $self->conf->param("${prefix}pass"),
                       -dbname => $self->conf->param("${prefix}dbname"),
                       -group  => $prefix, );

    if ( !defined( $self->conf->param("${prefix}host_dna") ) ) {
      # explicitely set the dnadb to itself - by default the Registry
      # assumes a group 'core' for this now
      $dba->dnadb($dba);
    } else {
      my $dna_dba =
        new Bio::EnsEMBL::DBSQL::DBAdaptor(
                   -host   => $self->conf->param("${prefix}host_dna"),
                   -port   => $self->conf->param("${prefix}port_dna"),
                   -user   => $self->conf->param("${prefix}user_dna"),
                   -pass   => $self->conf->param("${prefix}pass_dna"),
                   -dbname => $self->conf->param("${prefix}dbname_dna"),
                   -group  => $prefix, );
      $dba->dnadb($dna_dba);
    }

    $self->{'_dba'}->{$prefix} = $dba;
  } ## end unless ( $self->{'_dba'}->...)

  return $self->{'_dba'}->{$prefix};
} ## end sub get_DBAdaptor


sub get_production_DBAdaptor() {
  my ($self) = @_;
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                       -host   => $self->conf->param("productionhost"),
                       -port   => $self->conf->param("productionport"),
                       -user   => $self->conf->param("productionuser"),
                       -pass   => $self->conf->param("productionpass"),
                       -dbname => $self->conf->param("productiondbname"));
  return $dba;
}


sub cache_file_exists {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  my $cache_file = $self->cache_file($type);

  if (-e $cache_file) {
    $self->logger->info("Cache file found for $type.\n", 2);
    $self->logger->debug("Will read from $cache_file.\n", 2);
    return 1;
  } else {
    $self->logger->info("No cache file found for $type.\n", 2);
    $self->logger->info("Will build cache from db.\n", 2);
    return 0;
  }
}


sub cache_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  return $self->dump_path."/$type.object_cache.ser";
}


sub instance_file {
  my $self = shift;

  return $self->dump_path."/cache_instance.ser";
}


sub dump_path {
  my $self = shift;

  $self->{'dump_path'} ||= path_append($self->conf->param('basedir'), 'cache');

  return $self->{'dump_path'};
}


sub write_all_to_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  my $size = 0;
  $size += $self->write_to_file($type);
  $size += $self->write_instance_to_file;

  return parse_bytes($size);
}


sub write_to_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  unless ($self->{'cache'}->{$type}) {
    $self->logger->warning("No features found in $type. Won't write cache file.\n");
    return;
  }

  my $cache_file = $self->cache_file($type);

  eval { nstore($self->{'cache'}->{$type}, $cache_file) };
  if ($@) {
    throw("Unable to store $cache_file: $@\n");
  }

  my $size = -s $cache_file;
  return $size;
}


sub write_instance_to_file {
  my $self = shift;

  my $instance_file = $self->instance_file;

  eval { nstore($self->{'instance'}, $instance_file) };
  if ($@) {
    throw("Unable to store $instance_file: $@\n");
  }

  my $size = -s $instance_file;
  return $size;
}


sub read_from_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  my $cache_file = $self->cache_file($type);

  if (-s $cache_file) {
    
    #$self->logger->info("Reading cache from file...\n", 0, 'stamped');
    #$self->logger->info("Cache file $cache_file.\n", 1);
    eval { $self->{'cache'}->{$type} = retrieve($cache_file); };
    if ($@) {
      throw("Unable to retrieve cache: $@");
    }
    #$self->logger->info("Done.\n", 0, 'stamped');

  } else {
    $self->logger->warning("Cache file $cache_file not found or empty.\n");
  }


  return $self->{'cache'}->{$type};
}


sub read_and_merge {
  my $self = shift;
  my $dbtype = shift;

  unless ($dbtype eq 'source' or $dbtype eq 'target') {
    throw("Db type must be 'source' or 'target'.");
  }

  # read cache from single or multiple files, depending on caching strategy
  my $cache_method = $self->cache_method;
  if ($cache_method eq 'ALL') {
    $self->read_from_file("$dbtype.ALL");
  } elsif ($cache_method eq 'BY_SEQ_REGION') {
    foreach my $slice_name (@{ $self->slice_names($dbtype) }) {
      $self->read_from_file("$dbtype.$slice_name");
    }
  } else {
    throw("Unknown cache method: $cache_method.");
  }

  $self->merge($dbtype);

  # flag as being loaded
  $self->{'instance'}->{'loaded'}->{$dbtype} = 1;
}


sub merge {
  my $self = shift;
  my $dbtype = shift;

  unless ($dbtype eq 'source' or $dbtype eq 'target') {
    throw("Db type must be 'source' or 'target'.");
  }

  foreach my $type (keys %{ $self->{'cache'} || {} }) {
    next unless ($type =~ /^$dbtype/);

    foreach my $name (keys %{ $self->{'cache'}->{$type} || {} }) {
    
      foreach my $key (keys %{ $self->{'cache'}->{$type}->{$name} || {} }) {
        if (defined $self->{'cache'}->{$dbtype}->{$name}->{$key}) {
          # warning("Duplicate key in cache: $name|$dbtype|$key. Skipping.\n");
        } else {
          $self->{'cache'}->{$dbtype}->{$name}->{$key} =
            $self->{'cache'}->{$type}->{$name}->{$key};
        }

        delete $self->{'cache'}->{$type}->{$name}->{$key};
      }
      
      delete $self->{'cache'}->{$type}->{$name};
    }
    
    delete $self->{'cache'}->{$type};

  }
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
  my $self   = shift;
  my $dbtype = shift;

  throw("You must provide a db type (source|target).") unless $dbtype;

  my $dba = $self->get_DBAdaptor($dbtype);
  my $sa  = $dba->get_SliceAdaptor;

  my @slice_names = ();

  if ( $self->conf->param('chromosomes') ) {
    # Fetch the specified chromosomes.
    foreach my $chr ( $self->conf->param('chromosomes') ) {
      my $slice = $sa->fetch_by_region( 'chromosome', $chr );
      push @slice_names, $slice->name;
    }

  }
  elsif ( $self->conf->param('region') ) {
    # Fetch the slices on the specified regions.  Don't use
    # SliceAdaptor->fetch_by_name() since this will fail if assembly
    # versions are different for source and target db.
    my ( $cs, $version, $name, $start, $end, $strand ) =
      split( /:/, $self->conf->param('region') );

    my $slice = $sa->fetch_by_region( $cs, $name, $start, $end );

    push @slice_names, $slice->name;

  }
  else {
    # Fetch all slices that have genes on them.
    my $ga = $dba->get_GeneAdaptor;
    my $sa = $dba->get_SliceAdaptor;

    foreach my $srid ( @{ $ga->list_seq_region_ids } ) {
      my $slice = $sa->fetch_by_seq_region_id($srid);
      my $slices = $sa->fetch_by_region_unique( $slice->coord_system_name(), $slice->seq_region_name() );

      push( @slice_names, map { $_->name() } @{$slices} );
    }
  }

  return \@slice_names;
} ## end sub slice_names


sub logger {
  my $self = shift;
  $self->{'logger'} = shift if (@_);
  return $self->{'logger'};
}

sub conf {
  my $self = shift;
  $self->{'conf'} = shift if (@_);
  return $self->{'conf'};
}


sub cache_method {
  my $self = shift;
  $self->{'instance'}->{'cache_method'} = shift if (@_);
  return $self->{'instance'}->{'cache_method'};
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

  my $csid = join(':', $cs->name, $cs->version);

  $self->{'instance'}->{'ccs'}->{$csid} = 1;
}


sub is_common_cs {
  my $self = shift;
  my $csid = shift;

  return $self->{'instance'}->{'ccs'}->{$csid};
}


1;

