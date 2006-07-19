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
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Storable qw(nfreeze thaw nstore retrieve);


=head2 new

  Arg[1]      : 
  Example     : 
  Description : constructor
  Return type : Bio::EnsEMBL::IdMapping::Cache object
  Exceptions  : 
  Caller      : general

=cut

sub new {
  my $class = shift;

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


=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 
  Status      :

=cut

sub cache_file_exists {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type (old|new).") unless $type;

  my $cache_path = $self->cache_path($type);

  if (-s $cache_path) {
    $self->logger->log("Cache file found. Will read from $cache_path.\n", 1);
    return(1);
  } else {
    $self->logger->log("No cache file found. Will build cache from db.\n", 1);
    return(0);
  }
}


sub cache_path {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type (old|new).") unless $type;

  my $cache_path = ($self->conf->param('dumppath') || '.')."/$type.".
    ($self->conf->param('cachefile') || 'object_cache.ser');

  return $cache_path;
}


sub write_to_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type (old|new).") unless $type;

  # create dump directory if it doesn't exist
  if (my $dump_path = $self->conf->param('dumppath')) {
    unless (-d $dump_path) {
      $self->logger->log("Cache directory $dump_path doesn't exist. Will create it.\n", 1);
      system("mkdir -p $dump_path") == 0 or
        throw("Unable to create directory $dump_path.\n");
    }
  }
  
  my $cache_path = $self->cache_path($type);

  $self->logger->log("Will dump to $cache_path.\n", 1);
  $self->logger->log_stamped("Dumping cache...\n", 1);

  eval { nstore($self->{'_cache'}->{$type}, $cache_path) };
  if ($@) {
    throw("Unable to store $cache_path: $@\n");
  }
  my $size = -s $cache_path;
  $size = parse_bytes($size);
  
  $self->logger->log_stamped("Done (cache file is $size).\n", 1);
}


sub read_from_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type (old|new).") unless $type;

  my $cache_path = $self->cache_path($type);

  unless (-s $cache_path) {
    throw("No valid cache file found at $cache_path.");
  }

  $self->logger->log("Reading cache from $cache_path...\n");
  eval { $self->{'_cache'}->{$type} = retrieve($cache_path); };
  if ($@) {
    throw("Unable to retrieve cache: $@");
  }
  $self->logger->log_stamped("Done.\n");

  return $self->{'_cache'}->{$type};
}


sub fetch_genes {
  my $self = shift;
  my $dba = shift;

  $self->logger->log_stamped("Fetching genes...\n", 1);

  unless ($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    throw("You must provide a database adaptor.");
  }

  my $sa = $dba->get_SliceAdaptor;
  my $ga = $dba->get_GeneAdaptor;

  #
  # fetch genes, depending on filters to apply
  #
  my @all_genes = ();
  my $genes = [];

  if ($self->conf->param('chromosomes')) {
    # filter by chromosome
    $self->logger->log("Filtering by chromosome: ", 3);
    
    foreach my $chr ($self->conf->param('chromosomes')) {
      $self->logger->log("$chr ");
      my $slice = $sa->fetch_by_region('chromosome', $chr);
      push @all_genes, @{ $slice->get_all_Genes(undef, undef, 1) };
    }
    $self->logger->log("\n");
    
  } elsif ($self->conf->param('region')) {
    # filter by region (specific slice)
    $self->logger->log("Filtering by region: ".$self->conf->param('region')."\n", 3);

    my $slice = $sa->fetch_by_name($self->conf->param('region'));
    @all_genes = @{ $slice->get_all_Genes(undef, undef, 1) };
    
  } else {
    # fetch all genes
    @all_genes = @{ $ga->fetch_all };
  }

  # filter by biotype
  if ($self->conf->param('biotypes')) {
    $self->logger->log("Filtering by biotype: ", 3);

    foreach my $biotype ($self->conf->param('biotypes')) {
      $self->logger->log("$biotype ");
      push @$genes, grep { $_->biotype eq $biotype } @all_genes;
    }

    $self->logger->log("\n");

  } else {
    $genes = \@all_genes;
  }

  $self->logger->log_stamped("Done loading ".scalar(@$genes)." genes.\n\n", 1);

  return $genes;
}


sub build_cache_from_genes {
  my $self = shift;
  my $type = shift;
  my $genes = shift;
  
  $self->logger->log_stamped("Indexing...\n", 1);

  throw("You must provide a type (old|new).") unless $type;
  throw("You must provide a listref of genes.") unless (ref($genes) eq 'ARRAY');

  my $i = 0;
  my $num_genes = scalar(@$genes);
  $self->logger->init_progressbar('index_genes', $num_genes);

#use Data::Dumper;
#$Data::Dumper::Indent = 1;

  foreach my $gene (@$genes) {
    #$self->logger->log_progressbar('index_genes', ++$i, 2);
    $self->logger->log_progress($num_genes, ++$i, 20, 3, 1);
    
    # build gene caches
    $self->add($type, 'genes_by_id', $gene->dbID, $gene);
    $self->add($type, 'genes_by_stable_id', $gene->stable_id, $gene);
    
    # transcripts
    foreach my $tr (@{ $gene->get_all_Transcripts }) {
      # build transcript caches
      $self->add($type, 'transcripts_by_id', $tr->dbID, $tr);
      $self->add($type, 'transcripts_by_stable_id', $tr->stable_id, $tr);
      $self->add($type, 'genes_by_transcript_id', $tr->dbID, $gene);

      # translation (if there is one)
      if (my $tl = $tr->translation) {
        $self->add($type, 'translations_by_id', $tl->dbID, $tl);
        $self->add($type, 'translations_by_stable_id', $tl->stable_id, $tl);
        $self->add($type, 'translations_by_transcript_id', $tr->dbID, $tl);
      }

      # exons
      foreach my $exon (@{ $tr->get_all_Exons }) {
        # force sequence lazy-loading
        $exon->{'_seq_cache_dump'} =
          $exon->slice->subseq($exon->start, $exon->end, $exon->strand);
        
        $self->add($type, 'exons_by_id', $exon->dbID, $exon);
        $self->add($type, 'genes_by_exon_id', $exon->dbID, $gene);
        $self->add_list($type, 'transcripts_by_exon_id', $exon->dbID, $tr);
      }
    }
  }

  foreach my $gene (@$genes) {
    #
    # delete all unwanted data from objects
    #
    $gene->adaptor(undef);
    $gene->analysis(undef);
    $gene->display_xref(undef);
    $gene->slice->adaptor(undef);
    $gene->slice->coord_system->adaptor(undef);

    foreach my $tr (@{ $gene->get_all_Transcripts }) {
      $tr->adaptor(undef);
      $tr->display_xref(undef);
      $tr->analysis(undef);
      
      if (my $tl = $tr->translation) {
        $tl->adaptor(undef);
      }

      foreach my $exon (@{ $tr->get_all_Exons }) {
        $exon->adaptor(undef);
        $exon->analysis(undef);
      }
    }

    #warn Data::Dumper::Dumper($gene) if ($i < 2);;
  }

  $self->logger->log_stamped("Done building the index:\n", 2);
  foreach my $name (qw(genes transcripts translations exons)) {
    my $num = scalar(keys %{ $self->get_by_name($type, "${name}_by_id") });
    $self->logger->log(sprintf("%12.0f %-20s\n", $num, $name), 3);
  }
  $self->logger->log("\n");
}


sub add {
  my $self = shift;
  my $type = shift;
  my $name = shift;
  my $key = shift;
  my $val = shift;

  throw("You must provide a cache type (old|new).") unless $type;
  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  $self->{'_cache'}->{$type}->{$name}->{$key} = $val;

  return $self->{'_cache'}->{$type}->{$name}->{$key};
}

sub add_list {
  my $self = shift;
  my $type = shift;
  my $name = shift;
  my $key = shift;
  my @vals = @_;

  throw("You must provide a cache type (old|new).") unless $type;
  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  push @{ $self->{'_cache'}->{$type}->{$name}->{$key} }, @vals;

  return $self->{'_cache'}->{$type}->{$name}->{$key};
}

sub get_by_key {
  my $self = shift;
  my $type = shift;
  my $name = shift;
  my $key = shift;

  throw("You must provide a cache type (old|new).") unless $type;
  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  return $self->{'_cache'}->{$type}->{$name}->{$key};
}

sub get_by_name {
  my $self = shift;
  my $type = shift;
  my $name = shift;

  throw("You must provide a cache type (old|new).") unless $type;
  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  
  return $self->{'_cache'}->{$type}->{$name} || {};
}


sub flush_by_type {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type (old|new).") unless $type;

  undef $self->{'_cache'}->{$type};
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


1;

