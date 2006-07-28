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


sub write_all_to_file {
  my $self = shift;
  my $type = shift;

  throw("You must provide a cache type.") unless $type;

  foreach my $name (@{ $self->cache_names }) {
    $self->write_to_file($name, $type);
  }
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

  eval { nstore($self->{'_cache'}->{$name}->{$type}, $cache_file) };
  if ($@) {
    throw("Unable to store $cache_file: $@\n");
  }

  my $size = -s $cache_file;
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
  eval { $self->{'_cache'}->{$name}->{$type} = retrieve($cache_file); };
  if ($@) {
    throw("Unable to retrieve cache: $@");
  }
  $self->logger->log_stamped("Done.\n");

  return $self->{'_cache'}->{$name}->{$type};
}


sub slice_names {
  my $self = shift;
  my $dbtype = shift;

  throw("You must provide a db type (old|new).") unless $dbtype;

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

  # build cache
  my $type = "$dbtype.$slice_name";
  $self->build_cache_from_genes($type, $genes);
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
  
  throw("You must provide a type.") unless $type;
  throw("You must provide a listref of genes.") unless (ref($genes) eq 'ARRAY');

  #my $i = 0;
  #my $num_genes = scalar(@$genes);
  #$self->logger->init_progressbar('index_genes', $num_genes);

  foreach my $gene (@$genes) {
    #$self->logger->log_progressbar('index_genes', ++$i, 2);
    #$self->logger->log_progress($num_genes, ++$i, 20, 3, 1);
    
    # create lightweigt gene
    my $lgene = Bio::EnsEMBL::IdMapping::TinyFeature->new_fast([
        'g',
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
      my $ltr = Bio::EnsEMBL::IdMapping::TinyFeature->new_fast([
          'tr',
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
        my $ltl = Bio::EnsEMBL::IdMapping::TinyFeature->new_fast([
            'tl',
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
        my $lexon = Bio::EnsEMBL::IdMapping::TinyFeature->new_fast([
            'e',
            $exon->dbID,
            $exon->stable_id,
            $exon->start,
            $exon->end,
            $exon->strand,
            $exon->slice->seq_region_name,
            $exon->slice->coord_system_name,
            $exon->slice->coord_system->version,
            #$exon->slice->subseq($exon->start, $exon->end, $exon->strand),
        ]);

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

  $self->{'_cache'}->{$name}->{$type}->{$key} = $val;

  return $self->{'_cache'}->{$name}->{$type}->{$key};
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

  push @{ $self->{'_cache'}->{$name}->{$type}->{$key} }, @vals;

  return $self->{'_cache'}->{$name}->{$type}->{$key};
}

sub get_by_key {
  my $self = shift;
  my $name = shift;
  my $type = shift;
  my $key = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  throw("You must provide a cache key (e.g. a gene dbID).") unless $key;

  return $self->{'_cache'}->{$name}->{$type}->{$key};
}

sub get_by_name {
  my $self = shift;
  my $name = shift;
  my $type = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;
  throw("You must provide a cache type.") unless $type;
  
  return $self->{'_cache'}->{$name}->{$type} || {};
}


sub merge {
  my $self = shift;
  my $name = shift;

  throw("You must provide a cache name (e.g. genes_by_id.") unless $name;

  foreach my $type (keys %{ $self->{'_cache'}->{$name} || {} }) {
    (my $merged_type = $type) =~ s/^(\w+)\..+/$1/;
    
    foreach my $key (keys %{ $self->{'_cache'}->{$name}->{$type} || {} }) {
      if (defined $self->{'_cache'}->{$name}->{$merged_type}->{$key}) {
        warning("Duplicate key in cache: $name|$merged_type|$key. Skipping.\n");
      } else {
        $self->{'_cache'}->{$name}->{$merged_type}->{$key} =
          $self->{'_cache'}->{$name}->{$type}->{$key};
      }

      delete $self->{'_cache'}->{$name}->{$type}->{$key};
    }
    
    delete $self->{'_cache'}->{$name}->{$type};

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

