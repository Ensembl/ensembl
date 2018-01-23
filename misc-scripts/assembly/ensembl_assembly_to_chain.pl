#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=pod

=head1 CHAIN FILES

This code writes chain files from Ensembl assembly records. Chain files are defined at

  http://genome.ucsc.edu/goldenPath/help/chain.html

They can be used with liftover (http://genome.ucsc.edu/cgi-bin/hgLiftOver) and crossmap (http://crossmap.sourceforge.net). The file format says the following:

  chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
  size dt dq
  size

You can see this in the the example chain files for human mappings http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz and http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz.

dt & dq are meant to be the differences between the end of this current block and the start of the next in the target and query. The final line is just a size to indicate no further offset.

tSize and qSize are the total sizes of the target and reference sequences.

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::IO qw/work_with_file gz_work_with_file/;
use Scalar::Util qw/looks_like_number/;
use feature qw/say/;
use File::Path qw/mkpath/;
use File::Spec;
use JSON;

my $COMPRESS = 0;
my $UCSC = 0;
my %ucsc_name_cache;

sub get_options {
  my ($db_name, $db_host, $db_user, $db_pass, $db_port, $help, @species, $group, $release, $dir, $compress, $ucsc);
  $db_port = 3306;
  $group = 'core';
  $compress=0;
  $ucsc=0;
  $dir = File::Spec->curdir();

  GetOptions(
    "db_name|dbname|database=s"       => \$db_name,
    "db_host|dbhost|host=s"           => \$db_host,
    "db_user|dbuser|user|username=s"  => \$db_user,
    "db_pass|dbpass|pass|password=s"  => \$db_pass,
    "db_port|dbport|port=s"           => \$db_port,
    "species=s@"                      => \@species,
    "version|release=i"               => \$release,
    "directory|dir=s"                 => \$dir,
    "compress|gzip|gz!"               => \$compress,
    "ucsc!"                           => \$ucsc,
    "help|h!"                         => \$help,
  );

  my %args = (
    -HOST => $db_host, -PORT => $db_port, 
    -USER => $db_user
  );
  $args{-PASS} = $db_pass if $db_pass;

  # must be a DB we need
  if($db_name) {
    Bio::EnsEMBL::DBSQL::DBAdaptor->new(
      %args,
      -DBNAME => $db_name,
      -SPECIES => $db_name,
      -GROUP => 'core'
    );
  }
  else {
    Bio::EnsEMBL::Registry->load_registry_from_db(
      %args,
      -DB_VERSION => $release,
    );
  }
  
  $COMPRESS = $compress;
  $UCSC = $ucsc;
  my @dbas;
  my %final_dbas;
    
  if(@species) {
    say "Working against a restricted species list";
    foreach my $s (@species) {
      my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($s, $group);
      die "Cannot find a DBAdaptor for the species ${s}" unless $dba;
      push(@dbas, $dba);
    }
  }
  else {
    say "Dumping chain file for all available species";
    @dbas = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors(-GROUP => 'core')};
  }
  
  foreach my $dba (@dbas) {
    my $entries = get_liftover_mappings($dba);
    my $prod_name = $dba->get_MetaContainer()->get_production_name();
    $prod_name //= $dba->species();
    $final_dbas{$prod_name} = $dba if @{$entries};
    $dba->dbc->disconnect_if_idle();
  }
  
  # Return including sorting the species by name as it's just nicer that way
  return ($dir, map { $final_dbas{$_} } sort keys %final_dbas);
}

run();

sub run {
  my ($dir, @dbas) = get_options();
  foreach my $dba (@dbas) {
    run_on_dba($dir, $dba);
    $dba->dbc->disconnect_if_idle;
  }
  return;
}

sub run_on_dba { 
  my ($dir, $core_dba) = @_;
  my $prod_name = $core_dba->get_MetaContainer->get_production_name();
  $prod_name //= $core_dba->species();
  say "Processing ${prod_name}";
  my $liftovers = get_liftover_mappings($core_dba);
  foreach my $mappings (@{$liftovers}) {
    my ($asm_cs, $cmp_cs) = @{$mappings};
    say "\tWorking with $asm_cs to $cmp_cs";
    say "\t\tFetching mappings";
    my $asm_to_cmp_mappings = get_assembly_mappings($core_dba, $asm_cs, $cmp_cs);
    write_mappings($dir, $asm_cs, $cmp_cs, $prod_name, $asm_to_cmp_mappings, $core_dba);
    say "\t\tFetching reverse mappings";
    my $cmp_to_asm_mappings = get_reverse_assembly_mappings($core_dba, $asm_cs, $cmp_cs);
    write_mappings($dir, $cmp_cs, $asm_cs, $prod_name, $cmp_to_asm_mappings, $core_dba);
    say "\t\tFinished";
  }
  return;
}

#Parse mapping keys like chromosome:GRCh37#chromosome:NCBI36 to ['GRCh37','NCBI36']
sub get_liftover_mappings {
  my ($core_dba) = @_;
  my $mappings = $core_dba->get_MetaContainer()->list_value_by_key('liftover.mapping');
  my %unique_mappings;
  foreach my $mapping (@{$mappings}) {
      die "Can't parse mapping string for coord system version '${mapping}'!\n"
          unless $mapping =~ /.+:(.+)#.+:(.+)/;
      my $version_1 = $1;
      my $version_2 = $2;
      $unique_mappings{$version_1. ":". $version_2} = [$version_1, $version_2];
  }
  return [values %unique_mappings];
}

sub write_mappings {
  my ($dir, $source_cs, $target_cs, $prod_name, $mappings, $core_dba) = @_;
  my $file = "${source_cs}_to_${target_cs}.chain";
  $file .= '.gz' if $COMPRESS;
  my $target_dir = File::Spec->catdir($dir, $prod_name);
  mkpath($target_dir);
  my $path = File::Spec->catfile($target_dir, $file);
  say "\t\tBuilding chain mappings";
  my $chains = build_chain_mappings($mappings);
  my $ucsc_chains = [];
  if($UCSC) {
    say "\t\tBuilding UCSC mappings";
    $ucsc_chains = create_ucsc_chains($core_dba, $prod_name, $chains);
  }
  say "\t\tWriting mappings to $path";
  
  my $writer = sub {
    my ($fh) = @_;
    print_chains($fh, $chains);
    print_chains($fh, $ucsc_chains) if @{$ucsc_chains};
  };
  
  if($COMPRESS) {
    gz_work_with_file($path, 'w', $writer);
  }
  else {
    work_with_file($path, 'w', $writer);
  }
  
  return;
}

sub build_chain_mappings {
  my ($assembly_mappings) = @_;
  my @chain_mappings;
  
  my ($t_name, $t_size, $t_strand, $t_start, $t_end);
  my ($q_name, $q_size, $q_strand, $q_start, $q_end);
  my $chain_id = 1;
  my @chain_gaps;
  
  my $length = scalar(@{$assembly_mappings});
  for (my $i = 0; $i < $length; $i++) {

    my $current = $assembly_mappings->[$i];
    my $next = ($i+1 != $length) ? $assembly_mappings->[$i+1] : undef;

    my $ori = $current->{ori};
    my ($asm_diff, $cmp_diff);
    if($next) {
      $asm_diff = ($next->{asm_start} - $current->{asm_end})-1;
      # Rev strands means the next cmp region has a lower start than the 
      # current end (because it's running in reverse). Rember length in 1-based
      # coords always is (high-low)-1
      $cmp_diff = ($ori == 1) ? ($next->{cmp_start} - $current->{cmp_end})-1 : ($current->{cmp_start} - $next->{cmp_end})-1;
    }

    if(! $t_name) {
      # Reset variables to current
      @chain_gaps = ();
      $chain_id++;
      ($t_name, $t_size, $t_strand, $t_start) = ($current->{asm_name}, $current->{asm_length}, 1, $current->{asm_start});
      ($q_name, $q_size, $q_strand) = ($current->{cmp_name}, $current->{cmp_length}, $current->{ori});
      $q_start = ($ori == 1) ? $current->{cmp_start} : $current->{cmp_end};
    }

    # Block that decides we need to start a new chain definition
    #
    # Can mean we have run out of mappings, we are into a new chromsome (both source and target) 
    # or strand has swapped.
    # Final reason is we've had an out-of-order meaning a negative gap was produced 
    # (we're going backwards). In any situation this means the chain is finished
    if( ! defined $next || 
        $t_name ne $next->{asm_name} ||
        $q_name ne $next->{cmp_name} || # we can switch target chromosomes. e.g. cross chromsome mappings in mouse NCBI37->GRCm38
        $ori != $next->{ori} ||
        $asm_diff < 0 ||
        $cmp_diff < 0) {
      # Add the last gap on which is just the length of this alignment
      push(@chain_gaps, [$current->{length}]);
      # Set the ends of the chain since this is the last block
      $t_end = $current->{asm_end};
      $q_end = ($ori == 1) ? $current->{cmp_end} : $current->{cmp_start};

      #If strand was negative we need to represent all data as reverse complemented regions
      if($q_strand == -1) {
        # $t_start = ($t_size - $t_start)+1;
        # $t_end = ($t_size - $t_end)+1;
        $q_start = ($q_size - $q_start)+1;
        $q_end = ($q_size - $q_end)+1;
      }
      # Convert to UCSC formats (0-based half-open intervals and +/- strands)
      $t_start--;
      $q_start--;
      $t_strand = ($t_strand == 1) ? '+' : '-';
      $q_strand = ($q_strand == 1) ? '+' : '-';

      #Store the chain
      my $chain_score = 1;
      push(@chain_mappings, {
        header => ['chain', $chain_score, $t_name, $t_size, $t_strand, $t_start, $t_end, $q_name, $q_size, $q_strand, $q_start, $q_end, $chain_id],
        gaps => [@chain_gaps]
      });

      if(! defined $next) {
        last;
      }

      # Clear variables
      ($t_name, $t_size, $t_strand, $t_start) = ();
      ($q_name, $q_size, $q_strand, $q_start) = ();
    }

    push(@chain_gaps, [$current->{length}, $asm_diff, $cmp_diff]);
  }

  return \@chain_mappings;
}

sub get_assembly_mappings {
  my ($dba, $asm_version, $cmp_version) = @_;
  my $sql = get_sql(
    ['sr1.name as asm_name', 
    'sr1.length as asm_length', 
    'sr2.name as cmp_name', 
    'sr2.length as cmp_length', 
    'asm.asm_start', 'asm.asm_end', 'asm.cmp_start', 'asm.cmp_end'],
    'order by cs2.version, sr2.name, asm.cmp_start');
  return $dba->dbc->sql_helper->execute( -SQL => $sql, -PARAMS => [$asm_version, $cmp_version], -USE_HASHREFS => 1);
}

sub get_reverse_assembly_mappings {
  my ($dba, $asm_version, $cmp_version) = @_;
  # Reverse the columns to get the reverse mapping
  my $sql = get_sql(
  ['sr1.name as cmp_name', 
  'sr1.length as cmp_length', 
  'sr2.name as asm_name', 
  'sr2.length as asm_length', 
  'asm.cmp_start as asm_start', 'asm.cmp_end as asm_end', 'asm.asm_start as cmp_start', 'asm.asm_end as cmp_end'],
  'order by cs2.version, sr2.name, asm.cmp_start');
  return $dba->dbc->sql_helper->execute( -SQL => $sql, -PARAMS => [$asm_version, $cmp_version], -USE_HASHREFS => 1);
}

sub get_sql {
  my ($columns, $order_by) = @_;
  my $select = join(q{, }, @{$columns}, 'asm.ori', '(asm.asm_end - asm.asm_start)+1 as length');
  return <<SQL;
select ${select}
from coord_system cs1
join seq_region sr1 on (cs1.coord_system_id = sr1.coord_system_id)
join assembly asm on (sr1.seq_region_id = asm.asm_seq_region_id)
join seq_region sr2 on (sr2.seq_region_id = asm.cmp_seq_region_id)
join coord_system cs2 on (sr2.coord_system_id = cs2.coord_system_id)
where cs1.version =?
and cs2.version =?
$order_by
SQL
}

sub print_chains {
  my ($fh, $chains) = @_;
  while(my $chain_def = shift @{$chains}) {
    my $header = $chain_def->{header};
    my $gaps = $chain_def->{gaps};
    print $fh join(q{ }, @{$header});
    print $fh "\n";
    foreach my $gap (@{$gaps}) {
      print $fh join(q{ }, @{$gap});
      print $fh "\n";
    }
    print $fh "\n";
  }
}

sub create_ucsc_chains {
  my ($dba, $prod_name, $chains) = @_;
  my @new_chains;
  foreach my $chain_def (@{$chains}) {
    my $ensembl_name = $chain_def->{header}->[2];
    my $target_name = ensembl_to_ucsc_name($dba, $prod_name, $ensembl_name);
    next if $ensembl_name eq $target_name;
    my $new_chain_def = decode_json(encode_json($chain_def)); # quick clone
    # Substitute tName which in a GRCh37 -> GRCh38 mapping tName is GRCh37's code
    $new_chain_def->{header}[2] = $target_name;
    push(@new_chains, $new_chain_def);
  }
  return \@new_chains;
}

sub ensembl_to_ucsc_name {
  my ($dba, $prod_name, $ensembl_name) = @_;

  return $ucsc_name_cache{$prod_name}{$ensembl_name} if exists $ucsc_name_cache{$prod_name}{$ensembl_name};

  # Default is Ensembl name
  my $ucsc_name = $ensembl_name;
  my $slice = $dba->get_SliceAdaptor()->fetch_by_region(undef, $ensembl_name);
  my $synonyms = $slice->get_all_synonyms('UCSC');
  if(@{$synonyms}) {
    $ucsc_name = $synonyms->[0]->name();
  }
  else {
    if($slice->is_chromosome()) {
      #MT is a special case; it's chrM
      if($ensembl_name eq 'MT' ) {
        $ucsc_name = 'chrM';
      }
      # If it was a ref region add chr onto it (only check if we have an adaptor)
      elsif($slice->is_reference()) {
        $ucsc_name = 'chr'.$ensembl_name;
      }
    }
  }
  return $ucsc_name_cache{$prod_name}{$ensembl_name} = $ucsc_name;  
}
