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

# Convert release 18/19-era schemas to use new non-clone/contig schema

use strict;
use warnings;
use DBI;


package SeqStoreConverter::BasicConverter;

###############################################################################
# Constructor
###############################################################################

sub new {
  my ( $class, $user, $pass, $host, $source, $target, $schema, $vega_schema, $force, $verbose, $limit ) = @_;

  my $self = bless {}, $class;

  my $port;
  ($host, $port) = split(/:/, $host);
  $port ||= 3306;

  my $dbh = DBI->connect( "DBI:mysql:host=$host:port=$port", $user, $pass,
                          {'RaiseError' => 1});


  $self->verbose( $verbose );
  $self->dbh( $dbh );
  $self->force( $force );
  $self->source( $source );
  $self->target( $target );
  $self->schema( $schema );
  $self->vegaschema( $vega_schema);
  $self->host( $host );
  $self->password( $pass);
  $self->user($user);
  $self->port($port);
  $self->limit($limit);


  #check to see if the destination and source databases exist already.
  my %dbs = map {$_->[0] => 1} @{$dbh->selectall_arrayref('show databases')};
  if( !$dbs{$source} ) {
    die ("Source db $source does not exist" );
  }

  if( $dbs{$target} ) {
    if( $force ) {
      $dbh->do( "drop database $target" );
    } else {
      die("Target db $target already exists. Use -force option to overwrite.");
    }
  }

  $dbh->do( "create database ".$self->target() );

  $self->debug("Building schema for $target from $schema");
  die "Cannot open $schema" if (! -e $schema);
  my $cmd = "/usr/local/mysql/bin/mysql -u $user -p$pass -P $port -h $host $target < $schema";
  system ($cmd);

  if ($vega_schema) {
      $self->debug("Adding vega tables for $target");
      die "Cannot open vega creation script" if (! -e $vega_schema);
      my $cmd = "/usr/local/mysql/bin/mysql -u $user -p$pass -P $port -h $host $target < $vega_schema";
      system ($cmd);
  }

  $self->debug("Creating temporary tables");
  #create a temporary table to store the mapping of old ids to new ids
  $dbh->do
    ("CREATE TEMPORARY TABLE $target.tmp_cln_map (" .
     "old_id INT, new_id INT, INDEX new_idx (new_id))");

  #create a temp table which will store the mapping of old chromosome
  #identifiers to new identifiers
  $dbh->do("CREATE TEMPORARY TABLE $target.tmp_chr_map (" .
           "  old_id INT, new_id INT,".
           "  INDEX new_idx (new_id))");

  #create a temporary table to hold old supercontig name -> new id mappings
  $dbh->do("CREATE TEMPORARY TABLE $target.tmp_superctg_map (" .
           "name VARCHAR(255), new_id INT, ".
           "INDEX new_idx (new_id))");

  return $self;
}


###############################################################################
# Getter/Setters for converter properties
###############################################################################

sub force {
  my $self = shift;
  $self->{'force'} = shift if(@_);
  return $self->{'force'};
}


sub dbh {
  my $self = shift;
  $self->{'dbh'} = shift if(@_);
  return $self->{'dbh'};
}

sub user {
  my $self = shift;
  $self->{'user'} = shift if(@_);
  return $self->{'user'};
}
sub host {
  my $self = shift;
  $self->{'host'} = shift if(@_);
  return $self->{'host'};
}

sub port {
  my $self = shift;
  $self->{'port'} = shift if(@_);
  return $self->{'port'};
}

sub password {
  my $self = shift;
  $self->{'password'} = shift if(@_);
  return $self->{'password'};
}


sub verbose {
  my $self = shift;
  $self->{'verbose'} = shift if (@_);
  return $self->{'verbose'};
}

sub schema {
  my $self = shift;
  $self->{'schema'} = shift if (@_);
  return $self->{'schema'};
}

sub vegaschema {
  my $self = shift;
  $self->{'vega_schema'} = shift if (@_);
  return $self->{'vega_schema'};
}

sub source {
  my $self = shift;
  $self->{'source'} = shift if(@_);
  return $self->{'source'};
}


sub target {
  my $self = shift;
  $self->{'target'} = shift if(@_);
  return $self->{'target'};
}

sub limit {
  my $self = shift;
  $self->{'limit'} = shift if(@_);
  return $self->{'limit'};
}


###############################################################################
# Utility methods
###############################################################################

sub debug {
  my $self = shift;
  my $str = shift;
  print STDERR $str . "\n" if $self->verbose();
  return;
}

sub copy_tables {
  my ($self, @tables) = @_;

  foreach my $table (@tables) {
    $self->debug("Copying $table");
    
    my $source = $self->source();
    my $target = $self->target();

    eval {
      my $sth = $self->dbh()->prepare
        ("INSERT INTO $target.$table SELECT * FROM $source.$table");
      $sth->execute();
      $sth->finish();
    };

    if($@) {
      warn("Copy of table $table failed: $@\n");
    }
  }

  return;
}

sub get_coord_system_id {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift;

  my $target = $self->target();

  my @bind_vals = ($cs_name);
  my $sql = "SELECT cs.coord_system_id " .
            "FROM   $target.coord_system cs " .
            "WHERE  cs.name = ?";

  if($cs_version) {
    push(@bind_vals, $cs_version);
    $sql .= " AND cs.version = ?";
  }

  my $sth = $self->dbh()->prepare($sql);
  $sth->execute(@bind_vals);

  if($sth->rows() != 1) {
    die("Id for non-existant or ambiguous coord system requested " .
        "$cs_name:$cs_version");
  }

  my ($id) = $sth->fetchrow_array();

  $sth->finish();

  return $id;
}

sub get_default_assembly {
  my $self = shift;

  my $source = $self->source();

  my $sth = $self->dbh->prepare
    ("SELECT meta_value FROM $source.meta WHERE meta_key='assembly.default'");
  $sth->execute();

  if(!$sth->rows() == 1) {
    die("This species has an ambiguous or non-existant assembly.default" .
        " in the meta table");
  }

  my ($result) = $sth->fetchrow_array();
  
  return $result;
}



sub contig_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh     = $self->dbh();

  $target_cs_name ||= 'contig';

  $self->debug("Transforming contigs into $target_cs_name seq_regions");

  my $cs_id = $self->get_coord_system_id($target_cs_name);

  my $sth = $dbh->prepare
    ("INSERT INTO $target.seq_region " .
     "SELECT contig_id, name, $cs_id, length FROM $source.contig");

  $sth->execute();
  $sth->finish();

  return;
}


sub clone_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  # target coord_system will have a different ID
  $target_cs_name ||= "clone";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("Transforming clones into $target_cs_name seq_regions");

  my $select_sth = $dbh->prepare
    ("SELECT cl.clone_id,
             CONCAT(cl.embl_acc, '.', cl.embl_version),
             MAX(ctg.embl_offset+ctg.length-1)
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id GROUP BY ctg.clone_id");
  $select_sth->execute();

  my ($clone_id, $embl_acc, $length);
  $select_sth->bind_columns(\$clone_id, \$embl_acc, \$length);

  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES(?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_cln_map (old_id, new_id) VALUES (?, ?)");

  while ($select_sth->fetch()) {
    $insert_sth->execute("$embl_acc", $cs_id, $length);

    #store mapping of old -> new ids in temp table
    $tmp_insert_sth->execute($clone_id, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}



sub chromosome_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $target_cs_name ||= "chromosome";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("Transforming chromosomes into $target_cs_name seq_regions");


  my $select_sth = $dbh->prepare
    ("SELECT chromosome_id, name, length FROM $source.chromosome");

  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES (?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_chr_map (old_id, new_id) VALUES (?, ?)");

  $select_sth->execute();

  my ($chrom_id, $name, $length);
  $select_sth->bind_columns(\$chrom_id, \$name, \$length);

  while ($select_sth->fetch()) {
    #insert into seq_region table
    $insert_sth->execute($name, $cs_id, $length);
    #copy old/new mapping into temporary table
    $tmp_insert_sth->execute($chrom_id, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}



sub supercontig_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift || "supercontig";

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $self->debug("Transforming supercontigs into $target_cs_name seq_regions");

  my $cs_id = $self->get_coord_system_id($target_cs_name);

  my $select_sth = $dbh->prepare
    ("SELECT superctg_name, " .
     "MAX(superctg_end) AS length " .
     "FROM $source.assembly " .
     "GROUP BY superctg_name");

  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES (?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_superctg_map (name, new_id) VALUES (?, ?)");

  my ($name, $length);
  $select_sth->execute();
  $select_sth->bind_columns(\$name, \$length);

  while ($select_sth->fetch()) {
    $insert_sth->execute($name, $cs_id, $length);
    $tmp_insert_sth->execute($name, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}


sub assembly_contig_chromosome {
  my $self = shift;

  $self->debug("Building assembly table - contig/chromosome");

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $dbh->do(
	  "INSERT INTO $target.assembly " .
	  "SELECT tcm.new_id, " .	# asm_seq_region_id (old-new chromosome ID mapping)
	  "a.contig_id, " .	# cmp_seq_region_id
	  "a.chr_start, " .	# asm_start
	  "a.chr_end, " .		# asm_end
	  "a.contig_start, " .	# cmp_start
	  "a.contig_end, " .	# cmp_end
	  "a.contig_ori " .	# ori
	  "FROM $target.tmp_chr_map tcm, $source.assembly a, $source.contig c " .
	  "WHERE tcm.old_id = a.chromosome_id " .
	  "AND c.contig_id = a.contig_id "); # only copy assembly entries that 
                                       # refer to valid contigs (test db has
                                       # superfluous assembly entries) 

}

sub assembly_contig_clone {
  my $self = shift;
  
  $self->debug("Building assembly table - contig/clone");

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $dbh->do(
	  "INSERT INTO $target.assembly " .
	  "SELECT tcm.new_id, " .	# asm_seq_region_id (old-new clone ID mapping)
	  "ctg.contig_id, ".	# cmp_seq_region_id
	  "ctg.embl_offset, " .	# asm_start
	  "ctg.embl_offset+ctg.length-1, " . # asm_end
	  "1, " .			# cmp_start
	  "ctg.length, " .	# cmp_end
	  "1 " .	# ori - contig always positively oriented on the clone
	  "FROM $target.tmp_cln_map tcm, " .
	  "$source.clone cln, $source.contig ctg " .
	  "WHERE tcm.old_id = cln.clone_id " .
	  "AND cln.clone_id = ctg.clone_id");

}


sub assembly_contig_supercontig {
  my $self = shift;

  $self->debug("Building assembly table - contig/supercontig");

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $dbh->do(
	  "INSERT INTO $target.assembly " .
	  "SELECT tsm.new_id, " .	# asm_seq_region_id (superctg name-sr_id mapping)
	  "a.contig_id, " .	# cmp_seq_region_id
	  "a.superctg_start, " .	# asm_start
	  "a.superctg_end, " .	# asm_end
	  "a.contig_start, " .	# cmp_start
	  "a.contig_end, " .	# cmp_end
	  "a.contig_ori " .	# ori
	  "FROM $target.tmp_superctg_map tsm, $source.assembly a, $source.contig c ".
	  "WHERE tsm.name = a.superctg_name " .
	  "AND c.contig_id = a.contig_id "); # only copy assembly entries that 
                                       # refer to valid contigs (test db might
                                       # have these)

}



sub assembly_supercontig_chromosome {
  my $self = shift;

  $self->debug("Building assembly table - supercontig/chromosome");

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $dbh->do(
	  "INSERT INTO $target.assembly " .
	  "SELECT tcm.new_id, " .	# asm_seq_region_id (chr id)
	  "tsm.new_id, " .	# cmp_seq_region_id (supercontig id)
	  "min(a.chr_start), " .	# asm_start
	  "max(a.chr_end), " .	# asm_end
	  "min(a.superctg_start), " .	# cmp_start
	  "max(a.superctg_end), " .	# cmp_end
	  "a.superctg_ori " .	# ori
	  "FROM $target.tmp_superctg_map tsm, $target.tmp_chr_map tcm, " .
    "     $source.assembly a ".
	  "WHERE tsm.name = a.superctg_name " .
    "AND   tcm.old_id = a.chromosome_id " .
    "GROUP BY superctg_name");

}



###############################################################################
# Base class implementations of transfer methods. Can be overridden to
# create species specific behaviour
###############################################################################


sub create_coord_systems {
  my $self = shift;

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["chromosome" , $ass_def, "default_version"               ,1],
     ["supercontig", undef   , "default_version"               ,2],
     ["clone"      , undef   , "default_version"               ,3],
     ["contig", undef        , "default_version,sequence_level",4]);

  my @assembly_mappings =  ("chromosome:$ass_def|contig",
                            "clone|contig",
                            "supercontig|contig",
                            "supercontig|contig|clone",
                            "chromosome:$ass_def|contig|clone",
                            "chromosome:$ass_def|contig|supercontig");

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib,rank) VALUES (?,?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
  }
  $sth->finish();

  $self->debug("Adding assembly.mapping entries to meta table");

  $sth = $dbh->prepare("INSERT INTO $target.meta(meta_key, meta_value) " .
                       "VALUES ('assembly.mapping', ?)");

  foreach my $mapping (@assembly_mappings) {
    $sth->execute($mapping);
  }

  $sth->finish();

  return;
}


#
# populates the contents of the meta_coord table
# must be executed after all of the feature tables in the target database
# have already been populated
#

sub create_meta_coord {
  my $self = shift;

  $self->debug("Building meta_coord table");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my @feature_tables = qw(density_feature
                          dna_align_feature
                          exon
                          gene
                          karyotype
                          marker_feature
                          misc_feature
                          prediction_exon
                          prediction_transcript
                          protein_align_feature
                          repeat_feature
                          simple_feature
                          transcript);

  foreach my $ft (@feature_tables) {

    $dbh->do(qq{INSERT INTO $target.meta_coord(table_name, coord_system_id,
                                               max_length)
                SELECT '$ft', sr.coord_system_id,
                        MAX(f.seq_region_end - f.seq_region_start + 1)
                FROM $target.$ft f, $target.seq_region sr
                WHERE sr.seq_region_id = f.seq_region_id
                GROUP BY sr.coord_system_id});
  }

  # special case for assembly exception, features are created from both
  # sides of table

  $dbh->do(qq{INSERT INTO $target.meta_coord
              SELECT 'assembly_exception', sr.coord_system_id, 
              MAX(IF(ae.seq_region_end - ae.seq_region_start > ae.exc_seq_region_end - ae.exc_seq_region_start, ae.seq_region_end - ae.seq_region_start + 1, ae.exc_seq_region_end - ae.exc_seq_region_start + 1))
              FROM   $target.assembly_exception ae, $target.seq_region sr
              WHERE  sr.seq_region_id = ae.seq_region_id
              GROUP BY sr.coord_system_id});

}


sub create_seq_regions {
  my $self = shift;

  my $target = $self->target();
  my $dbh    = $self->dbh();

  #default behaviour is to simply copy all tables as they come

  $self->contig_to_seq_region('contig');
  $self->chromosome_to_seq_region();
  $self->supercontig_to_seq_region();
  $self->clone_to_seq_region();
  
  return;
}

sub create_assembly {
  my $self = shift;

  $self->assembly_contig_chromosome();
  $self->assembly_contig_clone();
  $self->assembly_contig_supercontig();

  return;
}

sub transfer_dna {
  my $self = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $self->debug("Building dna table");

  $dbh->do("INSERT INTO $target.dna " .
          "SELECT c.contig_id as seq_region_id, d.sequence as sequence " .
          "FROM   $source.dna d, $source.contig c " .
          "WHERE  c.dna_id = d.dna_id");
  return;
}


sub transfer_genes {
  my $self = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  #
  # Transfer the gene table
  #

  $self->debug("Building gene table");

  $dbh->do
    ("INSERT INTO $target.gene " .
     "SELECT g.gene_id, g.type, g.analysis_id, tcm.new_id, " .
     "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
     "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
     "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
     "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
     "       a.contig_ori*e.contig_strand as strand, " .
     "       g.display_xref_id " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.assembly a, $source.gene g, " .
     "       $target.tmp_chr_map tcm " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    e.contig_id = a.contig_id " .
     "AND    g.gene_id = t.gene_id " .
     "AND    a.chromosome_id = tcm.old_id " .
     "GROUP BY g.gene_id");


  # 
  # Transfer the transcript table
  #
  $self->debug("Building transcript table");

  $dbh->do
    ("INSERT INTO $target.transcript " .
     "SELECT t.transcript_id, t.gene_id, tcm.new_id, " .
     "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
     "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
     "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
     "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
     "       a.contig_ori*e.contig_strand as strand, " .
     "       t.display_xref_id " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.assembly a, $target.tmp_chr_map tcm " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    e.contig_id = a.contig_id " .
     "AND    a.chromosome_id = tcm.old_id " .
     "GROUP BY t.transcript_id");

  #
  # Transfer the exon table
  #
  $self->debug("Building exon table");

  $dbh->do
    ("INSERT INTO $target.exon " .
     "SELECT e.exon_id, tcm.new_id, " .
     "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
     "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
     "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
     "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
     "       a.contig_ori*e.contig_strand as strand, " .
     "       e.phase, e.end_phase " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.assembly a, $source.gene g, " .
     "       $target.tmp_chr_map tcm " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    e.contig_id = a.contig_id " .
     "AND    g.gene_id = t.gene_id " .
     "AND    a.chromosome_id = tcm.old_id " .
     "GROUP BY e.exon_id");

  # 
  # Transfer translation table
  # 

  $self->debug("Building translation table");

  $dbh->do
    ("INSERT INTO $target.translation " .
     "SELECT tl.translation_id, ts.transcript_id, tl.seq_start, " .
     "       tl.start_exon_id, tl.seq_end, tl.end_exon_id " .
     "FROM $source.transcript ts, $source.translation tl " .
     "WHERE ts.translation_id = tl.translation_id");

  return;
}

sub transfer_prediction_transcripts {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();
  
  # prediction_transcript / prediction_exon
  
  $self->debug( "Building prediction_exon table" );

  $dbh->do
    ("INSERT INTO $target.prediction_exon ".
	   "( prediction_transcript_id, seq_region_id, seq_region_start, " .
     "  seq_region_end, seq_region_strand, start_phase, score, p_value," .
     "  exon_rank ) " .
	   "SELECT prediction_transcript_id, contig_id, contig_start, contig_end, " .
     "       contig_strand, start_phase, score, p_value, exon_rank " .
	   "FROM   $source.prediction_transcript" );

  $self->debug("Building prediction_transcript table");

  $dbh->do
    ("INSERT INTO $target.prediction_transcript ".
	   "( prediction_transcript_id, seq_region_id, seq_region_start, " .
     "seq_region_end,  seq_region_strand, analysis_id ) " .
	   "SELECT prediction_transcript_id, contig_id, MIN(contig_start), " .
     "       MAX(contig_end), contig_strand, analysis_id ".
	   "FROM   $source.prediction_transcript " .
	   "GROUP BY prediction_transcript_id ");

  return;
}



sub transfer_features {
  my $self = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  my $limit = '';
  if($self->limit()) {
    $limit = ' limit ' . $self->limit();
  }


  #
  # Feature tables
  # Note that we can just rename contig_* to set_region_* since the
  # contig IDs were copied verbatim into seq_region
  #

  # For some reason mysql occasionally refuses to use the index on large 
  # tables following copies like the following.
  # So: drop the indexes first and then re-add them after

  # simple_feature
  $self->debug("Building simple_feature table");
  $dbh->do
    ("INSERT INTO $target.simple_feature (simple_feature_id, seq_region_id, ".
     "       seq_region_start, seq_region_end, seq_region_strand, " .
     "       display_label, analysis_id, score) " .
     "SELECT simple_feature_id, contig_id, contig_start, contig_end, " .
     "       contig_strand, display_label, analysis_id, score " .
     "FROM $source.simple_feature $limit");
  
  # repeat_feature
  $self->debug("Dropping indexes on repeat_feature");
  $dbh->do("ALTER TABLE $target.repeat_feature DROP INDEX seq_region_idx");
  $dbh->do("ALTER TABLE $target.repeat_feature DROP INDEX repeat_idx");
  $dbh->do("ALTER TABLE $target.repeat_feature DROP INDEX analysis_idx");

  $self->debug("Building repeat_feature table");
  $dbh->do
    ("INSERT INTO $target.repeat_feature (repeat_feature_id, seq_region_id, " .
     "    seq_region_start, seq_region_end, seq_region_strand, analysis_id, " .
     "    repeat_start, repeat_end, repeat_consensus_id, score) " .
     "SELECT repeat_feature_id, contig_id, contig_start, contig_end, " .
     "       contig_strand, analysis_id, repeat_start, repeat_end, " .
     "       repeat_consensus_id, score FROM $source.repeat_feature $limit");

  $self->debug("Reading indexes on repeat_feature");
  $dbh->do("ALTER TABLE $target.repeat_feature " .
           "ADD INDEX seq_region_idx( seq_region_id, seq_region_start)");
  $dbh->do("ALTER TABLE $target.repeat_feature " .
           "ADD INDEX repeat_idx( repeat_consensus_id )");
  $dbh->do("ALTER TABLE $target.repeat_feature " .
           "ADD INDEX analysis_idx(analysis_id)");

  # protein_align_feature
  $self->debug("Dropping indexes on protein_align_feature");
  $dbh->do("ALTER TABLE $target.protein_align_feature DROP INDEX hit_idx");
  $dbh->do( "ALTER TABLE $target.protein_align_feature " .
            "DROP INDEX seq_region_idx");

  $self->debug("Building protein_align_feature table");
  $dbh->do( "INSERT INTO $target.protein_align_feature " .
            "  (protein_align_feature_id, seq_region_id, seq_region_start, " .
            "   seq_region_end, seq_region_strand, analysis_id, hit_start, " .
            "   hit_end, hit_name, cigar_line, evalue, perc_ident, score) " .
            "SELECT protein_align_feature_id, contig_id, contig_start, " .
            "       contig_end, contig_strand, analysis_id, hit_start, " .
            "       hit_end, hit_name, cigar_line, evalue, perc_ident, score ".
            "FROM $source.protein_align_feature $limit");

  $self->debug("Reading indexes on protein_align_feature");
  $dbh->do( qq{ALTER TABLE $target.protein_align_feature
               ADD index  seq_region_idx(  analysis_id, seq_region_id,
                                           seq_region_start, score )});
  $dbh->do( "ALTER TABLE $target.protein_align_feature " .
            "ADD index hit_idx(hit_name)");

  # dna_align_feature
  $self->debug("Dropping indexes on dna_align_feature");
  $dbh->do( "ALTER TABLE $target.dna_align_feature DROP INDEX seq_region_idx");
  $dbh->do( "ALTER TABLE $target.dna_align_feature DROP INDEX hit_idx");

  $self->debug("Building dna_align_feature table");
  $dbh->do( "INSERT INTO $target.dna_align_feature " .
            "       (dna_align_feature_id, seq_region_id, seq_region_start, ".
            "        seq_region_end, seq_region_strand, analysis_id, " .
            "        hit_start, hit_end, hit_name, hit_strand, cigar_line, " .
            "        evalue, perc_ident, score) " .
            "SELECT dna_align_feature_id, contig_id, contig_start, " .
            "       contig_end, contig_strand, analysis_id, hit_start, " .
            "       hit_end, hit_name, hit_strand, cigar_line, evalue, " .
            "       perc_ident, score FROM $source.dna_align_feature $limit");


  $self->debug("Reading indexes on dna_align_feature");
  $dbh->do( qq{ALTER TABLE $target.dna_align_feature
               ADD INDEX seq_region_idx(seq_region_id, analysis_id,
                                        seq_region_start, score)});
  $dbh->do( "ALTER TABLE $target.dna_align_feature " .
            "ADD index hit_idx(hit_name)");

  # marker_feature
  $self->debug("Building marker_feature table");
  $dbh->do( "INSERT INTO $target.marker_feature " .
            "       (marker_feature_id, marker_id, seq_region_id, " .
            "        seq_region_start, seq_region_end, analysis_id, " .
            "        map_weight) " .
            "SELECT marker_feature_id, marker_id, contig_id, contig_start, ".
            "       contig_end, analysis_id, map_weight " .
            "FROM   $source.marker_feature $limit");
  
  # qtl_feature
  # Note this uses chromosome coords so we have to join with tmp_chr_map to 
  # get the mapping
  $self->debug("Building qtl_feature table");

  $dbh->do
    ("INSERT INTO $target.qtl_feature( seq_region_id, seq_region_start, " .
     "  seq_region_end, qtl_id, analysis_id) " .
     "SELECT tcm.new_id, " .
     "       q.start, q.end, q.qtl_id, q.analysis_id " .
     "FROM $target.tmp_chr_map tcm, $source.qtl_feature q " .
     "WHERE tcm.old_id = q.chromosome_id $limit");
  
  # These tables now have seq_region_* instead of chromosome_*
  
  $self->debug("Building karyotype table");
  $dbh->do(
           "INSERT INTO $target.karyotype " .
           "SELECT null, tcm.new_id, " .
           "       k.chr_start, k.chr_end, k.band, k.stain " .
           "FROM $target.tmp_chr_map tcm, $source.karyotype k " .
           "WHERE tcm.old_id = k.chromosome_id $limit");


  $self->debug("Building marker_map_location table");
  $dbh->do(
           "INSERT INTO $target.marker_map_location " .
           "SELECT mml.marker_id, mml.map_id, " .
           "       c.name, " .
           "       mml.marker_synonym_id, mml.position, mml.lod_score " .
           "FROM $source.chromosome c, $source.marker_map_location mml " .
           "WHERE c.chromosome_id = mml.chromosome_id $limit");

  $self->debug( "Building misc_feature table" );
  $dbh->do
    ("INSERT INTO $target.misc_feature( misc_feature_id, seq_region_id, " . 
     "            seq_region_start, seq_region_end, seq_region_strand ) " .
     "SELECT m.mapfrag_id, sr.seq_region_id, m.seq_start, m.seq_end, " .
     "       m.orientation " .
     "FROM   $source.mapfrag m, $target.seq_region sr, $source.dnafrag d " . 
     "WHERE  m.dnafrag_id = d.dnafrag_id " .
     "AND    d.name = sr.name $limit" );

  $self->debug( "Building misc_set table" );
  $dbh->do
    ("INSERT INTO $target.misc_set( misc_set_id, code, name, description, " .
     "                              max_length ) " .
     "SELECT mapset_id, code, name, description, max_length " . 
     "FROM $source.mapset ms" );

  $self->debug( "Building misc_attrib table" );
  $dbh->do
	 ("INSERT INTO $target.misc_attrib( misc_feature_id, attrib_type_id, " . 
    "                    value ) ". 
    "SELECT mapfrag_id, mapannotationtype_id, value " .
    "FROM $source.mapannotation" );

  $dbh->do
   ("INSERT INTO $target.misc_attrib( misc_feature_id, attrib_type_id, " .
    "                                 value) " .
    "SELECT mf.mapfrag_id, at.attrib_type_id, mf.name " .
    "FROM   $source.mapfrag mf, $target.attrib_type at " .
    "WHERE  at.code = 'name'");

   $dbh->do
   ("INSERT INTO $target.misc_attrib( misc_feature_id, attrib_type_id, " .
    "                                 value) " .
    "SELECT mf.mapfrag_id, at.attrib_type_id, mf.type " .
    "FROM   $source.mapfrag mf, $target.attrib_type at " .
    "WHERE  at.code = 'type'");

  $self->debug( "Building misc_feature_misc_set table" );
  $dbh->do
    ("INSERT INTO $target.misc_feature_misc_set(misc_feature_id, misc_set_id)".
     "SELECT mapfrag_id, mapset_id ".
     "FROM $source.mapfrag_mapset $limit" );

  return;
}


sub transfer_stable_ids {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->debug("Building stable id event tables");

  $self->copy_tables
    ("stable_id_event","mapping_session","gene_archive","peptide_archive");

  return;
}

sub transfer_vega_stable_ids {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->debug("Building vega_stable id event tables");


  $self->copy_tables
    ("stable_id_event","mapping_session","gene_archive","peptide_archive");

  return;
}

sub transfer_meta {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();

  my $dbh = $self->dbh();

  $dbh->do("INSERT INTO $target.meta (meta_key, meta_value) " .
           "SELECT m.meta_key, m.meta_value FROM $source.meta m " .
           "ORDER BY meta_id");

  return;
}


sub copy_other_tables {
  my $self = shift;

  #xref tables
  $self->copy_tables("xref",
                     "go_xref",
                     "identity_xref",
                     "object_xref",
                     "external_db",
                     "external_synonym",
  #marker/qtl related tables
                     "map",
                     "marker",
                     "marker_synonym",
                     "qtl",
                     "qtl_synonym",
  #misc other tables
                    "supporting_feature",
                    "analysis",
                    "exon_transcript",
                    "interpro",
                    "gene_description",
                    "protein_feature");
}


sub copy_repeat_consensus {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();

  my $dbh = $self->dbh();

  $self->debug("Converting repeat_consensus table.");

  $dbh->do("INSERT INTO $target.repeat_consensus " .
           "(repeat_consensus_id, repeat_name, repeat_class, repeat_type, ".
           " repeat_consensus) " .
           "SELECT repeat_consensus_id, repeat_name, repeat_class, " .
           "       '', repeat_consensus " .
           "FROM $source.repeat_consensus rc" );

  return;
}



sub create_attribs {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  #copy the attrib types from the mapannotation type table
  
  $dbh->do
    ("INSERT INTO $target.attrib_type( attrib_type_id, code, " .
     "                                 name, description ) " .
     "SELECT mapannotationtype_id, code, name, description " .
     "FROM $source.mapannotationtype " );

  $dbh->do
    ("INSERT INTO $target.attrib_type( code, name, description ) " .
     "VALUES ('name', 'Name',''), ('type', 'Type of feature','')");

  return;
}


#
# The process of actually identifying toplevel seq_regions using the info in
# the database is quite slow.  Make the assumption that the coordsystem with
# lowest rank value is going to have all of the toplevel seq_regions.
#
# This method must be overridden if alternate behaviour is required.
#

sub set_top_level {
  my $self = shift;

  my $target = $self->target();
  my $dbh = $self->dbh();

  my $attrib_type_id = $self->add_attrib_code();

  $self->debug("Setting toplevel attributes of seq_regions");

  my $sth = $dbh->prepare("DELETE FROM $target.seq_region_attrib " .
                          "WHERE attrib_type_id = ?");
  $sth->execute($attrib_type_id);
  $sth->finish();


  $sth = $dbh->prepare("SELECT coord_system_id FROM $target.coord_system " .
                       "ORDER BY RANK ASC LIMIT 1");
  $sth->execute();

  my ($cs_id) = $sth->fetchrow_array();

  $sth->finish();

  $sth = $dbh->prepare("INSERT INTO $target.seq_region_attrib " .
                      '(seq_region_id, attrib_type_id, value) ' .
                      "SELECT sr.seq_region_id, $attrib_type_id, 1 " .
                      "FROM $target.seq_region sr " .
                      "WHERE sr.coord_system_id = $cs_id");

  $sth->execute();
  $sth->finish();

}

sub add_attrib_code {
  my $self = shift;
  my $dbh = $self->dbh();
  my $target = $self->target();

  # add a toplevel code to the attrib_type table if it is not there already

  my $sth = $dbh->prepare("SELECT attrib_type_id " .
                         "FROM $target.attrib_type " .
                         "WHERE code = 'toplevel'");

  $sth->execute();

  if($sth->rows()) {
    my ($attrib_type_id) = $sth->fetchrow_array();
    $sth->finish();
    return $attrib_type_id;
  }
  $sth->finish();


  $sth = $dbh->prepare("INSERT INTO $target.attrib_type " .
                    "SET code = 'toplevel', " .
                    "name = 'Top Level', " .
                    "description = 'Top Level Non-Redundant Sequence Region'");

  $sth->execute();
  my $attrib_type_id = $sth->{'mysql_insertid'};
  $sth->finish();

  return $attrib_type_id;
}

1;
