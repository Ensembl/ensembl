use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::CaenorhabditisBriggsae;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("CaenorhabditisBriggsae Specific: creating clone, scaffold," .
              " and contig coordinate systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["scaffold" , $ass_def, "top_level,default_version"     ],
     ['clone'      , undef   , 'default_version'               ],
     ["contig"     , undef   , "default_version,sequence_level"]);

  my @assembly_mappings =  ("scaffold:$ass_def|contig",
                            "clone|contig");

  my %cs = (gene                  => 'scaffold',
            transcript            => 'scaffold',
            exon                  => 'scaffold',
            dna_align_feature     => 'contig',
            protein_align_feature => 'contig',
            marker_feature        => 'contig',
            simple_feature        => 'contig',
            repeat_feature        => 'contig',
            qtl_feature           => 'scaffold',
            misc_feature          => 'scaffold',
            prediction_transcript => 'contig',
            karyotype             => 'scaffold');

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib) VALUES (?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
  }
  $sth->finish();

  $self->debug("Building meta_coord table");
  $sth = $dbh->prepare("INSERT INTO $target.meta_coord VALUES (?, ?)");
  foreach my $val (keys %cs) {
    $sth->execute($val, $coord_system_ids{$cs{$val}});
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



sub create_seq_regions {
  my $self = shift;

  $self->debug("CaenorhabditisBriggsae Specific: creating contig, " .
               "clone, contig and scaffold seq_regions");

  $self->contig_to_seq_region();
  $self->clone_to_seq_region();
  $self->chromosome_to_seq_region('scaffold');
}


sub chromosome_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $target_cs_name ||= "chromosome";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("CaenorhabditisBriggsae Specific: Transforming " .
               "chromosomes into $target_cs_name seq_regions");


  ## For consistancy with mart and v19 we need to keep chr name the same for
  ## now, so the following section is commented out and replaced:
  ##strip off the leading 'cb25.' from the chromosome name
  #my $select_sth = $dbh->prepare
  #  ("SELECT chromosome_id,substring(name,6),length FROM $source.chromosome");

  my $select_sth = $dbh->prepare
    ("SELECT chromosome_id,name,length FROM $source.chromosome");


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


sub create_assembly {
  my $self = shift;

  $self->debug("CaenorhabditisBriggsae Specific: loading assembly data");

  $self->assembly_contig_chromosome();
  $self->assembly_contig_clone();
}

#
# Override contig_to_seq_region and clone_to_seq_region to provide 
# briggsae specific behaviour
#

# sub contig_to_seq_region {
#   my $self = shift;
#   my $target_cs_name = shift;

#   my $target = $self->target();
#   my $source = $self->source();
#   my $dbh     = $self->dbh();

#   $target_cs_name ||= 'contig';

#   $self->debug("CaenorhabditisBriggsae Specific: Transforming contigs into " .
#                "$target_cs_name seq_regions");

#   my $cs_id = $self->get_coord_system_id($target_cs_name);

#   #There are two types of contigs in briggsae:

#   #
#   # cosmids/clones
#   #
#   my $sth = $dbh->prepare
#     ("INSERT INTO $target.seq_region " .
#      "SELECT contig_id, name, $cs_id, length " .
#      "FROM $source.contig " .
#      "WHERE  name not like 'c%'");

#   $sth->execute();
#   $sth->finish();

#   #
#   # WGS contigs
#   #
#   $sth = $dbh->prepare
#     ("INSERT INTO $target.seq_region " .
#      "SELECT ctg.contig_id, cln.name, $cs_id, length " .
#      "FROM   $source.contig ctg, $source.clone cln " .
#      "WHERE  ctg.clone_id = cln.clone_id " .
#      "AND    ctg.name like 'c%'");

#   $sth->execute();
#   $sth->finish();

#   return;
# }



sub clone_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  # target coord_system will have a different ID
  $target_cs_name ||= "clone";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("CaenorhabditisBriggsae Specific:Transforming clones " .
               "into $target_cs_name seq_regions");

  #
  # We don't want to make clones out of the WGS contigs, only out of
  # the clones with proper embl accessions
  #
  my $select_sth = $dbh->prepare
    ("SELECT cl.clone_id,
             CONCAT(cl.embl_acc, '.', cl.embl_version),
             MAX(ctg.embl_offset)+ctg.length-1
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id
     AND    cl.embl_acc not like 'c%'
     GROUP BY ctg.clone_id");
 
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


1;
