use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::DanioRerio;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("DanioRerio Specific: creating scaffold coord system");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["scaffold" , $ass_def, "top_level,default_version,sequence_level"]);

  my %cs = (gene                  => 'scaffold',
            transcript            => 'scaffold',
            exon                  => 'scaffold',
            dna_align_feature     => 'scaffold',
            protein_align_feature => 'scaffold',
            marker_feature        => 'scaffold',
            simple_feature        => 'scaffold',
            repeat_feature        => 'scaffold',
            qtl_feature           => 'scaffold',
            misc_feature          => 'scaffold',
            prediction_transcript => 'scaffold',
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

  return;
}



#
# getter/setter needed by danio rerio to perform coord transforms
# of features after conversion of fake contigs in scaffolds
#
sub contig_mappings {
  my $self = shift;
  $self->{'contig_mappings'} = shift if(@_);
  return $self->{'contig_mappings'};
}


sub dna_mappings {
  my $self = shift;
  $self->{'dna_mappings'} = shift if(@_);
  return $self->{'dna_mappings'};
}

sub create_seq_regions {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->debug("DanioRerio Specific: creating scaffold seq_regions");

  #
  # retrieve fake contigs and group them into scaffolds
  #

  my $sth = $dbh->prepare
    ("SELECT c.contig_id, c.name, " .
     "a.contig_end - a.contig_start+1, dna_id "  .
     "FROM $source.contig c, $source.assembly a " .
     "WHERE c.contig_id = a.contig_id");
  $sth->execute();

  my %contig_lists;

  my($contig_id, $name, $length, $dna_id);
  $sth->bind_columns(\$contig_id, \$name, \$length, \$dna_id);

  while($sth->fetch()) {
    my $position;
    ($name, $position) = split(/\./, $name);
    $contig_lists{$name} ||= [];
    $contig_lists{$name}->[$position-1] = [$contig_id,$length,$dna_id];
  }
  
  $sth->finish();

  #
  # calculate the offsets of the fake contigs, and the lengths of the
  # scaffolds
  #

  my %old_new_mappings;
  my %dna_mappings;

  my %scaffolds;
  foreach my $ctgname (keys %contig_lists) {
    my $ctg_list = $contig_lists{$ctgname};
    $length = 0;
    my $new_id  = $ctg_list->[0]->[0];
    $dna_mappings{$new_id} = [];
    foreach my $contig (@$ctg_list) {
      my $old_id = $contig->[0];
      $old_new_mappings{$old_id} = [$new_id,$length+1];
      $length += $contig->[1];
      push(@{$dna_mappings{$new_id}}, [$contig->[2], $contig->[1]]); 
    }

    $scaffolds{$ctgname} = [$new_id, $length];
  }

  
  #
  # load the scaffold seq regions
  #
  $self->debug("DanioRerio Specific: loading scaffold seq regions");

  my $cs_id = $self->get_coord_system_id('scaffold');

  $sth = $dbh->prepare("INSERT INTO $target.seq_region (seq_region_id, " .
                          "name, length, coord_system_id) " .
                          "VALUES (?,?,?,?)");

  foreach my $scaf_name (keys %scaffolds) {
    my ($seq_reg_id, $len) = @{$scaffolds{$scaf_name}};
    $sth->execute($seq_reg_id, $scaf_name, $len, $cs_id);
  }
  $sth->finish();


  #
  # load temporary mapping of chromosomes to new seq_region ids.  Even though
  # there are not actually any chromosomes this info is used when transfering
  # genes to scaffold coords.
  #
  my $chr_select_sth = $dbh->prepare
    ("SELECT chromosome_id, name FROM $source.chromosome");
  $chr_select_sth->execute();

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_chr_map (old_id, new_id) VALUES (?,?)");

  while(my $row = $chr_select_sth->fetchrow_arrayref) {
    my ($chr_id, $chr_name) = @$row;
    $tmp_insert_sth->execute($chr_id, $scaffolds{$chr_name}->[0]);
  }
  
  #keep this mapping info for later processing of features & dna
  $self->dna_mappings(\%dna_mappings);
  $self->contig_mappings(\%old_new_mappings);

}



sub transfer_dna {
  my $self = shift;
  
  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->debug("DanioRerio Specific: building dna table");

  #in danio rerio est/estgene databases there is no dna so we may want to 
  #skip this function
  
  my $count_sth = $dbh->prepare("SELECT count(*) from $source.dna");
  $count_sth->execute();
  my ($count) = $count_sth->fetchrow_array();
  return if(!$count);

  my $dna_mappings = $self->dna_mappings();

  my $get_seq_sth = $dbh->prepare
    ("SELECT substring(sequence,1,?) from $source.dna " .
     "WHERE dna_id = ?");

  my $store_seq_sth = $dbh->prepare
    ("INSERT INTO $target.dna (seq_region_id, sequence) " .
     "VALUES (?,?)");

  my $update_seq_sth = $dbh->prepare
    ("UPDATE $target.dna " .
     "SET sequence = concat(sequence, ?) " .
     "WHERE seq_region_id = ?"); 

  foreach my $new_id (keys %$dna_mappings) {
    my $first = 1;
    foreach my $row (@{$dna_mappings->{$new_id}}) {
      my ($old_id, $length) = @$row;
      $get_seq_sth->execute($length, $old_id);

      if(!$get_seq_sth->rows() == 1) {
        warn("Missing dna sequence for dna id $old_id");
        next;
      }
      my $seq = $get_seq_sth->fetchrow_arrayref->[0];

      if($first) {
        #$self->debug("storing sequence of length " . length($seq) . 
        #             "for seq_region $new_id");
        $store_seq_sth->execute($new_id, $seq);
        $first = 0;
      } else {
        $update_seq_sth->execute($seq, $new_id);
        #$self->debug("concatting sequence to length " . length($seq));
      }
    }
  }

  $get_seq_sth->finish();
  $store_seq_sth->finish();
  $update_seq_sth->finish();
  return;
}


#
# override the transfer features so that features which are on fake contigs
# that were merged into scaffolds can have their coordinates adjusted
#

sub transfer_features {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->SUPER::transfer_features();

  my $contig_mappings = $self->contig_mappings();

  my @tables = qw(simple_feature
                  repeat_feature
                  protein_align_feature
                  dna_align_feature
                  marker_feature
                  prediction_transcript);
                  
  foreach my $table (@tables) {

    $self->debug("DanioRerio Specific: updating $table coordinates");

    my $sth = $dbh->prepare
      ("UPDATE $target.$table " .
       "SET    seq_region_id    = ?, " .
       "       seq_region_start = seq_region_start + ?, " .
       "       seq_region_end   = seq_region_end   + ? " . 
       "WHERE  seq_region_id = ?");

    foreach my $old_id (keys %$contig_mappings) {
      my ($new_id, $offset) = @{$contig_mappings->{$old_id}};
      if($offset > 1) {
        $sth->execute($new_id, $offset, $offset, $old_id)
      }
    }
  }
}



sub create_assembly {
  my $self = shift;

  $self->debug("DanioRerio Specific: no assembly data");

  return;
}



#
# override the contig_to_seqregion method so that contigs are given clone
# names instead
#
sub contig_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh     = $self->dbh();

  $target_cs_name ||= 'contig';

  $self->debug("DanioRerio Specific: Transforming contigs " .
               "into $target_cs_name seq_regions");

  my $cs_id = $self->get_coord_system_id($target_cs_name);

  my $sth = $dbh->prepare
    ("INSERT INTO $target.seq_region " .
     "SELECT ctg.contig_id, CONCAT(cln.embl_acc, '.', cln.embl_version), " .
     "       $cs_id, ctg.length " .
     "FROM   $source.contig ctg, $source.clone cln " .
     "WHERE  ctg.clone_id = cln.clone_id");

  $sth->execute();
  $sth->finish();

  return;
}



1;
