use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::RattusNorvegicus;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("RattusNorvegicus Specific: loading assembly data");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["chromosome" , $ass_def, "default_version", 1               ],
     ["supercontig", undef   , "default_version", 2               ],
     ["contig"     , undef   , "default_version,sequence_level", 3]);

  my @assembly_mappings =  ("chromosome:$ass_def|contig",
                            "chromosome:$ass_def|supercontig",
                            "chromosome:$ass_def|contig|supercontig");

  my %cs = (gene                  => 'chromosome',
            transcript            => 'chromosome',
            exon                  => 'chromosome',
            dna_align_feature     => 'contig',
            protein_align_feature => 'contig',
            marker_feature        => 'contig',
            simple_feature        => 'contig',
            repeat_feature        => 'contig',
            qtl_feature           => 'chromosome',
            misc_feature          => 'chromosome',
            prediction_transcript => 'contig',
            karyotype             => 'chromosome');

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib, rank) VALUES (?,?,?,?)");

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

  $self->debug("RattusNorvegicus Specific: creating contig, " .
               "chromosome and supercontig seq_regions");

  $self->contig_to_seq_region();
  $self->chromosome_to_seq_region();
  $self->supercontig_to_seq_region();
}

sub create_assembly {
  my $self = shift;

  $self->debug("RattusNorvegicus Specific: loading assembly data");

  $self->assembly_contig_chromosome();
  $self->assembly_supercontig_chromosome();
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

  $self->debug("RattusNorvegicus Specific: Transforming contigs " .
               "into $target_cs_name seq_regions");

  my $cs_id = $self->get_coord_system_id($target_cs_name);

  my $sth = $dbh->prepare
    ("INSERT INTO $target.seq_region " .
     "SELECT ctg.contig_id, cln.embl_acc, " .
     "       $cs_id, ctg.length " .
     "FROM   $source.contig ctg, $source.clone cln " .
     "WHERE  ctg.clone_id = cln.clone_id");

  $sth->execute();
  $sth->finish();

  return;
}



1;
