use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::ApisMellifera;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("ApisMellifera Specific: creating " .
               "contig and scaffold coord systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["scaffold" , $ass_def, "default_version", 1     ],
     ["contig",      undef   , "default_version,sequence_level", 2]);

  my @assembly_mappings =  ("chromosome:$ass_def|chunk");

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
            prediction_exon       => 'contig',
            karyotype             => 'scaffold');

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

  $self->debug( "ApisMellifera Specific: creating contig and " .
               "scaffold seq_regions");

  $self->contig_to_seq_region('contig');
  $self->chromosome_to_seq_region("scaffold");
}


sub create_assembly {
  my $self = shift;

  $self->debug("DrosophilaMelanogaster Specific: loading assembly data");

  $self->assembly_contig_chromosome();
}


1;
