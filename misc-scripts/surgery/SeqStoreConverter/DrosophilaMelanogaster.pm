use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::DrosophilaMelanogaster;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("DrosophilaMelanogaster Specific: creating " .
               "chunk and chromosome coord systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["chromosome" , $ass_def, "default_version", 1     ],
     ["chunk",      undef   , "default_version,sequence_level", 2]);

  my @assembly_mappings =  ("chromosome:$ass_def|chunk");

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib, rank) VALUES (?,?,?,?)");

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



sub create_seq_regions {
  my $self = shift;

  $self->debug("DrosophilaMelanogaster Specific: creating chunk and " .
               "chromosome seq_regions");

  $self->contig_to_seq_region('chunk');
  $self->chromosome_to_seq_region();
}


sub create_assembly {
  my $self = shift;

  $self->debug("DrosophilaMelanogaster Specific: loading assembly data");

  $self->assembly_contig_chromosome();
}


1;
