use strict;
use warnings;

use SeqStoreConverter::MusMusculus;

package SeqStoreConverter::vega::MusMusculus;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::MusMusculus);

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
		     "protein_feature",
  #vega tables
		     "gene_synonym",
		     "transcript_info",
		     "current_gene_info",
		     "current_transcript_info",
		     "author",
		     "gene_name",
		     "transcript_class",
		     "gene_remark",
		     "gene_info",
		     "evidence",
		     "transcript_remark",
		     "clone_remark",
		     "clone_info",
		     "clone_info_keyword",
		     "clone_lock",
		     "current_clone_info",
		     "keyword",
		     "job",
		     "job_status",
		     "input_id_analysis");
}

sub update_clone_info {
  my $self = shift;
  return;
}

sub remove_supercontigs {
    my $self = shift;
    
    my $target = $self->target();
    my $dbh    = $self->dbh();
    $self->debug("Vega mouse specific - removing supercontigs from $target");

    $dbh->do("DELETE FROM $target.meta ". 
	     "WHERE meta_value like '%supercontig%'");

    $dbh->do("DELETE FROM $target.coord_system ".
	     "WHERE name like 'supercontig'");
    
    $dbh->do("DELETE $target.assembly ".
	     "FROM $target.assembly a, $target.seq_region sr ". 
	     "WHERE sr.coord_system_id = 2 ".
	     "and a.asm_seq_region_id = sr.seq_region_id");

    $dbh->do("DELETE FROM $target.seq_region ".
	     "WHERE coord_system_id = 2");
}

1;
