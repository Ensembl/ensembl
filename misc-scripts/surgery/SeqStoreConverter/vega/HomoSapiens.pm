
use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::vega::HomoSapiens;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_attribs {
  my $self = shift;

  #
  # Human clones need their htg phase information copied
  #

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $self->SUPER::create_attribs();
  
  $self->debug("HomoSapiens specific: Creating HTG Phase seq_region attribs");

  $dbh->do
    ("INSERT INTO $target.attrib_type( code, name, description ) " .
     "VALUES ('htg_phase', 'HTG Phase', 'High Throughput Genome Phase')");


  $dbh->do
    ("INSERT INTO $target.seq_region_attrib( seq_region_id, attrib_type_id, " .
                                            "value) " .
     "SELECT tmp_cln.new_id, attrib_type.attrib_type_id, cln.htg_phase " .
     "FROM   $target.tmp_cln_map tmp_cln, $target.attrib_type attrib_type, " .
     "       $source.clone cln " .
     "WHERE  cln.clone_id = tmp_cln.old_id " .
     "AND    attrib_type.code = 'htg_phase'");

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

1;
