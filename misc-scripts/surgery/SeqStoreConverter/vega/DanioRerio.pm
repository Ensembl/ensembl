use strict;
use warnings;

use SeqStoreConverter::DanioRerio;

package SeqStoreConverter::vega::DanioRerio;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::DanioRerio);

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
#		     "current_clone_info",
		     "keyword",
		     "job",
		     "job_status",
		     "input_id_analysis");
$self->copy_current_clone_info;
}

sub copy_current_clone_info {
    my $self=shift;
    my $source = $self->source();
    my $target = $self->target();
    my $sth = $self->dbh()->prepare
        ("INSERT INTO $target.current_clone_info(clone_id,clone_info_id) SELECT * FROM $source.current_clone_info");
    $sth->execute();
    $sth->finish();    
}

sub update_clone_info {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  $self->debug("Danio_specific - Transforming clone_id into seq_region_id for clone_info and current_clone_info");

  foreach my $table_name ('clone_info','current_clone_info') {
      my $select_st1 = 
      "SELECT ctg.name, ctg.clone_id " .
      "FROM   $source.contig ctg, $source.$table_name ci " .
      "WHERE  ctg.clone_id = ci.clone_id " .
      "AND ctg.name not like 'ctg%' " .
      "AND ctg.name not like 'NA%'";
      
      my $query_results1 = $dbh->selectall_arrayref($select_st1);
      
      my $i = 0;
      foreach my $contig_name (@$query_results1) {
	  my $embl_acc = $contig_name->[0];
	  my $select_st2 = 
	  "SELECT sr.seq_region_id " .
	  "FROM $target.seq_region sr " . 
	  "WHERE sr.name = '$embl_acc'";
	  my @query_results2 = $dbh->selectrow_array($select_st2);
	  push @{$query_results1->[$i]},@query_results2;
	  $i++;
      }
      
      foreach my $clone (@$query_results1) {
	  my $seq_reg_id = $clone->[2];
	  my $clone_id = $clone->[1];
	  
	  my $update_query = 
	  "UPDATE $target.$table_name " .
	  "SET clone_id = '$seq_reg_id' " .
	  "WHERE clone_id = '$clone_id'";
	  $dbh->do($update_query);
      }
      my $alter_struct_1 = 
      "ALTER table $target.$table_name " .
      "CHANGE clone_id seq_region_id int(10) not null";
      my $alter_struct_2 = 
      "ALTER table $target.$table_name " .
      "add unique index (seq_region_id)";
      $dbh->do($alter_struct_1);
      $dbh->do($alter_struct_2);
  }
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
