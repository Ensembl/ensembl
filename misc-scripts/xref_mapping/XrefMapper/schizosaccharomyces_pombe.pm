package XrefMapper::schizosaccharomyces_pombe;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);

sub get_set_lists {

  return [["ExonerateGappedBest1", ["schizosaccharomyces_pombe","*"]]];

}



# Pombe is imported from PomBase. The gene and transcript stable IDs
# are the PomBase identifiers. The display_xref_ids for genes and
# transcripts are calculated directly rather than via the more complex
# priority-based method in BasicMapper.pm

sub build_display_xrefs {

  my ($self, $type, $external_db) = @_;

  print "Setting $type display_xrefs from $type stable IDs\n";
  my $dir = $self->core()->dir();

  my $sql = "UPDATE $type t, ${type}_stable_id s, xref x, external_db e SET t.display_xref_id=x.xref_id WHERE t.${type}_id=s.${type}_id AND s.stable_id=x.dbprimary_acc AND e.external_db_id=x.external_db_id AND e.db_name=\'${external_db}\'\n";

  open (SQL, ">$dir/${type}_display_xref.sql");

  print SQL $sql;

  close(SQL);

}


sub transcript_display_xref_sources {
    my $self     = shift;
    my $fullmode = shift;

    my @list = qw(
                PomBase_GENE
                PomBase_TRANSCRIPT
               );
    
    my %ignore;
    
    
    # Both methods
    
    if(!$fullmode){
	$ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
    }
    else{
	$ignore{"EntrezGene"} = 'select ox.object_xref_id from object_xref ox, dependent_xref dx, source s1, xref x1, source s2, xref x2 where ox.object_xref_id = dx.object_xref_id and dx.dependent_xref_id = x1.xref_id and x1.source_id = s1.source_id and s1.name = "EntrezGene" and x2.xref_id = dx.master_xref_id and x2.source_id = s2.source_id and (s2.name like "Refseq_dna_predicted" or s2.name like "RefSeq_peptide_predicted") and ox.ox_status = "DUMP_OUT"';
	
    }
    
    return [\@list,\%ignore];
}

sub gene_description_sources {
  return (
          "PomBase_GENE"
         );
}


sub gene_description_filter_regexps {

  return ();

}



1;
