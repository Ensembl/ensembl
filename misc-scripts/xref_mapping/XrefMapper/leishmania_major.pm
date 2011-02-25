package XrefMapper::leishmania_major;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub get_set_lists {

  return [["ExonerateGappedBest1", ["leishmania_major","*"]]];

}


sub transcript_display_xref_sources {
    my $self     = shift;
    my $fullmode = shift;

    my @list = qw(
                 Uniprot_genename
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


1;
