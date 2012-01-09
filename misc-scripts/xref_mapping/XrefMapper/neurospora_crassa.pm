package XrefMapper::neurospora_crassa;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub get_set_lists {

  return [["ExonerateGappedBest5", ["neurospora_crassa","RefSeq_mRNA"]],
	  ["ExonerateGappedBest5", ["neurospora_crassa","RefSeq_mRNA_predicted"]],
	  ["ExonerateGappedBest5", ["neurospora_crassa","RefSeq_ncRNA"]],
	  ["ExonerateGappedBest5", ["neurospora_crassa","RefSeq_ncRNA_predicted"]],
          ["ExonerateGappedBest1", ["neurospora_crassa","*"]]];

}


sub transcript_display_xref_sources {
    my $self     = shift;
    my $fullmode = shift;

    my @list = qw(
                 Uniprot_genename
               );
    
    my %ignore;
    
    return [\@list,\%ignore];
}


1;
