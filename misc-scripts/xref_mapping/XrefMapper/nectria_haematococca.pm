=head1 XrefMapper::nectria_haematococca

	Gets created in xref_mapper.pl like this:

		my $mapper = XrefMapper::BasicMapper->process_file($file, !$notverbose, $no_xref);

	the process_file method returns an instance of a species specific xref mapper like this:

		$mapper = "XrefMapper::$module"->new();

=cut
package XrefMapper::nectria_haematococca;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);

=head1 get_set_lists

	Gets called in XrefMapper::SubmitMapper:

	if($self->mapper->can('get_set_lists')){
		@lists =@{$self->mapper->get_set_lists()};  }
	else {
		@lists =@{$self->get_set_lists()};
	}

=cut
sub get_set_lists {

  return [["ExonerateGappedBest1", ["nectria_haematococca","*"]]];

}

=head1 transcript_display_xref_sources

	Gets called in XrefMapper::DisplayXrefs:

		if( $self->mapper->can("transcript_display_xref_sources") ){
			($presedence, $ignore) = @{$self->mapper->transcript_display_xref_sources(1)}; # FULL update mode pass 1
		}
		else{
			($presedence, $ignore) = @{$self->transcript_display_xref_sources(1)}; # FULL update mode pass 1
		}
=cut
sub transcript_display_xref_sources {
    my $self     = shift;
    my $fullmode = shift;

    my @list = qw(
                 Uniprot_genename
                 TRNASCAN_SE
                 RNAMMER
                 RFAM
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
