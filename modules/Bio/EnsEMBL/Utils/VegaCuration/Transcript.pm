package Bio::EnsEMBL::Utils::VegaCuration::Transcript;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::VegaCuration::Gene;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Utils::VegaCuration::Gene);


=head2 find_non_overlaps

   Args       : arrayref of B::E::Transcripts
   Example    : find_non_overlaps($all_transcripts)
   Description: identifies any gaps between transcripts
   Returntype : array refs of stable IDs
   Exceptions : none

=cut

sub find_non_overlaps {
	my $self = shift;
	my ($all_transcripts) = @_;
	my $non_overlaps = [];
	foreach my $transcript1 (@{$all_transcripts}) {
#		warn $transcript1->stable_id;
		foreach my $transcript2 (@{$all_transcripts}) {
#			warn "\t".$transcript2->stable_id;
			if ($transcript1->end < $transcript2->start) {
				push @{$non_overlaps}, $transcript1->stable_id;
				push @{$non_overlaps}, $transcript2->stable_id;
			}
		}
	}
	return $non_overlaps;
}

=head2 update names

   Arg[1]     : B::E::Gene
   Arg[2]     : FH
   Arg[3]     : counter 1
   Arg[4]     : counter 2
   Example    : $support->update_names($gene,$fh,\$c1,\$c2)
   Description: - patches transcripts with identical names according to presence
                of CDS and length. Only does so if there are no fragmented remarks
                - adds remark attribute to gene
                - writes IDs of previously seen genes to file
   Returntype : true | false, counter1, counter2

=cut

sub update_names {
	my $self = shift;
	my ($gene,$k_flist_fh,$c1,$c2) = @_;
	my $aa  = $gene->adaptor->db->get_AttributeAdaptor;
	my $dbh = $gene->adaptor->db->dbc->db_handle;

	#get list of IDs that have previously been sent to annotators
	my $seen_genes = $self->get_havana_comments;

	my $gsi    = $gene->stable_id;
	my $gid    = $gene->dbID;
	my $g_name;
	eval {
		$g_name = $gene->display_xref->display_id;
	};	
	if ($@) {
		$g_name = $gene->get_all_Attributes('name')->[0]->value;
	}
	my $gene_remark = 'This locus has been annotated as fragmented because either there is not enough evidence covering the whole locus to identify the exact exon structure of the transcript, or because the transcript spans a gap in  the assembly';
	my $attrib = [
		Bio::EnsEMBL::Attribute->new(
			-CODE => 'remark',
			-NAME => 'Remark',
			-DESCRIPTION => 'Annotation remark',
			-VALUE => $gene_remark,
		) ];
	#get existing gene and transcript remarks
	my %remarks;
	foreach my $type ('remark','hidden_remark') {
		$remarks{$type}->{'gene'} = [ map {$_->value} @{$gene->get_all_Attributes($type)} ];
		foreach my $trans (@{$gene->get_all_Transcripts()}) {
			my $tsi = $trans->stable_id;
			push @{$remarks{$type}->{'transcripts'}}, map {$_->value} @{$trans->get_all_Attributes('remark')};
		}
	}

	#see if any of the remarks identify this gene as being known by Havana as being fragmented
	if ( (grep {$_ eq 'fragmented_locus'} @{$remarks{'hidden_remark'}->{'gene'}})
			 || (grep {$_ =~ /fragmen/} @{$remarks{'remark'}->{'transcripts'}})
				 || (grep {$_ =~ /fragmen/} @{$remarks{'hidden_remark'}->{'transcripts'}})
			 ) {
		if (grep { $_ eq $gene_remark} @{$remarks{'remark'}->{'gene'}}) {
			$self->log("Fragmented loci annotation remark for gene $gid already exists\n");
		}
		#add gene_attrib 
		else {
			if (! $self->param('dry_run') ) {
				$aa->store_on_Gene($gid,$attrib);
			}			
			$self->log("Added correctly formatted fragmented loci annotation remark for gene $gsi\n");
		}
		return (0,$c1,$c2);
	}
	#log if it's been reported before since the gene should should have a remark.
	elsif ($seen_genes->{$gsi} eq 'fragmented') {
		$self->log_warning("PREVIOUS: Added correctly formatted fragmented loci annotation remark for gene $gsi (has previously been OKeyed by Havana as being fragmented but has no Annotation remark, please add one!)\n");
		print $k_flist_fh "$gsi\n";
		#add gene_attrib anyway.
		if (! $self->param('dry_run') ) {
			$aa->store_on_Gene($gid,$attrib);
		}
		return (0,$c1,$c2);
	}

	#patch transcript names if no remarks found
	else {
		$c1++;
		my @trans = $gene->get_all_Transcripts();

		#separate coding and non_coding transcripts
		my $coding_trans = [];
		my $noncoding_trans = [];
		foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
			if ($trans->translate) {
				push @$coding_trans, $trans;
			}
			else {
				push @$noncoding_trans, $trans;
			}
		}

		#sort transcripts coding > non-coding, then on length
		my $c = 0;
		$self->log("\nPatching names according to CDS and length:\n",1);
		foreach my $array_ref ($coding_trans,$noncoding_trans) {
			foreach my $trans ( sort { $b->length <=> $a->length } @$array_ref ) {
				my $tsi = $trans->stable_id;
				my $t_name;
				eval {
					$t_name = $trans->display_xref->display_id;
				};	
				if ($@) {
					$t_name = $trans->get_all_Attributes('name')->[0]->value;
				}
				$c++;
				my $ext = sprintf("%03d", $c);
				my $new_name = $g_name.'-'.$ext;
				$self->log(sprintf("%-20s%-3s%-20s", "$t_name ", "-->", "$new_name")."\n",1);
				if (! $self->param('dry_run')) {

					# update transcript display xref
					$dbh->do(qq(UPDATE  xref x, external_db edb
                                SET     x.display_label  = "$new_name"
                                WHERE   x.external_db_id = edb.external_db_id
                                AND     x.dbprimary_acc  = "$tsi"
                                AND     edb.db_name      = "Vega_transcript"));
				}
				$c2++;
			}
		}
	}
	return (1,$c1,$c2);
}

#check for duplicated loutre names and distinguish between overlapping and non overlapping
#transcripts

sub check_names_and_overlap {
	my $self = shift;
	my ($transcript_info,$gene,$n_flist_fh) = @_;
	my $ta  = $gene->adaptor->db->get_TranscriptAdaptor;
	my $gsi = $gene->stable_id;
	my $g_name = $gene->get_all_Attributes('name')->[0]->value;
#	warn Dumper($transcript_info);# if ($gsi eq 'OTTHUMG00000019604');
	foreach my $set (values %{$transcript_info} ) {
		next if (scalar @{$set} == 1);
		my $transcripts = [];
		my $all_t_names;

		#check for identical names
		my $duplicate = 0;
#		if ($gsi eq 'OTTHUMG00000008082') {
#			warn Dumper($set);
#		}
		foreach my $id1 (@{$set}) {
			my ($name1,$tsi1) = split /\|/, $id1;
#			if ($gsi eq 'OTTHUMG00000008082') {
#				my $t = $ta->fetch_by_stable_id($tsi1);
#				warn "$tsi1--",ref($t);
#			}
			push @{$transcripts} , $ta->fetch_by_stable_id($tsi1);
			$all_t_names .= "$tsi1 ";
			foreach my $id2 (@{$set}) {
				my ($name2,$tsi2) = split /\|/, $id2;
				next if ($tsi1 eq $tsi2);
				$duplicate = 1 if ( $name1 eq $name2);
			}
		}
		if ($duplicate) {				
			$self->log_warning("IDENTICAL: Gene $gsi ($g_name) has transcripts with identical loutre names, please fix\n");
		}

		my $non_overlaps;
		eval {
			$non_overlaps = $self->find_non_overlaps($transcripts);
		};
		if ($@) {
			$self->log_warning("Problem looking for overlapping transcripts for gene $gsi (is_current = 0 ?). Skipping this bit\n");
		}
		elsif (@{$non_overlaps}) {
			my $tsi_string = join ' ', @{$non_overlaps};
	
			#log gsi (to be sent to Havana) if the transcripts don't overlap				
			$self->log_warning("NEW: Non-overlapping: $gsi ($g_name) has non-overlapping transcripts ($tsi_string) with duplicated names, and it has no \'Annotation_remark- fragmented_loci\' on the gene or \'\%fragmen\%\' remark on any transcripts. Neither has it been OKeyed by Havana before. Transcript names are being patched but this needs checking by Havana.\n");
			print $n_flist_fh "$gsi\n";
		}
		elsif ($self->param('verbose')) {
			$self->log_warning("NEW: Overlapping: $gsi ($g_name) has overlapping transcripts ($all_t_names) with duplicated names and it has no \'Annotation_remark- fragmented_loci\' on the gene or \'\%fragmen\%\' remark on any transcripts. Neither has it been OKeyed by Havana before. Transcript names are being patched but this could be checked by Havana if they were feeling keen.\n");
			print $n_flist_fh "$gsi\n";
		}
	}
}		

sub get_havana_comments {
	my $seen_genes;
	while (<DATA>) {
		next if /^\s+$/ or /#+/;
		my ($obj,$comment) = split /=/;
		$obj =~ s/^\s+|\s+$//g;
		$comment =~ s/^\s+|\s+$//g;
		$seen_genes->{$obj} = $comment;
	}
	return $seen_genes;
}



#details of genes with duplicated transcript names that have already been reported to Havana
#identified as either fragmented or as being OK to patch
__DATA__

OTTMUSG00000005478 = fragmented
OTTMUSG00000001936 = fragmented
OTTMUSG00000017081 = fragmented
OTTMUSG00000011441 = fragmented
OTTMUSG00000013335 = fragmented
OTTMUSG00000011654 = fragmented
OTTMUSG00000001835 = fragmented
OTTMUSG00000012302 =
OTTMUSG00000013368 =
OTTMUSG00000015766 =
OTTMUSG00000016025 =
OTTMUSG00000001066 =
OTTMUSG00000016331 =
OTTMUSG00000006935 =
OTTMUSG00000007263 =
OTTMUSG00000000304 =
OTTMUSG00000009150 =
OTTMUSG00000008023 =
OTTMUSG00000017077 =
OTTMUSG00000003440 =
OTTMUSG00000016310 =
OTTMUSG00000026199 =
OTTMUSG00000028423 =
OTTMUSG00000007427 =
