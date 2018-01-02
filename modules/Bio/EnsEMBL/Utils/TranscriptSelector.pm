=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::TranscriptSelector - Finds canonical transcripts

=head1 SYNOPSIS

  my $selector = Bio::EnsEMBL::Utils::TranscriptSelector->new($ccds_dba);
  my $canonical_transcript = $selector->select_canonical_transcript_for_Gene($gene);

=head1 DESCRIPTION

  The decision process for choosing a canonical transcript of a given Gene is
  an involved process. This package converts transcript attributes into
  numeric values, sorts the values and returns the favourite transcript.
  
  The canonical order of precedence is as follows:
    longest translation of transcript present in CCDS that is reference sequence
    longest translation of protein-coding transcript
    longest translation of transcript marked nonsense-mediated-decay
    longest translation of any other transcript (premature stop codon translations have an effective length of 0)
    longest non-coding transcript
    first stable ID in alphabetical order
    
  The last condition is to give consistent behaviour when everything is else is equal.
  It selects the "older" stable ID, thus preventing new IDs supplanting old ones that
  remain correct.
=cut

package Bio::EnsEMBL::Utils::TranscriptSelector;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;

=head2 new

  Arg [1]    : Optional - CCDS database adaptor - needed for species with CCDS only
  Arg [2]    : Optional - Boolean verbose flag. Turn on to fill your logs
  Description: Constructor
  Returntype : Bio::EnsEMBL::Utils::TranscriptSelector

=cut

sub new {
    my $class = shift;
    my $self = {
                'ccds_dba' => shift,
                'verbose' => shift,
    };
    bless $self, $class;
    if (not defined ($self->{'ccds_dba'}) ) { warning ("Running without CCDS DB");}
    return $self;
}

=head2 select_canonical_transcript_for_Gene

 Arg 1      : Bio::EnsEMBL::Gene 
 Example    : $canonical_transcript = $selector->select_canonical_transcript_for_Gene
 Description: Sorts the Transcripts of this Gene into order of suitability,
              and returns the favourite Transcript.
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : 

=cut

sub select_canonical_transcript_for_Gene {
    my $self = shift;
    my $gene = shift;

    if (!$gene) {throw ('Cannot select canonical transcript without a gene.');}

    my $transcript_array = $gene->get_all_Transcripts;
    my @transcripts;
    if ($transcript_array) {
        @transcripts = @$transcript_array;
    } else {
        warning('No transcripts attached to gene '.$gene->stable_id);
        return;
    }
    my @encoded; # array of encoded transcripts
    
    foreach my $transcript (@transcripts) {
        my $encoded_transcript = $self->encode_transcript($transcript); 
        push(@encoded, $encoded_transcript );
        if ($self->{'verbose'}) { 
            printf "%s encoded to: [%s,%s,%s,%s,%s,%s,%s]\n",$transcript->stable_id, @$encoded_transcript;
        }
    }

    my $sorted_ids = $self->sort_into_canonical_order(\@encoded);
    if ($self->{'verbose'}) { 
        print "Sorted order: \n";
        foreach my $dbID (@$sorted_ids) {
            print $dbID."\n";
        }
    }
    my $canonical_dbID = $sorted_ids->[0];
    
    foreach my $transcript (@transcripts) {
        if ($transcript->dbID == $canonical_dbID) {
            if ($self->{'verbose'}) {print 'Chosen transcript: '.$transcript->stable_id."\n";}
            return $transcript;
        }
    }
    throw ('Run out of transcripts without finding selected canonical dbID.')
}


# Constants for doing numerical sorting of transcripts based upon which ones we prefer. Lowest is best. 
my %source_priority = ('ccds' => 1,
                       'merged' => 2,
                       'other' => 3);
my %biotype_priority = ('protein_coding' => 1,
                        'nonsense_mediated_decay' => 2,
                        'non_stop_decay' => 2,
                        'polymorphic_pseudogene' => 2,
                        'other' => 3,
);

=head2 encode_transcript

 Arg 1      : Transcript
 Description: Converts a transcript into a list of encoded values for sorting
              Priorities are defined immediately above
              Unimportant biotypes and sources are classed as 'other'
 Returntype : Listref of encoded attributes
=cut

sub encode_transcript {
    my $self = shift;
    my $transcript = shift;

    my $type;
    if ( $self->{'ccds_dba'} && $transcript->slice->is_reference() 
            && $self->check_Ens_trans_against_CCDS($transcript) ) {
        $type = 'ccds';
    } elsif ($transcript->analysis->logic_name eq 'ensembl_havana_transcript') {
        $type = 'merged';
    } else {
        $type = 'other';
    }
        
    my $translation = $transcript->translate;
    my $translation_length = 0;
    if ($translation) { 
        $translation_length = $translation->length;
        # Translations containing premature stops are undesirable.
        if ($translation->seq =~ /\*/) {$translation_length = 0;}
    }
    # Zero-length/non-existent translations are ok. We sort by coding and non-coding first
    my $translates = 0;
    if ($translation_length != 0) {$translates = 1;}
    
    my $transcript_length = $transcript->length;
    
    my $biotype;
    if (   $transcript->biotype() ne 'protein_coding' 
        && $transcript->biotype() ne 'nonsense_mediated_decay'
        && $transcript->biotype() ne 'non_stop_decay'
        && $transcript->biotype() ne 'polymorphic_pseudogene' ) {
        $biotype = 'other';
    } else { $biotype = $transcript->biotype(); }

    return [$transcript->dbID,
            $translates,
            $source_priority{$type}, 
            $biotype_priority{$biotype}, 
            $translation_length, 
            $transcript_length,
            $transcript->stable_id];
}


=head2 sort_into_canonical_order

  Arg 1      : 2D array reference of numerically encoded values
                        0              1        2         3             4                  5               6
              ( [transcript dbID, translates, source , biotype, translation length, transcript length, stable ID],
                ...
              ) 
  Description: see Schwartzian transform for method in the following madness:
               sort the 6-column array by the last 5 columns, then map the first elements
               into a list of dbIDs, now in canonical order. 
  Returntype : Listref of ensembl dbIDs
  Caller     : select_canonical_transcript_for_Gene
=cut
 
sub sort_into_canonical_order {
    my $self = shift;
    my $encoded_list_ref = shift;
    
    my @sorted_ids = map { $_->[0] }    
        sort {
            # [0] contains ID
          $b->[1] <=> $a->[1] ||    # translates
          $a->[2] <=> $b->[2] ||    # source
          $a->[3] <=> $b->[3] ||    # biotype
          $b->[4] <=> $a->[4] ||    # translation length (largest is best, $a and $b reversed)
          $b->[5] <=> $a->[5] ||    # transcript length               "
          $a->[6] cmp $b->[6]       # stable id. All other things being equal, 'lowest' transcript ID wins
        } @{$encoded_list_ref};
    
    return \@sorted_ids;
}

=head2 check_Ens_trans_against_CCDS

  Arg 1      : Transcript 
  Description: Attempts to find a matching transcript in CCDS by comparing Exon
               composition. Returns true if one is found, or silently ends. 
  Returntype : Boolean
  Caller     : encode_transcript
=cut

sub check_Ens_trans_against_CCDS {
    my ( $self ,$transcript ) = @_;
   
    my @translateable_exons = @{ $transcript->get_all_translateable_Exons };
    
    my $seq_region_name = $transcript->slice->seq_region_name;
    my $seq_region_start = $transcript->seq_region_start;
    my $seq_region_end = $transcript->seq_region_end;
    my $ccds_dba = $self->{'ccds_dba'};
    
    my $ext_slice = $ccds_dba->get_SliceAdaptor->fetch_by_region(
                                       'toplevel',
                                       $seq_region_name,
                                       $seq_region_start,
                                       $seq_region_end );

    EXT_GENE: foreach my $ext_gene ( @{ $ext_slice->get_all_Genes } ) {
        EXT_TRANS: foreach my $ext_trans ( @{ $ext_gene->get_all_Transcripts } ) {
            my @ext_exons = @{ $ext_trans->get_all_Exons };
    
            if ( scalar(@translateable_exons) == scalar(@ext_exons) ) {
                for ( my $i = 0 ; $i < scalar(@translateable_exons) ; $i++ ) {
                    if ( $translateable_exons[$i]->coding_region_start($transcript) 
                        != $ext_exons[$i]->seq_region_start
                        || $translateable_exons[$i]->strand 
                        != $ext_exons[$i]->strand
                        || $translateable_exons[$i]->coding_region_end($transcript) 
                        != $ext_exons[$i]->seq_region_end 
                    ) {
                        next EXT_TRANS;
                    }
                }
                print "Ensembl transcript "
                    . $transcript->display_id
                    . " found match "
                    . $ext_gene->display_id
                    . " in CCDS DB.\n" if $self->{'verbose'};
                if ($ext_gene->stable_id !~ /^CCDS/) {
                    throw (sprintf("Database does not appear to contain CCDS IDs. Possible configuration problem with CCDS source. Found ID %s", $ext_gene->stable_id()));
                }
                return 1;
            }
        }    # end foreach EXT_TRANS
    }    # end foreach EXT_GENE
} ## end sub check_Ens_trans_against_CCDS

1;
