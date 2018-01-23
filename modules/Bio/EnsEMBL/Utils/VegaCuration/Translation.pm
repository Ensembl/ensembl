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
=head1 SYNOPSIS
=head1 DESCRIPTION
=head1 METHODS
=cut

package Bio::EnsEMBL::Utils::VegaCuration::Translation;
use strict;
use warnings;
use vars qw(@ISA);
use Data::Dumper;
use Bio::EnsEMBL::Utils::VegaCuration::Transcript;
@ISA = qw(Bio::EnsEMBL::Utils::VegaCuration::Transcript);

=head2 check_CDS_start_end_remarks

   Args       : B::E::Transcript
   Example    : my $results = $support->check_CDS_end_remarks($transcript)
   Description: identifies incorrect 'CDS end...' transcript remarks in a
                otter-derived Vega database
   Returntype : hashref

=cut

sub check_CDS_start_end_remarks {
  my $self = shift; 
  my $trans = shift;
  # info for checking
  my @remarks = @{$trans->get_all_Attributes('remark')};
  my $coding_end   = $trans->cdna_coding_end;
  my $coding_start = $trans->cdna_coding_start;
  my $trans_end    = $trans->length;
  my $trans_seq    = $trans->seq->seq;
  my $stop_codon   = substr($trans_seq, $coding_end-3, 3);
  my $start_codon  = substr($trans_seq, $coding_start-1, 3);
  #hashref to return results
  my $results;
  
  #extra CDS end not found remarks
  if (grep {$_->value eq 'CDS end not found'} @remarks) {
    if (   ($coding_end != $trans_end) 
	     && ( grep {$_ eq $stop_codon} qw(TGA TAA TAG) ) ) {
      $results->{'END_EXTRA'} = 1;
    }
  }
  #missing CDS end not found remark
  if ( $coding_end == $trans_end ) {
    if (! grep {$_->value eq 'CDS end not found'} @remarks) {
      if (grep {$_ eq $stop_codon} qw(TGA TAA TAG)) {
	$results->{'END_MISSING_2'} = 1;
      }
      else {
	$results->{'END_MISSING_1'} = $stop_codon;
      }
    }
  }
  #extra CDS start not found remark
  if (grep {$_->value eq 'CDS start not found'} @remarks) {
    if (   ($coding_start != 1) 
	     && ($start_codon eq 'ATG') ) {
      $results->{'START_EXTRA'} = 1;
    }
  }
  #missing CDS start not found remark
  if ( $coding_start == 1) {
    if ( ! grep {$_->value eq 'CDS start not found'} @remarks) {
      if ($start_codon eq 'ATG') {
	$results->{'START_MISSING_2'} = 1;
      } else {
	$results->{'START_MISSING_1'} = $start_codon;
      }
    }
  }
  return $results;
}

=head2 check_CDS_start_end_remarks_loutre

   Args       : B::E::Transcript
   Example    : my $results = $support->check_CDS_end_remarks($transcript)
   Description: identifies incorrect 'CDS end...' transcript attribs in a loutre
                of a loutre-derived Vega database.
   Returntype : hashref

=cut

sub check_CDS_start_end_remarks_loutre {
  my $self = shift;
  my $trans = shift;

  # info for checking
  my @stops = qw(TGA TAA TAG);
  my %attributes;
  foreach my $attribute (@{$trans->get_all_Attributes()}) {
    push @{$attributes{$attribute->code}}, $attribute;
  }
#  warn $trans->stable_id;
#  warn Data::Dumper::Dumper(\%attributes);
  my $coding_end   = $trans->cdna_coding_end;
  my $coding_start = $trans->cdna_coding_start;
  my $trans_end    = $trans->length;
  my $trans_seq    = $trans->seq->seq;
  my $stop_codon_offset = 3 + $trans->translation->end_Exon->end_phase;
  my $initial_exon_phase = $trans->translation->start_Exon->phase;
  my $stop_codon   = substr($trans_seq, $coding_end-3, 3);
  my $start_codon  = substr($trans_seq, $coding_start-1, 3);

  my $start_codon_incorrect = 1;
  if ($start_codon eq 'ATG' ) {
    $start_codon_incorrect = 0;
  }
  elsif ($start_codon eq 'CTG') {
    foreach my $attrib (@{$attributes{'remark'}}) {
      if ($attrib->value =~ /non[- ]ATG start/) {
	$start_codon_incorrect = 0;
      }
    }
  }

#  warn "$start_codon -- $initial_exon_phase -- $coding_start -- $start_codon_incorrect";

  #hashref to return results
  my $results;

  #extra CDS end not found remarks
  if ($attributes{'cds_end_NF'}) {
    if ( ($attributes{'cds_end_NF'}->[0]->value == 1)
	   && ($coding_end != $trans_end) 
	   && ( grep {$_ eq $stop_codon} @stops) ) {
#      warn $trans->stable_id.": $coding_end--$trans_end--$stop_codon";
#      warn $trans->translation->end_Exon->end_phase;
      $results->{'END_EXTRA'} = $stop_codon;
    }
  }
  #missing CDS end not found remark
  if ( $coding_end == $trans_end ) {
    if ($attributes{'cds_end_NF'}) {
      if ($attributes{'cds_end_NF'}->[0]->value == 0 ) {
	if (! grep {$_ eq $stop_codon} @stops) {
#	  warn $trans->translation->end_Exon->end_phase;
	  $results->{'END_MISSING'}{'WRONG'} = $stop_codon;
	}
      }
    }
    elsif (! grep {$_ eq $stop_codon} @stops) {
      $results->{'END_MISSING'}{'ABSENT'} = $stop_codon;
    }
  }
  #extra CDS start not found remark 
  if ( $attributes{'cds_start_NF'}) {
    if (   ($attributes{'cds_start_NF'}->[0]->value == 1 )
        && (! $start_codon_incorrect)) {
      unless ( ($coding_start == 1) && ( $initial_exon_phase > 0)) {
	$results->{'START_EXTRA'} = $start_codon;
      }
    }
  }
  #missing CDS start not found remark
  if ( $coding_start == 1) {
    if ( $attributes{'cds_start_NF'} ) {
      if ( $attributes{'cds_start_NF'}->[0]->value == 0 ) {
	if ($start_codon_incorrect) {
	  $results->{'START_MISSING'}{'ABSENT'} = $start_codon;
	}
	elsif ($initial_exon_phase > 0) {
	  $results->{'START_MISSING'}{'INITIAL_PHASE'} = $initial_exon_phase;
	}
      }
    }
    elsif ($start_codon ne 'ATG') {
      $results->{'START_MISSING'}{'ABSENT'} = $start_codon;
    }

  }
  return $results;
}

=head2 get_havana_seleno_comments

   Args       : none
   Example    : my $results = $support->get_havana_seleno_comments
   Description: parses the HEREDOC containing Havana comments in this module
   Returntype : hashref

=cut

sub get_havana_seleno_comments {
  my $seen_translations;
  while (<DATA>) {
    next if /^\s+$/ or /#+/;
    my ($obj,$comment) = split /=/;
    $obj =~ s/^\s+|\s+$//g;
    $comment =~ s/^\s+|\s+$//g;
    # We add the origin as now "seen" can come from a number of places, and have
    # a number of consequences in different cases, not just discounted Secs from this method. -- ds23
    $seen_translations->{$obj} = [ $comment,"notsec-havana" ];
  }
  return $seen_translations;
}

sub check_for_stops {
  my $support = shift;
  my ($gene,$seen_transcripts,$log_object) = @_;
  my $transcripts;
  my $has_log_object=defined($log_object);
  if($has_log_object){
    my @help = $log_object->species_params->get_trans($gene->stable_id);
    $transcripts=\@help;
  }else{
    $log_object=$support;
    $transcripts=$gene->get_all_Transcripts;
  }

  my $gname = $gene->get_all_Attributes('name')->[0]->value;
  my $gsi = $gene->stable_id;
  my $scodon = 'TGA';
  my $mod_date = $support->date_format( $gene->modified_date,'%d/%m/%y' );
  my $hidden_remak_ttributes;
 TRANS:
  foreach my $trans (@$transcripts) {
    my $tsi = $trans->stable_id;
    my $tID = $trans->dbID;
    my $tname = $trans->get_all_Attributes('name')->[0]->value;
    if($has_log_object){
      $hidden_remak_ttributes=$log_object->species_params->get_attributes->{$tsi}->{'hidden_remark'};
    }else{
      $hidden_remak_ttributes=$trans->get_all_Attributes('hidden_remark');
    }
    foreach my $rem (@$hidden_remak_ttributes) {
      if ($rem->value =~ /not_for_Vega/) {
        #$support->log_verbose("Skipping transcript $tname ($tsi) since 'not_for_Vega'\n",1);
        $log_object->_save_log('log_verbose', '', $gsi, '', $tsi, '', "Skipping transcript $tname ($tsi) since 'not_for_Vega'");
        next TRANS;
      }
    }

    # go no further if there is a ribosomal framshift attribute
    foreach my $attrib (@{$trans->get_all_Attributes('_rib_frameshift')}) {
      if ($attrib->value) {
        $log_object->_save_log('log', '', $gsi, '', $tsi, '', "Skipping $tsi ($tname) since it has a ribosomal frameshift attribute");
        next TRANS;
      }
    }
    foreach my $attrib (@{$trans->get_all_Attributes('hidden_remark')}) {
      if ($attrib->value eq 'ribosomal_frameshift') {
        $log_object->_save_log('log', '', $gsi, '', $tsi, '', "Skipping $tsi ($tname) since it has a ribosomal frameshift hidden_remark");
        next TRANS;
      }
    }

    #$support->log_verbose("Studying transcript $tsi ($tname, $tID)\n",1);
    $log_object->_save_log('log_verbose', '', $gsi, '', $tsi, '', "Studying transcript $tsi ($tname, $tID)");
    my $peptide;
		
    # go no further if the transcript doesn't translate or if there are no stops
    next TRANS unless ($peptide = $trans->translate);
    my $pseq = $peptide->seq;
    my $orig_seq = $pseq;
    # (translate method trims stops from sequence end)

    next TRANS unless ($pseq =~ /\*/);
    next TRANS if ($trans->biotype eq 'polymorphic_pseudogene' or
                   $trans->biotype =~ /processed_pseudogene$/);

    # warn sprintf("Stop codon is '%s'\n",substr($trans->translateable_seq,-3));
    #$support->log_verbose("Stops found in $tsi ($tname)\n",1);
    $log_object->_save_log('log_verbose', '', $gsi, '', $tsi, '', "Stops found in $tsi ($tname)");

    # find out where and how many stops there are
    my @found_stops;
    my $mrna   = $trans->translateable_seq;
    my $offset = 0;
    my $tstop;
    while ($pseq =~ /^([^\*]*)\*(.*)/) {
      my $pseq1_f = $1;
      $pseq = $2;
      my $seq_flag = 0;
      $offset += length($pseq1_f) * 3;
      my $stop = substr($mrna, $offset, 3);
      my $aaoffset = int($offset / 3)+1;
      push(@found_stops, [ $stop, $aaoffset ]);
      $tstop .= "$aaoffset ";
      $offset += 3;
    }
		
    # are all stops TGA...?
    my $num_stops = scalar(@found_stops);
    my $num_tga = 0;
    my $positions;
    foreach my $stop (@found_stops) {
      $positions .= $stop->[0]."(".$stop->[1].") ";
      if ($stop->[0] eq $scodon) {
	$num_tga++;
      }
    }
    my $source = $gene->source;
    #...no - an internal stop codon error in the database...
    if ($num_tga < $num_stops) {
      if ($source eq 'havana') {
        #$support->log_warning("INTERNAL STOPS HAVANA: Transcript $tsi ($tname) from gene $gname has non \'$scodon\' stop codons [$mod_date]:\nSequence = $orig_seq\nStops at $positions)\n\n");
        $log_object->_save_log('log_warning', '', $gsi, 'TRANSCRIPT', $tsi, 'VQCT_internal_stop', "INTERNAL STOPS HAVANA: Transcript $tsi ($tname) from gene $gname has non \'$scodon\' stop codons [$mod_date]: Sequence = $orig_seq Stops at $positions)");
      }
      else {
        #$support->log_warning("INTERNAL STOPS EXTERNAL: Transcript $tsi ($tname) from gene $gname has non \'$scodon\' stop codons[$mod_date]:\nSequence = $orig_seq\nStops at $positions)\n\n");
        $log_object->_save_log('log_warning', '', $gsi, 'TRANSCRIPT', $tsi, 'VQCT_internal_stop', "INTERNAL STOPS EXTERNAL: Transcript $tsi ($tname) from gene $gname has non \'$scodon\' stop codons[$mod_date]: Sequence = $orig_seq Stops at $positions)");  
      }
    }
    #...yes - check remarks
    else {
      my $flag_remark  = 0; # 1 if word seleno has been used
      my $flag_remark2 = 0; # 1 if existing remark has correct numbering
      my $alabel       = 'Annotation_remark- selenocysteine ';
      my $alabel2      = 'selenocysteine ';
      my $annot_stops;
      my $remarks;
      my $att;
      #get both hidden_remarks and remarks
      foreach my $remark_type ('remark','hidden_remark') {
        if($has_log_object){
          $att=$log_object->species_params->get_attributes->{$trans->stable_id}->{$remark_type};
        }else{
          $att=$trans->get_all_Attributes($remark_type)
        }
        foreach my $attrib ( @$att) {
          push @{$remarks->{$remark_type}}, $attrib->value;
        }
      }
      #parse remarks to check syntax for location of edits
      while (my ($attrib,$remarks)= each %$remarks) {
        foreach my $text (@{$remarks}) {					
          if ( ($attrib eq 'remark') && ($text=~/^$alabel(.*)/) ){
            #$support->log_warning("seleno remark for $tsi stored as Annotation_remark not hidden remark) [$mod_date]\n");
            $log_object->_save_log('log_warning', '', $gsi, '', $tsi, 'VQCT_wrong_selC_coord', "seleno remark for $tsi stored as Annotation_remark not hidden remark) [$mod_date]");        
            $annot_stops=$1;
          }
          elsif ($text =~ /^$alabel2(.*)/) {
            my $maybe = $1;
            if($maybe =~ /^\s*\d+(\s+\d+)*\s*$/) {
              $annot_stops=$maybe;
            } else {
              $log_object->_save_log('log', '', $gene->stable_id, '', $tsi, '', "Maybe annotated stop in incorrect format, maybe just a remark that happens to begin '$alabel2'".
                                                                                " -- might need to investigate: '$alabel2$maybe' [$mod_date]");
            }
          }
        }
      }

      #check the location of the annotated edits matches actual stops in the sequence
      my @annotated_stops;
      if ($annot_stops){
        my $i = 0;
        foreach my $offset (split(/\s+/, $annot_stops)) {
          #OK if it matches a known stop
          if (
            defined($found_stops[$i]) && defined($found_stops[$i]->[1]) && ($found_stops[$i]->[1] == $offset)) {
            push  @annotated_stops, $offset;
          }
          # catch old annotations where number was in DNA not peptide coordinates
          elsif (defined($found_stops[$i]) && defined($found_stops[$i]->[1]) && (($found_stops[$i]->[1] * 3) == $offset)) {
            $log_object->_save_log('log_warning', '', $gene->stable_id, 'DNA', $tsi, 'VQCT_wrong_selC_coord', "DNA: Annotated stop for transcript tsi ($tname) is in DNA not peptide coordinates) [$mod_date]");
          }
          # catch old annotations where number off by one
          elsif (defined($found_stops[$i]) && defined($found_stops[$i]->[1]) && (($found_stops[$i]->[1]) == $offset+1)) {
            $log_object->_save_log('log_warning', '', $gene->stable_id, 'PEPTIDE', $tsi, 'VQCT_wrong_selC_coord', "PEPTIDE: Annotated stop for transcript $tsi ($tname) is out by one) [$mod_date]");
          }
          elsif (defined($offset)  && ($offset=~/^\d+$/)){
            if ($offset == length($orig_seq)+1) {
              if($seen_transcripts->{$tsi} && $seen_transcripts->{$tsi}->[1] eq 'known-tga-stop') {
                $log_object->_save_log('log', '', $gene->stable_id, 'TRANSCRIPT', $tsi, '', "Annotated stop for transcript $tsi ($tname) known to be a stop codon. Ok. [$mod_date]");
              } elsif($seen_transcripts->{$tsi} && $seen_transcripts->{$tsi}->[1] eq 'known-terminal-sec') {
                $log_object->_save_log('log', '', $gene->stable_id, 'TRANSCRIPT', $tsi, '', "Annotated stop for transcript $tsi ($tname) known to be a terminal Sec. Ok. [$mod_date]");
              } else {
                $log_object->_save_log('log_warning', '', $gene->stable_id, 'TRANSCRIPT', $tsi, '', "Annotated stop for transcript $tsi ($tname) \"$offset\" matches actual stop codon yet has no entry in script config to disambiguate it. Please investigate and add appropriate entry to config arrays in add_selcys.pl. [$mod_date]");
              }
            }
            else {
              $log_object->_save_log('log_warning', '', $gene->stable_id, 'TRANSCRIPT', $tsi, 'VQCT_wrong_selC_coord', "Annotated stop for transcript $tsi ($tname) \"$offset\" does not match a TGA codon) [$mod_date]");
              push  @annotated_stops, $offset;
            }
          }						
          $i++;
        }
      }

      #check location of found stops matches annotated ones - any new ones are reported
      foreach my $stop ( @found_stops ) {
        my $pos = $stop->[1];
        my $seq = $stop->[0];
        unless ( grep { $pos == $_} @annotated_stops) {
          if ($seen_transcripts->{$tsi} && $seen_transcripts->{$tsi}->[1] eq 'notsec-havana') {
            #$support->log_verbose("Transcript $tsi ($tname) has potential selenocysteines but has been discounted by annotators:\n\t".$seen_transcripts->{$tsi}.") [$mod_date]\n");
            $log_object->_save_log('log_verbose', '', $gene->stable_id, '', $tsi, 'VQCT_pot_selC', "Transcript $tsi ($tname) has potential selenocysteines but has been discounted by annotators: ".$seen_transcripts->{$tsi}->[0].") [$mod_date]");
          }
          else {
            #$support->log("POTENTIAL SELENO ($seq) in $tsi ($tname, gene $gname) found at $pos [$mod_date]\n");
            $log_object->_save_log('log', '', $gene->stable_id, '', $tsi, 'VQCT_pot_selC', "POTENTIAL SELENO ($seq) in $tsi ($tname, gene $gname) found at $pos [$mod_date]");
          }
        }
      }
    }
  }
}
sub _save_log{
  my $self=shift;
  my $log_type = shift;
  my $chrom_name=shift || '';
  my $gsi=shift || '';
  my $type=shift || '';  
  my $tsi=shift || '';  
  my $tag=shift || '';
  my $txt=shift || '';
  $self->$log_type($txt."\n");
}

#details of annotators comments
__DATA__
OTTHUMT00000144659 = FIXED- changed to transcript
OTTHUMT00000276377 = FIXED- changed to transcript
OTTHUMT00000257741 = FIXED- changed to nmd
OTTHUMT00000155694 = NOT_FIXED- should be nmd but external annotation but cannot be fixed
OTTHUMT00000155695 = NOT_FIXED- should be nmd but external annotation but cannot be fixed
OTTHUMT00000282573 = FIXED- changed to unprocessed pseudogene
OTTHUMT00000285227 = FIXED- changed start site
OTTHUMT00000151008 = FIXED- incorrect trimming of CDS, removed extra stop codon
OTTHUMT00000157999 = FIXED- changed incorrect stop
OTTHUMT00000150523 = FIXED- incorrect trimming of CDS
OTTHUMT00000150525 = FIXED- incorrect trimming of CDS
OTTHUMT00000150522 = FIXED- incorrect trimming of CDS
OTTHUMT00000150521 = FIXED- incorrect trimming of CDS
OTTHUMT00000246819 = FIXED- corrected frame
OTTHUMT00000314078 = FIXED- corrected frame
OTTHUMT00000080133 = FIXED- corrected frame
OTTHUMT00000286423 = FIXED- changed to transcript
OTTMUST00000055509 = FIXED- error
OTTMUST00000038729 = FIXED- corrected frame
OTTMUST00000021760 = FIXED- corrected frame
OTTMUST00000023057 = FIXED- corrected frame
OTTMUST00000015207 = FIXED- corrected frame
OTTMUST00000056646 = FIXED- error
OTTMUST00000059686 = FIXED- corrected frame
OTTMUST00000013426 = FIXED- corrected frame
OTTMUST00000044412 = FIXED- error
