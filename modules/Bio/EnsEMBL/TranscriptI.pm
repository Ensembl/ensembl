#
# EnsEMBL module for TranscriptI 
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright EMBL-EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

TranscriptI - Interface definition of EnsEMBL Transcript objects. 

=head1 SYNOPSIS

Used as a common interface definition to PredictionTranscript and Transcript
to allow these modules to be used interchangeably.

=head1 DESCRIPTION

Creation:

Should not be instantiated! Should be instantiated through implementing 
subclasses.   

=head1 CONTACT

=head1 AUTHOR - Graham McVicker

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::TranscriptI;
use strict;

=head2 get_all_Exons

 Title   : get_all_Exons
 Usage   : foreach $exon ( $transI->get_all_Exons)  
 Function: Returns an array of exons in the transcript
           in order, ie the first exon is the 5' most exon
           in the transcript (the one closest to the 5' UTR).
           Currently this is the only defined method in this interface
	   because it is the only one that is needed.  All methods common
	   to PredictionTranscripts and Transcripts should probably be
	   added however.
 Example : my @exons = $tr->get_all_Exons
 Returns : An array of exon objects
 Args    : none

=cut

sub get_all_Exons {
  warn("TranscriptI->get_all_Exons must be overridden by implementing " .
       "subclass\n");

  return ();
}


sub coding_start {
  warn("TranscriptI->coding_start not implemented by subclass\n");
  return undef;
}

sub coding_end {
  warn("TranscriptI->coding_end not implemented by subclass\n");
  return undef;
}

sub start {
  warn("TranscriptI->start not implemented by subclass\n");
  return undef;
}

sub end {
  warn("TranscriptI->end not implemented by subclass\n");
  return undef;
}

1;
