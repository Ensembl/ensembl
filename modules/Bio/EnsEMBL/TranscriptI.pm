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

A common interface definition to PredictionTranscript and Transcript
to allow these modules to be used interchangeably.

=head1 DESCRIPTION

Creation:

Should not be instantiated. Should only be instantiated through implementing 
subclasses.   

=head1 CONTACT

=head1 AUTHOR - Graham McVicker

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

=head1 APPENDIX

=cut


# Let the code begin...

package Bio::EnsEMBL::TranscriptI;

use strict;

sub get_all_Exons {}

sub coding_start {}

sub coding_end {}

sub start {}

sub end {}

1;
