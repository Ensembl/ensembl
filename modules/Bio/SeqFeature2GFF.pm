
#
# BioPerl module for Bio::SeqFeature2GFF
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::SeqFeature2GFF - Provides mappings from Bioperl's SeqFeature objects to Tim/Richard's GFF objects

=head1 SYNOPSIS

  @genefeatures = Bio::SeqFeature2GFF->convert(@seqfeatures);

=head1 DESCRIPTION

Factory method which translates SeqFeature objects across to 
GFF::GeneFeature objects.

This module will become less necessary (if at all) as the
SeqFeature and GFF stuff become aligned.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::SeqFeature2GFF;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;



=head2 convert

 Title   : convert
 Usage   : @gf = Bio::SeqFeature2GFF->convert(@sf)
 Function: Converts Bio::SeqFeatureI type objects to GFF::GeneFeature
           objects
 Example :
 Returns : 
 Args    :


=cut

sub convert{
   my ($class,@seqfeatures) = @_;

   my @ret;
   foreach my $sf ( @seqfeatures ) {
       my $gf = new GFF::GeneFeature;

       $gf->feature( $sf->primary_tag() );
       $gf->source ( $sf->source_tag() );

       $gf->start($sf->start());
       $gf->end($sf->end());
       if( $sf->strand == 1 ) {
	   $gf->strand('+');
       } elsif ( $gf->strand == -1 ) {
	   $gf->strand('-');
       } else {
	   $gf->strand('.');
       }

       if( $sf->can('score') ) {
	   $gf->score($sf->score());
       } else {
	   $gf->score('.');
       }

       if( $sf->can('frame') ) {
	   $gf->frame($sf->frame());
       } else {
	   $gf->frame('.');
       }

       push(@ret,$gf);
   }
   return @ret;
}

1;
