
#
# BioPerl module for Protein_FeaturePair
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Protein_FeaturePair.pm - DESCRIPTION of Object

=head1 SYNOPSIS

my $feature = new Bio::EnsEMBL::Protein_FeaturePair(-feature1 => $feat1,
						    -feature2 => $feat2,);

=head1 DESCRIPTION

This object inherits from Bio::EnsEMBL::FeaturePair. This extension has been implemented to work with the Protein object. Each Protein Feature should be stored in a Protein_FeaturePair object.

=head1 CONTACT

mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Protein_FeaturePair;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::FeaturePair;

@ISA = qw(Bio::EnsEMBL::FeaturePair);



=head2 to_FTHelper

 Title   : to_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub to_FTHelper{
   my ($self) = @_;

   # Make new FTHelper, and fill in the key
   my $fth = Bio::SeqIO::FTHelper->new;

#The key should not be always domain, its currently true because we store in Protein Features only Interpro hits but we should get the key information from the analysis table...but these is no column where this key could be stored...
   $fth->key('Domain');
   
   # Add location line
   my $g_start = $self->start;
   my $g_end   = $self->end;
   my $loc = "$g_start..$g_end";
   if ($self->strand == -1) {
        $loc = "complement($loc)";
    }
   $fth->loc($loc);
   
   # Add note describing similarity
   my $type    = $self->hseqname;
   my $r_start = $self->hstart;
   my $r_end   = $self->hend;
   $fth->add_field('note', "$type: matches $r_start to $r_end");
   $fth->add_field('note', "score=".$self->score);
   
   my $subject = $self->hseqname;

   $fth->add_field('description', $subject);

   return $fth;
}


