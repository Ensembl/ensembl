
#
# BioPerl module for ProteinFeature
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

ProteinFeature.pm - DESCRIPTION of Object

=head1 SYNOPSIS

my $feature = new Bio::EnsEMBL::ProteinFeature(-feature1 => $feat1,
					       -feature2 => $feat2,);

=head1 DESCRIPTION

This object inherits from Bio::EnsEMBL::FeaturePair. This extension has been implemented to work with the Protein object. Each Protein Feature should be stored in a Protein_FeaturePair object.

This object was formerly named Protein_FeaturePair.

=head1 CONTACT

mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::ProteinFeature;
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
  
#This information (description of the protein feature, eg: Domain, exon, ...) is stored in gff_feature. Obviously shouldn't be but here waiting for a new schema.
 
   my $desc = $self->analysis->gff_feature;
  
   
   $fth->key($desc);
   
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

=head2 hdbname

 Title   : hdbname
 Usage   : $obj->hdbname($newval)
 Function: 
 Returns : value of hdbname
 Args    : newvalue (optional)


=cut

sub hdbname{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'hdbname'} = $value;
    }
    return $obj->{'hdbname'};

}

=head2 idesc

 Title   : idesc
 Usage   : $obj->idesc($newval)
 Function: 
 Returns : value of idesc
 Args    : newvalue (optional)


=cut

sub idesc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'idesc'} = $value;
    }
    return $obj->{'idesc'};

}

=head2 interpro_ac

 Title   : interpro_ac
 Usage   : $obj->interpro_ac($newval)
 Function: 
 Returns : value of interpro_ac
 Args    : newvalue (optional)


=cut

sub interpro_ac{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'interpro_ac'} = $value;
    }
    return $obj->{'interpro_ac'};

}


=head2 intron_length

 Title   : intron_length
 Usage   : my $length = $intronfeature->intron_length
 Function: Return the length of an intron by calling its starting point and end point in global coordinates
 Example :
 Returns : Length of a given intron
 Args    : Nothing


=cut

sub intron_length{
   my ($self) = @_;
   my $start = $self->feature2->start;
   my $end = $self->feature2->end;
   my $length = $end - $start;
   return $length;
}

=head2 intron_position

 Title   : intron_position
 Usage   :
 Function: Return the position of the intron on the amino acid sequence
 Example :
 Returns : 
 Args    :


=cut

sub intron_position{
   my ($self,@args) = @_;
   my $pos = $self->feature1->start;
   
   return $pos;

}

1;
