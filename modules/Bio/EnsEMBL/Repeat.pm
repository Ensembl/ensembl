#
# BioPerl module for Bio::SeqFeature::Homol
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Repeat

=head1 SYNOPSIS

=head1 DESCRIPTION

Extends Bio::EnsEMBL::FeaturePair to store
genomic repeat features.  

Creation:

    my $rep = new Bio::EnsEMBL::Repeat(-start => $start,
				       -end    => $end,
				       -strand => $strand,
				       -source => $source,
				       -primary=> $primary,
				       );

Manipulation:

    my $start = $rep->start;
    my $end   = $rep->end;

    my $repname  = $rep->hseqname;
    my $repstart = $rep->hstart;
    my $repend   = $rep->hend;


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Repeat;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Homol

use Bio::AnnSeqIO::FTHelper;
use Bio::EnsEMBL::FeaturePair;


@ISA = qw(Bio::EnsEMBL::FeaturePair);

# new() is inherited from Bio::SeqFeature::Generic

sub _initialize {
  my($self,@args) = @_;
  my $make = $self->SUPER::_initialize(@args);

  return $make;
}

=head1 Methods unique to Repeat

=head2 to_FTHelper

Called by
C<Bio::AnnSeqIO::FTHelper::from_SeqFeature>,
which is called by
C<Bio::AnnSeqIO::EMBL::write_annseq>.

Returns a C<Bio::AnnSeqIO::FTHelper> object which
will print the repeat in the Sanger Centre
format.

=cut

sub to_FTHelper {
    my( $rep ) = @_;
    
    # Make new FTHelper, and fill in the key
    my $fth = Bio::AnnSeqIO::FTHelper->new;
       $fth->key('repeat_region');
    
    # Add location line
    my $g_start = $rep->start;
    my $g_end   = $rep->end;
    my $loc     = "$g_start..$g_end";

    if ($rep->strand == -1) {
        $loc = "complement($loc)";
    }

    $fth->loc($loc);
    
    # Add note describing repeat
    my $type    = $rep->hseqname;
       $type =~ s/^Motif://;
    my $r_start = $rep->hstart;
    my $r_end   = $rep->hend;

    $fth->add_field('note', "$type: matches $r_start to $r_end of consensus");
    
    return $fth;
}


1;

__END__
