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

Bio::EnsEMBL::Analysis::Repeat

=head1 SYNOPSIS

=head1 DESCRIPTION

Extends the bioperl Bio::SeqFeature::Homol to store
genomic repeat features.  

Creation:

    my $rep = new Bio::EnsEMBL::Analysis::Repeat(-start => $start,
						 -end    => $end,
						 -strand => $strand,
						 -source => $source,
						 -primary=> $primary,
						 );

Manipulation:

    my $start = $rep->start;
    my $end   = $rep->end;

    my $repname  = $rep->homol_SeqFeature->seqname;
    my $repstart = $rep->homol_SeqFeature->start;
    my $repend   = $rep->homol_SeqFeature->end;


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::Repeat;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Homol

use Bio::AnnSeqIO::FTHelper;
use Bio::SeqFeature::Homol;


@ISA = qw(Bio::SeqFeature::Homol);

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
    my $loc = "$g_start..$g_end";
    if ($rep->strand == -1) {
        $loc = "complement($loc)";
    }
    $fth->loc($loc);
    
    # Add note describing repeat
    my $type    = $rep->homol_SeqFeature->seqname;
    $type =~ s/^Motif://;
    my $r_start = $rep->homol_SeqFeature->start;
    my $r_end   = $rep->homol_SeqFeature->end;
    $fth->add_field('note', "$type: matches $r_start to $r_end of consensus");
    
    return $fth;
}


1;

__END__
