
=head1 NAME - Bio::EnsEMBL::Intron

=head1 SYNOPSIS

    # Given @exons, and array of Bio::EnsEMBL::Exon
    # objects, you can do:
    my $intron = Bio::EnsEMBL::Intron->new;
    $intron->upstream_Exon  (@exons[0]);
    $intron->downstream_Exon(@exons[1]);
    
    my $start  = $intron->start;
    my $end    = $intron->end;
    my $length = $intron->length;
    my $seq    = $intron->seq;

=head1 DESCRIPTION

Bio::EnsEMBL::Intron inherits most of its methods
from Bio::SeqFeature::Generic.

It requires two Bio::EnsEMBL::Exon objects, which
are the upstream and downsteam exons in the
transcript.  These exons must be attached to the
same Contig (or VirtualContig) object.

It overrides the location method in
Bio::SeqFeature::Generic, so that it can
automatically calculate its start and end.

It is not possible to set the start or end
directly.  It remains to be seen whether this is
a feature or a bug.

=over 4

=cut


package Bio::EnsEMBL::Intron;

use strict;
use Bio::SeqFeature::Generic;
use vars '@ISA';

@ISA = ('Bio::SeqFeature::Generic');

=head2 new

    my $intron = Bio::EnsEMBL::Intron->new;

Returns a new B<Bio::EnsEMBL::Intron> object. 
Does not take any arguments.  Sets the
C<primary_tag> to B<intron> and the C<source_tag>
to B<EnsEMBL>.

=cut

sub new {
    my $pkg = shift;
    
    #my $self = $pkg->SUPER::new;
    my $self = bless {}, $pkg;
    
    $self->throw("new() doesn't take any arguments") if @_;
    
    $self->primary_tag('intron');
    $self->source_tag('EnsEMBL');
    
    return $self;
}

=head2 upstream_Exon

    $intron->upstream_Exon($exon1);
    my $exon = $intron->upstream_Exon;

Set or get method for the exon which defines the
upstream border of this intron in the transcript.

=head2 downstream_Exon

    $intron->downstream_Exon($exon2);
    my $exon = $intron->downstream_Exon;

Set or get method for the exon which defines the
downstream border of this intron in the transcript.

=cut

sub upstream_Exon {
    my( $self, $exon ) = @_;
    
    if ($exon) {
        $self->{'_intron_location'} = undef;
        $self->throw("'$exon' is not a Bio::EnsEMBL::Exon")
            unless $exon->isa('Bio::EnsEMBL::Exon');
        $self->{'_upstream_exon'} = $exon;
    }
    return $self->{'_upstream_exon'};
}

sub downstream_Exon {
    my( $self, $exon ) = @_;
    
    if ($exon) {
        $self->{'_intron_location'} = undef;
        $self->throw("'$exon' is not a Bio::EnsEMBL::Exon")
            unless $exon->isa('Bio::EnsEMBL::Exon');
        $self->{'_downstream_exon'} = $exon;
    }
    return $self->{'_downstream_exon'};
}

sub location {
    my( $self ) = @_;
    
    unless ($self->{'_intron_location'}) {
        my $loc = Bio::Location::Simple->new;
    
        # Get the upstream and downstream exons
        my $up_exon = $self->upstream_Exon
            or $self->throw("Missing upstream_Exon");
        my $down_exon = $self->downstream_Exon
            or $self->throw("Missing downstream_Exon");
        
        # Get the PrimarySeqs attached to both and check it is the same sequence
        my $up_seq   = $up_exon  ->entire_seq;
        my $down_seq = $down_exon->entire_seq;
        unless ($up_seq == $down_seq) {
            $self->throw("upstream and downstream exons are attached to different sequences\n"
                . "'$up_seq' and '$down_seq'");
        }
        
        # Check that the exons are on the same strand.  (Do I need to bother?)
        my $up_strand   = $up_exon  ->strand;
        my $down_strand = $down_exon->strand;
        unless ($up_strand == $down_strand) {
            $self->throw("upstream and downstream exons are on different strands "
                . "('$up_strand' and '$down_strand')");
        }
        $loc->strand($up_strand);
        
        #   $exon_end is the  end  of the exon which is 5' of the intron on the genomic sequence.
        # $exon_start is the start of the exon which is 3' of the intron on the genomic sequence.
        my( $exon_end, $exon_start );
        if ($up_strand == 1) {
            $exon_end   = $up_exon  ->end;
            $exon_start = $down_exon->start;
        } else {
            $exon_end   = $down_exon->end;
            $exon_start = $up_exon  ->start;
        }
        unless ($exon_end < $exon_start) {
            $self->throw("Intron gap begins after '$exon_end' and ends before '$exon_start'");
        }
        $loc->start($exon_end   + 1);
        $loc->end  ($exon_start - 1);
        
        # Attach the sequence and location objects to the intron
        $self->attach_seq($up_seq);
        $self->{'_intron_location'} = $loc;
        
    }
    return $self->{'_intron_location'};
}

1;

__END__

=back

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

