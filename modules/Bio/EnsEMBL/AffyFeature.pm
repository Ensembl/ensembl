#
# Ensembl module for Bio::EnsEMBL::AffyFeature
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AffyFeature - a module to represent a affy probe hitting the genomic sequence.


=head1 SYNOPSIS

use Bio::EnsEMBL::AffyFeature;

$feature = Bio::EnsEMBL::AffyArray->new 
    ( -probe => $probe,
      -slice => $slice,
      -start => 23,
      -end => 30,
      -strand => 1,
      -probeset => 'some setname'
)

=head1 DESCRIPTION

AffyFeatures represent an oligo (probe) on an AffyArray that matches the genome. 
Its possible to get a probe for them, for performance reasons the probeset 
(which is the more important information) is settable and
retrievable without the probe information.

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::AffyFeature;

use Bio::EnsEMBL::Feature;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw );



use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);

sub new {
  my $caller = shift;

  #allow constructor to be called as class or object method
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my($probe, $probeset, $mismatchcount ) =
      rearrange(['PROBE', 'MISMATCHCOUNT' ],
		@_);

  $self->{'probe'} = $probe;
  $self->{'mismatchcount'} = $mismatchcount;

  return $self;
}

sub new_fast {
  my ( $class, $hashref )  = @_;

  return bless( $hashref, $class );
}



# if the probe is set, take it from there
sub probeset {
    my $self = shift;
    $self->{'probeset'} = shift if( @_ );
    if( $self->{'probe'} ) {
	$self->{'probeset'} = $self->probe()->probeset();
    }
    return $self->{'probeset'};
}

sub mismatchcount {
    my $self = shift;
    $self->{'mismatchcount'} = shift if @_;
    return $self->{'mismatchcount'};
}

sub probelength {
    my $self = shift;
    $self->length();
}



# lazy load this
sub probe {
    my $self = shift;
    if( @_ ) {
	my $probe = shift;
	if( $probe ) {
	    throw( "Need Bio::EnsEMBL::AffyProbe as probe" ) unless
		(ref( $probe ) && $probe->isa( "Bio::EnsEMBL::AffyProbe" ));
	}
	$self->{'probe'} = $probe;
    } else {
	if( !defined $self->{'probe'} && 
	    defined $self->adaptor() && 
	    $self->dbID() ) {
	    $self->{'probe'} = $self->adaptor()->db()->
		get_AffyProbeAdaptor()->fetch_by_AffyFeature( $self );
	}
    }
    return $self->{'probe'};
}







