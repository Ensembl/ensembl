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


=head2 new

  Arg [PROBE] : 
               Every AffyFeature needs an AffyProbe on construction. This should be already
               stored if you plan to store this feature.
  Arg [MISMATCHCOUNT] :
               How many mismatches over the length of the probe? (0,1) 
  Example    : my $feature = Bio::EnsEMBL::AffyFeature->new(
                 -PROBE => $affyProbe,
                 -MISMATCHCOUNT => 0,
                 -SLICE => $chr_1_slice,
                 -START => 1_000_000,
                 -END => 1_000_024,
                 -STRAND => -1 ); 
  Description: Constructor for AffyFeature objects. Mainly based on the Feature 
               constructor.
  Returntype : Bio::EnsEMBL::AffyFeature
  Exceptions : none
  Caller     : general


=cut


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


=head2 new_fast

  Args       : hashref with all internal attributes set
  Example    : none
  Description: Fast and dirty version of new. Only works because the code is 
               very diciplined :-)
  Returntype : Bio::EnsEMBL::AffyFeature
  Exceptions : none
  Caller     : general


=cut


sub new_fast {
  my ( $class, $hashref )  = @_;

  return bless( $hashref, $class );
}



=head2 probeset

  Args       : none 
  Example    : my $probes = $array->get_all_AffyProbes();
  Description: The probeset for this feature. Shortcut for $feature->probe->probeset().
               Possibly not needed.
  Returntype : string
  Exceptions : none
  Caller     : general


=cut

sub probeset {
    my $self = shift;
    $self->{'probeset'} = shift if( @_ );
    if( $self->{'probe'} ) {
	$self->{'probeset'} = $self->probe()->probeset();
    }
    return $self->{'probeset'};
}



=head2 mismatchcount

  Args       : int $mismatchcount
  Example    : none
  Description: getter setter for attribute mismatchcount
  Returntype : int
  Exceptions : none
  Caller     : general


=cut

sub mismatchcount {
    my $self = shift;
    $self->{'mismatchcount'} = shift if @_;
    return $self->{'mismatchcount'};
}


=head2 probelength

  Args       : none 
  Example    : none
  Description: The probelength is supposed to be the length of the match.
               Doesnt really belong here, but iss kind of too much to store 
               for each probe (probably 25 anyway)..
  Returntype : int
  Exceptions : none
  Caller     : general


=cut


sub probelength {
    my $self = shift;
    $self->length();
}



=head2 probe

  Args       : Bio::EnsEMBL::AffyProbe $probe
  Example    : none
  Description: getter / setter / lazy load of the probe attribute. Feature are retrieved
               from the database without attached probes. They are lazy loaded on demand.
               That means that retrieving probe information from Feature creates an SQL-query.
  Returntype : Bio::EnsEMBL::AffyProbe
  Exceptions : none
  Caller     : general


=cut

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







