# BioPerl module for Bio::EnsEMBL::Species
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 06.04.2000
# Last modified : 07.04.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Species

=head1 SYNOPSIS


=head1 DESCRIPTION

    A simple species object with 2 hardcoded species (mouse and
    human). Chromosome creation is also done here, to localize the
    hardcoding to this file.

=head1 CONTACT


    Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
    Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _
    
=cut


# Let the code begin...


package Bio::EnsEMBL::Species;
use vars qw(@ISA);
use strict;
use Bio::Root::Object;
use Bio::EnsEMBL::Chromosome;

@ISA = qw( Bio::Root::Object );

my @SpeciesList = ( [
		     'human', '9606', 
		     [ '1', '2', '3', '4', '5', '6', '7', '8', '9',
		       '10','11','12','13','14','15','16','17','18',
		       '19','20','21','22','X', 'Y','unknown' ],
		     [ 1..25 ]
		     ],
		    ['mouse','10090',
		     [ '1', '2', '3', '4', '5', '6', '7', '8', '9',
		       '10','11','12','13','14','15','16','17','18',
		       '19','X', 'Y','unknown' ],
		     [ 26..47 ] 
		     ] );

# numbers in species list are hardcoded EnsEMBL-ids for the chromosomes.
# so unknowns are 25 and 47 for homo and mus respectivly

sub _initialize {
    my $self = shift;
    my $nickname = shift;
    
    $self->SUPER::_initialize( @_ );
    foreach my $species (@SpeciesList) {
	if(( $species->[0] cmp $nickname ) == 0 ) {
	    $self->{_list} = $species;
	    return $self;
	}
    }
    $self->throw( "Unknown nickname in species constructor!" );
}
    

=head2 nickname
    
    Title   : nickname
    Usage   : $nickname = $species->nickname;

    Function: get the nickname of the species. 
    Example : -

    Returns : -
    Args    : -
    
=cut


sub nickname {
    my $self = shift;
    $self->{_list}->[0];
}

=head2 taxonomy_id
    
    Title   : taxonomy_id
    Usage   : $tid = $species->taxonomy_id;

    Function: get the NCBI taxonomy id of the species. 
    Example : -

    Returns : -
    Args    : -
    
=cut


sub taxonomy_id {
    my $self = shift;
    $self->{_list}->[1];
}

=head2 chromosomes
    
    Title   : chromosomes
    Usage   : @chromosomes = $species->chromosomes
    Function: get the chromosome-objects of the species. 
    Example : -

    Returns : a list of objects
    Args    : -
    
=cut


sub chromomsomes {
    my $self = shift;
    my @result = ();
    for ( 0..$#{$self->{_list}->[2]}) {
	my $chr = Bio::EnsEMBL::Chromosome->new
	    ( $self->{_list}->[0],
	      $self->{_list}->[2]->[$_],
	      $self->{_list}->[3]->[$_] );
	push( @result, $chr );
    }
    return @result;
}

=head2 chromosome_by_id
    
    Title   : chromosome_by_id
    Usage   : $chr = $species->chromosome_by_id( $id )
    Function: get a chromosome-object. The id is the one you should as well
              use when stroing a contig in SQL.
    Example : -

    Returns : a chromosome object or undef if id is not valid.
    Args    : -
    
=cut

sub chromosome_by_id {
    my $self = shift;
    my $id = shift;
    
    for my $species (@SpeciesList) {
	for my $idx (0.. $#{$species->[2]}) {
	    if( $species->[3]->[$idx] == $id ) {
		my $chr = Bio::EnsEMBL::Chromosome->new
		    ( $species->[0],
		      $species->[2]->[$idx],
		      $species->[3]->[$idx] );
		return $chr;
	    }
	}
    }
    return undef;
}
=head2 chromosome_by_name
    
    Title   : chromosome_by_name
    Usage   : $chr = $species->chromosome_by_name( $nickname, $chrname )
    Function: get a chromosome-object by species and chromosome name.
    Example : -

    Returns : a chromosome object or undef if names are not valid.
    Args    : -
    
=cut

sub chromosome_by_name {
    my $self = shift;
    my $nickname = shift;
    my $chrname = shift;

    for my $species (@SpeciesList) {
	next, if(( $species->[0] cmp $nickname ) != 0 );
	for my $idx (0..$#{$species->[3]}) {
	    if(( $species->[2]->[$idx] cmp $chrname ) == 0 ) {
		my $chr = Bio::EnsEMBL::Chromosome->new
		    ( $species->[0],
		      $species->[2]->[$idx],
		      $species->[3]->[$idx] );
		return $chr;
	    }
	}
    }
    return undef;
}



# compiled successfull

1;
