# BioPerl module for Bio::EnsEMBL::Chromosome
#
# Creator: Arne Stabenau <stabenau@ebi.ac.uk>
# Date of creation: 07.04.2000
# Last modified : 09.04.2000 by Arne Stabenau
#
# Copyright EMBL-EBI 2000
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Chromosome

=head1 SYNOPSIS


=head1 DESCRIPTION

    Contains very basic information of a chromosome + a function call to
    get the contigs which belong to the chromosome. A special unknown
    chromosome will hold contigs where the chromosome information is not
    available.

=head1 CONTACT


    Contact Arne Stabenau on implemetation/design detail: stabenau@ebi.ac.uk
    Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _
    
=cut


# Let the code begin...


package Bio::EnsEMBL::Chromosome;
use vars qw(@ISA);
use strict;
use Bio::Root::Object;
require Bio::EnsEMBL::Species;



@ISA = qw( Bio::Root::Object );

# new Chromosome( species-nickname, chromosome-name, chromosome-id
# dont use it, try get_by_name to get one.

sub _initialize {
    my $self = shift;
    my $nickname = shift;
    my $chrName = shift;
    my $dbId = shift;

    $self->SUPER::_initialize( @_ );
    $self->{_nickname} = $nickname;
    $self->{_chrname} = $chrName;
    $self->{_dbId} = $dbId;
    $self;
}

=head2 name
    
    Title   : name
    Usage   : $name = $chromomsome->name;

    Function: get the name of the chromosome. 
    Example : -

    Returns : -
    Args    : -
    
=cut


sub name {
    my $self = shift;
    $self->{_chrname};
}

=head2 get_by_name
    
    Title   : name
    Usage   : $chromsome = $chromomsome->get_by_name( $nickname, $chrname );

    Function: gets a chromosome object by nickname of species and its name.
              Use 'unknown' as chrname for contigs with not known
              chromosome.
    Example : -

    Returns : a Chromsome or undef, when the name is not known.
    Args    : -
    
=cut


sub get_by_name {
    # first parameter is ignored anyway
    shift;
    my $nickname = shift;
    my $name = shift;
    my $chr = Bio::EnsEMBL::Species->chromosome_by_name
	( $nickname, $name );
    return $chr;
}

sub get_by_id {
    # first parameter is ignored anyway
    shift;
    my $id = shift;
    my $chr = Bio::EnsEMBL::Species->chromosome_by_id
	( $id );
    return $chr;
}

=head2 species
    
    Title   : species
    Usage   : $species = $chromosome->species;

    Function: get the species object for the chromosome
    Example : -

    Returns : -
    Args    : -
    
=cut


sub species {
    my $self = shift;
  Bio::EnsEMBL::Species->new( $self->{_nickname} );
}

=head2 contigs
    
    Title   : contigs
    Usage   : @contigs = $chromsome->contigs( $dbobj );
    Function: get the contigs belonging to this chromosome from the give obj. 
    Example : -
    Returns : a list of names.
    Args    : -
    
=cut


sub contigs {
    my $self = shift;
    my $dbobj = shift;
    # please implement this in $obj
    my @contigs;
    eval {
	@contigs = $dbobj->get_Contigs_by_Chromosome( $self );
    };
    if( $@ ) {
	$self->throw( "Function contig not implemented in ".ref($self)."!" );
    }
    return @contigs;
}

=head2 get_db_id
    
    Title   : get_db_id
    Usage   : $dbid = $chr->get_db_id
    Function: returns the id of this chromosome for use in the SQL db.
    Example : -
    Returns : a number
    Args    : -
    
=cut

sub get_db_id {
    my $self = shift;
    $self->{_dbId};
}

# compiled successfull

1;




