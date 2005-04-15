#
# Ensembl module for Bio::EnsEMBL::Mapper::IndelPair
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::IndelPair

=head1 SYNOPSIS

=head1 DESCRIPTION

Two regions mapped between different coordinate systems are
each represented by a Bio::EnsEMBL::Mapper::Unit and joined
together as a Bio::EnsEMBL::Mapper::IndelPair, when one of the 
regions is an indel.

=head1 AUTHOR - Ewan Birney

This module is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::IndelPair;
use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Mapper::Pair);

=head2 new

    Arg[1]      : int $from. From where the Indel starts
    Arg[2]      : int $to. Where the Indel ends
    Arg[3]      : int $orientation (1/-1)
    Example     : my $IndelPair = Bio::EnsEMBL::Mapper::IndelPair->new(1000,1005,1)
    Description : Creator of the object. Creates a new Pair object, but adding 
                  the flag 'indel' to differentiate from the Pair object
    ReturnType  : Bio::EnsEMBL::Mapper::IndelPair
    Exceptions  : none
    Caller      : general
                  
=cut

sub new {
    my $caller = shift;
    my $class = ref($caller) || $caller;

    my $self = $class->SUPER::new(@_); #create the Pair object
    $self->{'indel'} = 1; #and add the Indel flag
    return $self;
}

1;
