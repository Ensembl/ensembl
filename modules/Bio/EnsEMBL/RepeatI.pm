
#
# BioPerl module for Bio::EnsEMBL::RepeatI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::RepeatI - Interface definition of a repeat class for Ensembl

=head1 SYNOPSIS

   # this class is only used to provide a bridge between 
   # C and Perl extensions

=head1 DESCRIPTION

This class is inherieted by the Perl and C extensions for Repeats.
It doesn't actually even define any more abstract methods over and
above the Bio::EnsEMBL::SeqFeatureI

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::RepeatI;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::SeqFeatureI;

# Object preamble - inherits from Bio::Root::Object

@ISA = qw(Bio::EnsEMBL::SeqFeatureI);


# no code

1;
