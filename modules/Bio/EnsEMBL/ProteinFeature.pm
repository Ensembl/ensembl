
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

ProteinFeature - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::ProteinFeature;
use vars qw($AUTOLOAD @ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::SeqFeature);

sub new {
    my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;
    
    # set stuff in self from @args
    return $self; # success - we hope!
}

=head2 translation

 Title   : translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub translation{
    my ($self,$value) = @_;
    
    if(defined($value)) {
	#Check if translation exists!
	$self->{'_translation'} = $value;
    } 

    return $self->{'_translation'};
}

