#
# BioPerl module for Bio::EnsEMBL::TimDB::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::SlimTim::Obj - Object representing Tims directory structure

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::SlimTim::Obj;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DB::ObjI;
use Bio::EnsEMBL::SlimTim::Clone;

use Bio::EnsEMBL::Analysis::LegacyParser;
use Bio::EnsEMBL::Analysis::ensConf qw(UNFIN_ROOT
				       UNFIN_DATA_ROOT
				       UNFIN_DATA_ROOT_CGP
				       CONFIRMED_EXON_FASTA
				       );
use NDBM_File;
use Fcntl qw( O_RDONLY );

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ObjI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args)=@_;

    my $make = $self->SUPER::_initialize(@args);
    my ($dir) = $self->_rearrange([qw(TIMDIR)],@args);
    
    $self->_dir($dir);
    
    $self->{'_species'}='human';
    $self->{'_contig_hash'} = {};
    
    
    return $make; # success - we hope!
}


=head2 get_Clone

 Title   : get_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :
    
=cut
    
sub get_Clone {
    my ($self,$id) = @_;
    
    # create clone object
    my $clone = new Bio::EnsEMBL::SlimTim::Clone(-diskid => $id,
					       -dbobj => $self);
    
    return $clone;
}

=head2 _dir

 Title   : _dir
 Usage   : $obj->_dir($newval)
 Function: 
 Returns : value of _dir
 Args    : newvalue (optional)


=cut

sub _dir{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_dir'} = $value;
    }
    return $obj->{'_dir'};

}

1;









