
#
# Object for converting file extensions to msp types
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::MSPtype - converting file extensions to MSP types

=head1 SYNOPSIS


=head1 DESCRIPTION



=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::MSPType;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

use Bio::Root::Object;

# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called


# This is by no means perfect - the data should be static for a start
#                           Program     database   type    extension               type   dbver progver gff_source
my $msptype =  [['swir_p',  'blastp',  'swir',     'pep', '.blastp_swir.msptmp',   'msp' ,  1 ,1,'similarity'],
		['ce_p',    'tblastn', 'ce',       'dna', '.tblastn_ce.msptmp',    'msp' ,  1 ,1,'similarity'],
		['vert_p',  'tblastn', 'vert',     'dna', '.tblastn_vert.msptmp',  'msp' ,  1 ,1,'similarity'],
		['sh_p',    'tblastn', 'sh',       'dna', '.tblastn_sh.msptmp',    'msp' ,  1 ,1,'similarity'],
		['dbest_p', 'tblastn', 'dbest',    'dna', '.tblastn_dbest.msptmp', 'msp' ,  1 ,1,'similarity'],
		['pfam_p',  'hmmpfam', 'PfamFrag', 'pep', '.hmmpfam_frag',         'pfam',  1 ,1,'pfam_prediction'],
		['repeat_n','RepeatMasker', ''   , 'dna', '.RepMask.out.gff',      'repeat','','042199','similarity'],
		['genewise_p','genewise', 'swir'   , 'pep', '.swir.msptmp',      'msp',1,1,'similarity'],
		];
  
sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;

  return $self; # success - we hope!
}


sub each_MSPType {
    my ($self) = @_;

    return @$msptype;
}


sub extension2MSPType {
    my ($self,$ext) = @_;

    foreach my $arr (each_MSPType) {
	if ($arr->[4] eq $ext) {
	    return $arr;
	}
    }

    print("Can't convert $ext\n");
    
}








