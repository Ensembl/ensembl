
#
# BioPerl module for Bio::EnsEMBL::DB::EmblVirtualContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::EmblVirtualContig - Virtual Contig for EMBL dumping

=head1 SYNOPSIS

    # EmblVirtualContigs are usually made from Clones.
    # if you want to Dump EMBL/GenBank format, you might
    # want to make EmblVirtualContigs 'by hand'. Check out the
    # CloneI code.

    $clone = $dbobj->get_Clone('AC000056');
    $vc = $clone->virtualcontig();

    # is-a ContigI
    foreach $gene ( $vc->get_all_Genes() ) {
       # do something
    }
    
    # can dump nicely to EMBL flat files
    $seqio = Bio::SeqIO->new( '-format' => 'embl' , -fh => \*STDOUT );
    $seqio->write_seq($vc);
    

=head1 DESCRIPTION

EmblVirtualContigs are a type of VirtualContig which has accessor for
all the annoying details for EMBL flat file dumping, such as dates,
divisions, phases, keywords and other aspects. We've put this into a
derived object so that the main VC code is not cluttered up with 
trivialities 

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
 to one of the Bioperl mailing lists.
Your participation is much appreciated.

  vsns-bcd-perl@lists.uni-bielefeld.de          - General discussion
  vsns-bcd-perl-guts@lists.uni-bielefeld.de     - Technically-oriented discussion
  http://bio.perl.org/MailList.html             - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via email or the web:

  bioperl-bugs@bio.perl.org
  http://bio.perl.org/bioperl-bugs/

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::EmblVirtualContig;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DB::VirtualContig;
use Bio::Annotation;


@ISA = qw(Bio::EnsEMBL::DB::VirtualContig);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  $self->{'date'} = [];
  $self->annotation( Bio::Annotation->new());
  $self->{'additional_seqf'} = [];

 # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 desc

 Title   : desc
 Usage   : $obj->desc($newval)
 Function: 
 Example : 
 Returns : value of desc
 Args    : newvalue (optional)


=cut

sub desc{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'desc'} = $value;
    }
    return $obj->{'desc'};

}

=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: 
 Example : 
 Returns : value of species
 Args    : newvalue (optional)


=cut

sub species{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'species'} = $value;
    }
    return $obj->{'species'};

}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}

=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->htg_phase($newval)
 Function: 
 Example : 
 Returns : value of htg_phase
 Args    : newvalue (optional)


=cut

sub htg_phase{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};

}

=head2 sv

 Title   : sv
 Usage   : $obj->sv($newval)
 Function: 
 Example : 
 Returns : value of sv
 Args    : newvalue (optional)


=cut

sub sv{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'sv'} = $value;
    }
    return $obj->{'sv'};

}


=head2 embl_id

 Title   : embl_id
 Usage   : $obj->embl_id($newval)
 Function: 
 Example : 
 Returns : value of embl_id
 Args    : newvalue (optional)


=cut

sub embl_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}

=head2 project_name

 Title   : project_name
 Usage   : $obj->project_name($newval)
 Function: 
 Example : 
 Returns : value of project_name
 Args    : newvalue (optional)


=cut

sub project_name{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'project_name'} = $value;
    }
    return $obj->{'project_name'};

}

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: 
 Example : 
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}

=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_SeqFeature{
   my ($self,$sf) = @_;

   push(@{$self->{'additional_seqf'}},$sf);
}


=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures{
   my ($self) = @_;
   my @sf;

   @sf = $self->SUPER::top_SeqFeatures();
   push(@sf,@{$self->{'additional_seqf'}});
   return @sf;
}
