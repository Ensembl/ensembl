
#
# BioPerl module for Bio::EnsEMBL::DB::CloneI
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::CloneI - Abstract Interface of a clone object

=head1 SYNOPSIS

    # get a clone object somehow

    @contigs = $clone->get_all_Contigs();

    @genes   = $clone->get_all_Genes();

    # dumping EMBL format

    $ostream = Bio::AnnSeqIO->new( -format => 'EMBL' , -fh => \*STDOUT );
    $annseq  = $clone->get_AnnSeq();
    $ostream->write_annseq($annseq);

    
=head1 DESCRIPTION

Ewan Birney

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::CloneI;
use vars qw($AUTOLOAD @ISA $CONTIG_SPACING);

use strict;
use POSIX;
#use Bio::EnsEMBL::DB::EmblVirtualContig;
use Bio::EnsEMBL::Virtual::EmblClone;


# Object preamble - inheriets from Bio::Root::Object

$CONTIG_SPACING = 800;

@ISA = qw();
# new() is inherited from Bio::Root::Object


=head2 id

 Title   : id
 Usage   : this is the primary id for ensembl. General the accession number from embl.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented this yet!");

}

=head2 embl_id

 Title   : embl_id
 Usage   : this is the embl_id for this clone, to generate nice looking files
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_id {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented embl_id yet!");

}


=head2 sv

 Title   : sv
 Function: returns the version number (not the acc.version, just verision).
 Example :
 Returns : 
 Args    :


=cut

sub sv {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented this yet!");

}


=head2 embl_version

 Title   : embl_version
 Usage   : $clone->embl_version()
 Function: Gives the value of the EMBL version, i.e. the data version
 Example : $clone->embl_version()
 Returns : version number
 Args    : none


=cut

sub embl_version {
    my ($self,@args) = @_;

    $self->warn("Base class has not implemented embl_version yet!");

}


=head2 seq_date

 Title   : seq_date
 Usage   : $clone->seq_date()
 Function: loops over all $contig->seq_date, throws a warning if they are different and 
           returns the first unix time value of the dna created datetime field, which indicates
           the original time of the dna sequence data
 Example : $clone->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date {
    my ($self,@args) = @_;

    $self->warn("Base class has not implemented seq_date yet!");
}


=head2 version

 Title   : version
 Function: Schema translation
 Example :
 Returns : 
 Args    :


=cut

sub version {
   my ($self,@args) = @_;

   $self->warn("Called version without implementation. Probably an old object. Called sv instead");

   return $self->sv();
}


=head2 htg_phase

 Title   : htg_phase
 Usage   : this is the phase being 1,2,3,4 (4 being finished).
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub htg_phase {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented htg_phase yet!");

}


=head2 created

 Title   : created
 Usage   : $clone->created()
 Function: Gives the unix time value of the created datetime field, which indicates
           the first time this clone was put in ensembl
 Example : $clone->created()
 Returns : unix time
 Args    : none


=cut

sub created {
    my ($self,@args) = @_;
    
    $self->warn("Base class has not implemented created  yet!");
}


=head2 modified

 Title   : modified
 Usage   : $clone->modified()
 Function: Gives the unix time value of the modified datetime field, which indicates
           the last time this clone was modified in ensembl
 Example : $clone->modified()
 Returns : unix time
 Args    : none


=cut

sub modified{
    my ($self,@args) = @_;

    $self->warn("Base class has not implemented modified yet!");
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig {
   my ($self,@args) = @_;

   $self->warn("Base class has not implemented get_Contig yet!");

}


=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   : 
 Function: gets all the contigs in this clone
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs {
   my ($self) = @_;

   $self->warn("Base class has not implemented get_all_Contigs yet!");

}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function: gets all the genes that overlap with this clone
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes {
   my ($self) = @_;

   $self->warn("Base class has not implemented get_all_Genes yet!");
}

=head2 get_all_ContigOverlaps 

 Title   : get_all_ContigOverlaps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ContigOverlaps {
    my ($self) = @_;
    $self->throw("Base class has not implemented get_all_ContigOverlaps");
}

=head2

Decorating functions. You do not need to implement these functions

=cut

=head2 virtualcontig

 Title   : virtualcontig
 Usage   : $vc = $clone->virtualcontig();
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub virtualcontig {
   my ($self) = @_;
   
   my $vc = Bio::EnsEMBL::Virtual::EmblClone->new($self);
   $vc->id($self->id);
   $vc->sv($self->embl_version);
   
   my $created = $self->_set_embl_date_format($self->created)." (CREATION DATE)";
   my $modified= $self->_set_embl_date_format($self->modified)." (LAST MODIFICATION DATE)";
   
   $vc->add_date($created);
   $vc->add_date($modified);
   return $vc;
}


sub _set_embl_date_format {

    my ($self,$unixtime)=@_;

    my ($seconds,$minutes,$hours,$day,$month,$year)=localtime ($unixtime);
    
    $month=$month+1;    
    $year=$year+1900;
   
    my %months = (
		1=>'JAN' , 2=>'FEB', 3=>'MAR', 4=>'APR',
		5=>'MAY', 6=>'JUN', 7=>'JUL', 8=>'AUG',
		9=>'SEP' , 10=>'OCT', 11=>'NOV', 12=>'DEC'
		);
   
    $month=$months{$month};  
   
    my $embl_date=$day."-".$month."-".$year;

    return $embl_date;
}


1;




