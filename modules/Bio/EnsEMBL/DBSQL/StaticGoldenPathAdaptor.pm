
#
# Ensembl module for Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor - Database adaptor for static golden path

=head1 SYNOPSIS

    # get a static golden path adaptor from the obj

    $adaptor = $dbobj->get_StaticGoldenPathAdaptor();

    # these return sorted lists:

    @rawcontigs = $adaptor->fetch_RawContigs_by_fpc_name('ctg123');

    @rawcontigs = $adaptor->fetch_RawContigs_by_chr('chr2');

    # can throw an exception: Not on Same Chromosome
    @rawcontigs = $adaptor->fetch_RawContigs_between_RawContigs($start_rc,$end_rc);

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Virtual::StaticContig;

@ISA = qw(Bio::Root::RootI);

# new() is written here 

sub new {
  my($class,@args) = @_;
  
  my $self = {};
  bless $self,$class;
  
  my ($dbobj) = $self->_rearrange([qw( DBOBJ)],@args);

  if( !defined $dbobj) {
      $self->throw("got no dbobj. Aaaaah!");
  }

  $self->dbobj($dbobj);

# set stuff in self from @args
  return $self;
}

=head2 fetch_RawContigs_by_fpc_name

 Title   : fetch_RawContigs_by_fpc_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_RawContigs_by_fpc_name{
   my ($self,$fpc) = @_;
   
   my $type = $self->dbobj->static_golden_path_type();
   
   # very annoying. DB obj wont make contigs by internalid. doh!
   my $sth = $self->dbobj->prepare("select c.id from static_golden_path st,contig c where c.internal_id = st.raw_id AND st.fpcctg_name = '$fpc' AND  st.type = '$type' ORDER BY st.fpcctg_start");
   $sth->execute;
   my @out;
   my $cid;
   while( ( my $cid = $sth->fetchrow_arrayref) ) {
       my $rc = $self->dbobj->get_Contig($cid->[0]);
       push(@out,$rc);
   }

   return @out;
}

=head2 fetch_RawContigs_by_chr_name

 Title   : fetch_RawContigs_by_chr_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_RawContigs_by_chr_name{
   my ($self,$chr) = @_;

   my $type = $self->dbobj->static_golden_path_type();
   
   # very annoying. DB obj wont make contigs by internalid. doh!
   my $sth = $self->dbobj->prepare("select c.id from static_golden_path st,contig c where c.internal_id = st.raw_id AND st.chr_name = '$chr' AND  st.type = '$type' ORDER BY st.fpcctg_start");
   $sth->execute;
   my @out;
   my $cid;
   while( ( my $cid = $sth->fetchrow_arrayref) ) {
       my $rc = $self->dbobj->get_Contig($cid->[0]);
       push(@out,$rc);
   }

   return @out;
}


=head2 VirtualContig_by_fpc_name

 Title   : VirtualContig_by_fpc_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub VirtualContig_by_fpc_name{
   my ($self,$name) = @_;
   
   return Bio::EnsEMBL::Virtual::StaticContig->new($self->fetch_RawContigs_by_fpc_name($name));
}

=head2 VirtualContig_by_chr_name

 Title   : VirtualContig_by_chr_name
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub VirtualContig_by_chr_name{
   my ($self,$name) = @_;

   return Bio::EnsEMBL::Virtual::StaticContig->new($self->fetch_RawContigs_by_chr_name($name));


}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}
