
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

Bio::EnsEMBL::TimDB::Obj - Object representing Tim's directory structure

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

package Bio::EnsEMBL::TimDB::Obj;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DB::ObjI;
use Bio::EnsEMBL::TimDB::Clone;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ObjI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_gene_hash'} = {};
  $self->{'_contig_hash'} = {};

  # set stuff in self from @args
  # (nothing)

  # unfinished analysis is in humpub, so use humconf to get location
  use humConf qw(HUMPUB_ROOT);
  
  # in order to access the flat file db, check that we can see the master dbm file
  # that will tell us where the relevant directory is
  my $unfinished_root="$HUMPUB_ROOT/th/unfinished_ana";
  $self->{'_unfinished_root'}=$unfinished_root;
  my $clone_dbm_file="$HUMPUB_ROOT/th/unfinished_ana/unfinished_clone.dbm";
  my %unfin_clone;
  unless(dbmopen(%unfin_clone,$clone_dbm_file,0666)){
      $self->throw("Error opening clone dbm file");
  }
  $self->{'_clone_dbm'}=\%unfin_clone;

  # define a few other important files
  my $exon_file="$HUMPUB_ROOT/blast/confirmed_exon";
  my $transcript_file="$HUMPUB_ROOT/th/unfinished_ana/unfinished_ana.transcript.lis";
  my $gene_file="$HUMPUB_ROOT/th/unfinished_ana/unfinished_ana.gene.lis";
  if(!-e $exon_file){
      $self->throw("Could not access exon file");
  }
  if(!-e $transcript_file){
      $self->throw("Could not access transcript file");
  }
  if(!-e $gene_file){
      $self->throw("Could not access gene file");
  }
  # only exon file needs to be saved as it contains more information than in following mappings
  $self->{'_exon_file'}=$exon_file;

  # build mappings from these flat files
  # FIXME - this should be moved to the pipeline so that this information
  # is stored in DBM files.
  # (better to do it here once than each time we need the information!)
  my %contig2exon;
  my %exons;
  my %exon2transcript;
  my %transcriptExons;
  my %transcript2gene;
  my %geneTranscripts;
  $self->{'_contig2exon'}=\%contig2exon;
  $self->{'_exons'}=\%exons;
  $self->{'_exon2transcript'}=\%exon2transcript;
  $self->{'_transcriptExons'}=\%transcriptExons;
  $self->{'_transcript2gene'}=\%transcript2gene;
  $self->{'_geneTranscripts'}=\%geneTranscripts;

  return $make; # success - we hope!
}

=head2 get_Gene

 Title   : get_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene{
   my ($self,$geneid) = @_;

   $self->throw("Tim has not reimplemented this function");

   $self->{'_gene_hash'}->{$geneid} || $self->throw("No gene with $geneid stored in this in-memory TimDB");
   return $self->{'_gene_hash'}->{$geneid};
}

=head2 get_Clone

 Title   : get_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Clone{
   my ($self,$id) = @_;

   # check to see if clone exists, and extract relevant items from dbm record
   # cgp is the clone category (SU, SF, EU, EF)
   my($line,$cdate,$type,$cgp,$acc,$sv);
   if($line=$self->{'_clone_dbm'}->{$id}){
       ($cdate,$type,$cgp,$acc,$sv)=split(/,/,$line);
   }else{
       $self->throw("$id is not a valid sequence in this database");
   }

   # create clone object
   my $clone = new Bio::EnsEMBL::TimDB::Clone(-id => $id,
					      -dbobj => $self,
					      -cgp => $cgp);
   return $clone;
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
    my ($self,$contigid)= @_;

    $self->throw("Tim has not reimplemented this function");

    $self->{'_contig_hash'}->{$contigid} || $self->throw("No contig with $contigid stored in this in-memory TimDB");
    return $self->{'_contig_hash'}->{$contigid};
}


=head2 write_Gene

 Title   : write_Gene
 Usage   : $obj->write_Gene($gene)
 Function: writes a particular gene into the database
           
 Example :
 Returns : 
 Args    :


=cut

sub write_Gene{
   my ($self,$gene) = @_;
   $self->throw("Cannot write to a TimDB");
}

=head2 write_Contig

 Title   : write_Contig
 Usage   : $obj->write_Contig($contigid,$dna)
 Function: writes a contig and its dna into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Contig {
   my ($self,$contig) = @_;
   $self->throw("Cannot write to a TimDB");
}

# simple internal methods (hardwired things that would have been done by AUTOLOAD)
#sub _db_handle{
#   my($self,$value)=@_;
#   $self->{'_db_handle'}=$value if(defined $value);
#   return $self->{'_db_handle'};
#}

#sub AUTOLOAD {
#    my $self = shift;
#    my $class = ref($self) || $self->throw("\'$self\' is not an object of mine!");
#    my $name = $AUTOLOAD;
#
#    # don't propagate DESTROY messages...
#
#    $name =~ /::DESTROY/ && return;
#
#    unless (exists $self->{'_permitted'}->{$name}) {
#	$self->throw("In type $class, can't access $name - probably passed a wrong variable");
#    }
#    if (@_) {
#	return $self->{$name} = shift;
#    } else {
#	return $self->{$name}
#    }
#}

# close the dbm clone file
sub DESTROY{
   my ($obj) = @_;
   if( $obj->{'_clone_dbm'} ) {
       dbmclose(%{$obj->{'_clone_dbm'}});
       $obj->{'_clone_dbm'} = undef;
   }
}



