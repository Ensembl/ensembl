
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
use Bio::EnsEMBL::Analysis::LegacyParser;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ObjI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$clone,@args) = @_;

  # DEBUG
  # second parameter is for debugging to avoid reading entire list of objects
  if($clone){
      $self->warn("DEBUG: only exon/transcript/gene objects associated with $clone are read");
  }

  my $make = $self->SUPER::_initialize;

  $self->{'_gene_hash'} = {};
  $self->{'_contig_hash'} = {};

  # set stuff in self from @args
  # (nothing)

  # unfinished analysis is in humpub, so use humconf to get location
  #use humConf qw(HUMPUB_ROOT);
  # NO - should hard code it here or provide it as a paramter. Not a humConf thing...

  my $HUMPUB_ROOT = '/nfs/disk100/humpub/';
  
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
  # (better to do it here once than each time we need the information!)
  # FIXME - should this be moved to the pipeline so that this information
  # is stored in DBM files - currently in legacy parser
  my $p=Bio::EnsEMBL::Analysis::LegacyParser->new($gene_file,$transcript_file,$exon_file);
  $p->map_all($self,$clone);

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
    $self->{'_gene_hash'}->{$geneid} || 
	$self->throw("No gene with $geneid stored in TimDB");
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

    $self->{'_contig_hash'}->{$contigid} || 
	$self->throw("No contig with $contigid stored in this in-memory TimDB");
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

# close the dbm clone file

sub DESTROY{
    my ($obj) = @_;
    if( $obj->{'_clone_dbm'} ) {
	dbmclose(%{$obj->{'_clone_dbm'}});
	$obj->{'_clone_dbm'} = undef;
    }
}



