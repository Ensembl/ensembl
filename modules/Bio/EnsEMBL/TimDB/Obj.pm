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

Bio::EnsEMBL::TimDB::Obj - Object representing Tims directory structure

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
use Bio::EnsEMBL::Analysis::ensConf qw(UNFIN_ROOT
				       UNFIN_DATA_ROOT
				       CONFIRMED_EXON_FASTA
				       );

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ObjI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$raclones,$noacc,$test,$part,@args) = @_;

  # DEBUG
  # second parameter is for debugging to avoid reading entire list of objects
  if($raclones){
      $self->warn("DEBUG: only exon/transcript/gene objects associated with clone list are read");
  }

  my $make = $self->SUPER::_initialize;

  $self->{'_gene_hash'} = {};
  $self->{'_contig_hash'} = {};

  # clone->acc translation for all timdb operations, unless $noacc
  $self->{'_byacc'}=1 unless $noacc;

  # set stuff in self from @args
  # (nothing)

  # in order to access the flat file db, check that we can see the master dbm file
  # that will tell us where the relevant directory is
  # NOTE FIXME it is not very clever to have this open DBM file hanging, even if 
  # it is only for reading (cannot open readonly) since to certainly of locking
  # or dataconsistency
  my $unfinished_root="$UNFIN_ROOT";
  my $exon_file;
  if($test){
      $unfinished_root.="/test";
      $self->{'_test'}=1;
      $exon_file="$unfinished_root/test_confirmed_exon";
  }else{
      $exon_file="$CONFIRMED_EXON_FASTA";
  }
  $self->{'_unfinished_root'}=$unfinished_root;
  my $clone_dbm_file="$unfinished_root/unfinished_clone.dbm";
  my %unfin_clone;
  unless(dbmopen(%unfin_clone,$clone_dbm_file,0666)){
      $self->throw("Error opening clone dbm file");
  }
  $self->{'_clone_dbm'}=\%unfin_clone;

  # if going to do things !$noacc then need to open this dbm file too
  if(!$noacc){
      my $accession_dbm_file="$unfinished_root/unfinished_accession.dbm";
      my %unfin_accession;
      unless(dbmopen(%unfin_accession,$accession_dbm_file,0666)){
	  $self->throw("Error opening accession dbm file");
      }
      $self->{'_accession_dbm'}=\%unfin_accession;
  }

  # define a few other important files
  my $file_root;
  if($part){
      $file_root="$unfinished_root/partial";
  }else{
      $file_root="$unfinished_root";
  }
  my $transcript_file="$file_root/unfinished_ana.transcript.lis";
  my $gene_file="$file_root/unfinished_ana.gene.lis";
  my $contig_order_file="$file_root/unfinished_ana.contigorder.lis";
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
  my $p=Bio::EnsEMBL::Analysis::LegacyParser->new($gene_file,$transcript_file,
						  $exon_file,$contig_order_file);

  # doing conversion acc->id->acc or id->acc, need it here too
  $p->map_all($self,$raclones);

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

    my($disk_id,$cgp,$sv,$emblid,$htgsp);
    ($id,$disk_id,$cgp,$sv,$emblid,$htgsp)=$self->get_id_acc($id);

    # create clone object
    my $clone = new Bio::EnsEMBL::TimDB::Clone(-id => $id,
					       -disk_id => $disk_id,
					       -dbobj => $self,
					       -cgp => $cgp,
					       -sv=>$sv,
					       -emblid=>$emblid,
					       -htgsp=>$htgsp,
					       -byacc => $self->{'_byacc'},
					       );
    return $clone;
}


=head2 get_all_Clone_id

 Title   : get_all_Clone_id
 Usage   : @cloneid = $obj->get_all_Clone_id
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    :

Note: for speed this does not return the ensembl_id but the disk_id


=cut

sub get_all_Clone_id{
   my ($self,$fall) = @_;
   my($key,$val);
   my @list;
   my $nc=0;
   my $nisv=0;
   my $nsid=0;
   while(($key,$val)=each %{$self->{'_clone_dbm'}}){
       $nc++;
       my($cdate,$type,$cgp,$acc,$sv,$emblid,$htgsp)=split(/,/,$val);
       if($key ne $acc){
	   $nsid++;
       }
       if($sv!~/^[1234]$/){
	   $nisv++;
	   next unless $fall;
       }
       push(@list,$key);
   }
   print STDERR "$nc clones in database\n";
   print STDERR "$nsid have cloneid rather than accession numbers\n";
   print STDERR "$nisv have invalid SV numbers";
   if($fall){
       print STDERR " and are included\n";
   }else{
       print STDERR " and are excluded\n";
   }
   return @list;
}


=head2 get_id_acc

 Title   : get_id_acc
 Usage   : @array=$self->$id;
 Function: returns id (id or acc dependent on _byacc flag) and other parameters associated with a clone
 Example :
 Returns : 
 Args    :

=cut

sub get_id_acc{
    my($self,$id)=@_;
    # check to see if clone exists, and extract relevant items from dbm record
    # cgp is the clone category (SU, SF, EU, EF)
    my($line,$cdate,$type,$cgp,$acc,$sv,$id2,$fok,$emblid,$htgsp);
    if($line=$self->{'_clone_dbm'}->{$id}){
	# first straight forward lookup
	($cdate,$type,$cgp,$acc,$sv,$emblid,$htgsp)=split(/,/,$line);
	# translate to $acc if output requires this
	if($self->{'_byacc'}){
	    $id2=$id;
	    if(!$acc){
		$self->throw("Accession number is unknown for $id");
	    }
	    $id=$acc;
	}else{
	    $id2=$id;
	}
	$fok=1;
    }elsif(($self->{'_byacc'}) && ($id2=$self->{'_accession_dbm'}->{$id})){
	# lookup by accession number, if valid
	if($line=$self->{'_clone_dbm'}->{$id2}){
	    ($cdate,$type,$cgp,$acc,$sv,$emblid,$htgsp)=split(/,/,$line);
	    if($acc ne $id){
		$self->throw("$id maps to $id2 but does not map back correctly ($acc)");
	    }else{
		$fok=1;
	    }
	}
    }
    if(!$fok){
	$self->throw("$id is not a valid sequence in this database");
    }
    # return $id = name in ensembl (determined by _byacc); $id2 = name on disk
    return $id,$id2,$cgp,$sv,$emblid,$htgsp;
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
    if( $obj->{'_accession_dbm'} ) {
	dbmclose(%{$obj->{'_accession_dbm'}});
	$obj->{'_accession_dbm'} = undef;
    }
}



