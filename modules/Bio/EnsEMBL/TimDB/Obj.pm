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
use NDBM_File;
use Fcntl qw( O_RDONLY );

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ObjI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$raclones,$noacc,$test,$part,$live,@args) = @_;

  # DEFAULT IS NOW LIVE
  $live=1;

  # if we have a list of clones and live not test, write a lock file
  if($live && !$test && $raclones){
      $self->_write_lockfile($raclones);
  }

  # DEBUG
  # second parameter is for debugging to avoid reading entire list of objects
  #if($raclones){
  #    $self->warn("DEBUG: only exon/transcript/gene objects associated with clone list are read");
  #}

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
  # FIXME
      # MACHINE SPECIFIC CONFIG
      if($ENV{'HOST'} eq 'sol28'){
	  $UNFIN_ROOT="/net/nfs0/vol0/home/elia/unfinished_ana";
	  $UNFIN_DATA_ROOT=$UNFIN_ROOT;
      # ADD NEW MACHINES HERE
      }elsif($ENV{'HOST'} eq 'croc'){
      	  $UNFIN_ROOT="/nfs/disk89/michele/pogdir/timdb/";
         $UNFIN_DATA_ROOT=$UNFIN_ROOT;
      }elsif($ENV{'HOST'} eq 'humsrv1'){
      	  $UNFIN_ROOT="/nfs/disk100/humpub/th/unfinished_ana/tmp1";
	  $UNFIN_DATA_ROOT=$UNFIN_ROOT;
      }
      $unfinished_root="$UNFIN_ROOT";
      $unfinished_root.="/test";
      $self->{'_test'}=1;
      $exon_file="$unfinished_root/test_confirmed_exon";
  }else{
      $exon_file="$CONFIRMED_EXON_FASTA";
  }
  $self->{'_unfinished_root'}=$unfinished_root;
  my $clone_dbm_file="$unfinished_root/unfinished_clone.dbm";
  my %unfin_clone;
  unless(tie(%unfin_clone,'NDBM_File',$clone_dbm_file,O_RDONLY,0644)){
      $self->throw("Error opening clone dbm file");
  }
  $self->{'_clone_dbm'}=\%unfin_clone;

  # if going to do things !$noacc then need to open this dbm file too
  my %unfin_accession;
  if(!$noacc){
      my $accession_dbm_file="$unfinished_root/unfinished_accession.dbm";
      unless(tie(%unfin_accession,'NDBM_File',$accession_dbm_file,O_RDONLY,0644)){
	  $self->throw("Error opening accession dbm file");
      }
      $self->{'_accession_dbm'}=\%unfin_accession;
#      print "ACC: $accession_dbm_file: ".scalar(keys %{$self->{'_accession_dbm'}})."\n";
  }

  # clone update file access
  my $clone_update_dbm_file="$unfinished_root/unfinished_clone_update.dbm";
  my %unfin_clone_update;
  unless(tie(%unfin_clone_update,'NDBM_File',$clone_update_dbm_file,O_RDONLY,0644)){
      $self->throw("Error opening clone update dbm file");
  }
  $self->{'_clone_update_dbm'}=\%unfin_clone_update;

  # define a few other important files, depending on options
  my $file_root;
  if($test){
      $self->warn("Using -test: fake test data");
      $file_root="$unfinished_root";
  }elsif($part){
      $self->warn("Using -part: to take g/t/co files from test_gtc/ [development option]");
      $file_root="$unfinished_root/test_gtc";
  }elsif($live){
      # this is now default behaviour - warning turned off
      #$self->warn("Using -live to access live version: may be data inconsistencies");
      $file_root="$unfinished_root";
  }else{
      $self->warn("Using current stable release version of e/t/g/co files");
      $file_root="$unfinished_root/release/current";
      $exon_file="$unfinished_root/release/current/confirmed_exon";
  }

  my $transcript_file="$file_root/unfinished_ana.transcript.lis";
  my $gene_file="$file_root/unfinished_ana.gene.lis";
  my $contig_order_file="$file_root/unfinished_ana.contigorder.lis";
  if(!-e $exon_file){
      $self->throw("Could not access exon file");
  }
  $self->{'_exon_file'}=$exon_file;

  if(!-e $transcript_file){
      $self->throw("Could not access transcript file");
  }
  $self->{'_transcript_file'}=$transcript_file;

  if(!-e $gene_file){
      $self->throw("Could not access gene file");
  }
  $self->{'_gene_file'}=$gene_file;

  if(!-e $contig_order_file){
      $self->throw("Could not access contig order file");
  }
  $self->{'_contig_order_file'}=$contig_order_file;

  # build mappings from these flat files
  # (better to do it here once than each time we need the information!)
  # FIXME - should this be moved to the pipeline so that this information
  # is stored in DBM files - currently in legacy parser
  my $p=Bio::EnsEMBL::Analysis::LegacyParser->new($self->{'_gene_file'},
						  $self->{'_transcript_file'},
						  $self->{'_exon_file'},
						  $self->{'_contig_order_file'});
  $self->{'_parser_object'}=$p;

  # need a full list if $raclones not set,
  # but also need to check clones in list provided to see if they are valid
  # list needs to include ones with invalid SV's as might be called
  # for dumping before having accession numbers
  my @clones=$self->get_all_Clone_id(1);
  if(!$raclones){
      $raclones=\@clones;
  }else{
      my %clones=map {$_,1} @clones;
      my @okclones;
      foreach my $clone (@$raclones){
	  next unless $clone;
	  my $fok;
	  if($clones{$clone}){
	      # see if maps directly
	      push(@okclones,$clone);
	      $fok=1;
	  }else{
	      # see if maps via a translation
	      my($clone2)=$self->get_id_acc($clone,1);
	      next if $clone2 eq 'unk';
	      if($clones{$clone2}){
		  push(@okclones,$clone2);
		  $fok=1;
	      }
	  }
	  if(!$fok){
	      $self->warn("Clone $clone is not recognised or locked");
	  }
      }
      $raclones=\@okclones;
  }

  # keep a record of which clones have been loaded
  %{$self->{'_active_clones'}}=map{$_,1} @$raclones;

  # doing conversion acc->id->acc or id->acc, need it here too
  $p->map_contigorder($self,\@clones);

  #$self->map_etg;

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

    $self->map_etg unless $self->{'_mapped'};
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

sub get_Clone {
    my ($self,$id) = @_;

    # translate incoming id to ensembl_id, taking into account nacc flag
    my($disk_id,$cgp,$sv,$emblid,$htgsp);
    ($id,$disk_id,$cgp,$sv,$emblid,$htgsp)=$self->get_id_acc($id);
    if($id eq 'unk'){
	$self->throw("Cannot get accession for $disk_id");
    }

    # if its already been created, get it from the hash
    if($self->{'_clone_array'}->{$id}){
	return $self->{'_clone_array'}->{$id};
    }

    # else, have to build it

    # can only build it, if it was 'loaded' in the initial call to timdb
    # (i.e. that it is in the active list)
    unless($self->{'_active_clones'}->{$id}){
	$self->throw("$id fetched - not loaded in original object call - locked?");
    }

    # test if clone is not locked (for safety); don't check for valid SV's
    # (probably overkill, as locking at Obj.pm level)
    my($flock,$fsv,$facc)=$self->_check_clone_entry($disk_id);
    if($flock){
	$self->throw("$id is locked by TimDB");
    }
    if($facc){
	$self->throw("$id does not have an accession number");
    }

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
    # save it to hash
    $self->{'_clone_array'}->{$id}=$clone;

    return $clone;
}


=head2 get_all_Clone_id

 Title   : get_all_Clone_id
 Usage   : @cloneid = $obj->get_all_Clone_id($flag)
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    : if $flag set, returns all clones regardless of invalid SV

Note: for speed this does not return the ensembl_id but the disk_id.


=cut

sub get_all_Clone_id{
    my ($self,$fall) = @_;
    return &_get_Clone_id($self,$fall);
}


# Following is not currently coded

=head2 get_updated_Objects

 Title   : get_updated_Objects
 Usage   : @cloneid = $obj->get_updated_Clone_id($last,$now_offset,flag)
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    : if $flag set, returns all clones regardless of invalid SV

=cut

sub get_updated_Objects{
    my($self)=@_;
    $self->throw("call not supported - use get_updated_Clone_id instead");
}


=head2 get_updated_Clone_id

 Title   : get_updated_Clone_id
 Usage   : @cloneid = $obj->get_updated_Clone_id($last,$now_offset,flag)
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    : if $flag set, returns all clones regardless of invalid SV

=cut

sub get_updated_Clone_id{
    my ($self,$last,$now_offset,$fall) = @_;
    my @objs;

    # FIXME - time offset of 30mins = 30x60 seconds
    my $offset_time = 1800;

    $last = $last - $offset_time;
    
    # get list of updated clones
    my @clones;

    my($val,$key);

    while(($key,$val) = each %{$self->{'_clone_update_dbm'}}){

	my($date2)=split(',',$val);
	# make list of updatable clones

	if($date2 > $last && $date2 < $now_offset){
	    push(@clones,$key);
	}
    }

    # check validity of clones selected, then fetch
    return &_get_Clone_id($self,$fall,\@clones);
}


=head2 _get_Clone_id, _check_clone_entry

 Title   : get_Clone_id, _check_clone_entry
 Usage   : private methods
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub _get_Clone_id{
   my ($self,$fall,$ralist) = @_;
   my @list;
   my $nc=0;
   my $nisv=0;
   my $nsid=0;
   my $nlock=0;
   if($ralist){
       # loop over list of clones supplied [unknown]
       foreach my $id (@$ralist){

	   # translate incoming id to ensembl_id, disk_id, taking into account nacc flag
	   my $disk_id;
	   ($id,$disk_id)=$self->get_id_acc($id);

	   # in cases where nacc flag set to return ensembl_id and emsembl_id missing
	   next if $id eq 'unk';

	   my($flock,$fsv,$facc)=$self->_check_clone_entry($disk_id,\$nc,\$nsid,
							   \$nisv,\$nlock);
	   if(!$flock && ($fall || !$fsv) && !$facc){
	       push(@list,$id);
	   }
       }
   }else{
       # loop over whole dbm file [disk_id]
       my($id,$val,$disk_id);
       while(($disk_id,$val)=each %{$self->{'_clone_dbm'}}){

	   my($flock,$fsv,$facc)=$self->_check_clone_entry($disk_id,\$nc,\$nsid,
							   \$nisv,\$nlock);
	   if(!$flock && ($fall || !$fsv) && !$facc){

	       # translate incoming disk_id to ensembl_id, taking into account nacc flag
	       # at this stage know there is an ensembl_id set, so translation is possible
	       ($id)=$self->get_id_acc($disk_id);
	       if($id eq 'unk'){
		   $self->throw("Unexpected failure of get_id_acc");
	       }
	       push(@list,$id);
	   }
       }
   }
   print STDERR "TimDB status:\n";
   if($ralist){
       print STDERR "  $nc clones have been updated\n";
   }else{
       print STDERR "  $nc clones in database\n";
   }
   print STDERR "    $nsid have cloneid rather than accession numbers [Sanger]\n";

   print STDERR "    $nlock clones are locked for reading and are excluded\n";
   print STDERR "    $nisv have invalid SV numbers";
   if($fall){
       print STDERR " and are included\n";
   }else{
       print STDERR " and are excluded\n";
   }
   print STDERR "  Total of ".scalar(@list)." clones are in final list\n";

   # !! no idea why sorting is necessary !!
   return sort @list;
}

# checks a clone entry in TimDB
# posibilities are:
# - clone is locked/unlocked
# - clone has SV/no SV
# - clone has no ACC
# - clone different from accession (information only)
# returns lock and sv state to allow external decision about accepting clone
# can increment counters
sub _check_clone_entry{
    my($self,$disk_id,$rnc,$rnsid,$rnisv,$rnlock)=@_;
    my $val;
    unless($val=$self->{'_clone_dbm'}->{$disk_id}){
	$self->throw("ERROR: $disk_id not in clone DBM");
    }
    $$rnc++;

    my($cdate,$type,$cgp,$acc,$sv,$emblid,$htgsp)=split(/,/,$val);
    # count cases where cloneid is not accession (for information purposes)
    if($disk_id ne $acc){
	$$rnsid++;
    }

    # flag where accession missing - not to go to ensembl
    my $facc;
    if(!$acc && $self->{'_byacc'}){
	$facc=1;
    }

    # count where sv is invalid
    my $fsv=0;
    if($sv!~/^\d+$/){
	$$rnisv++;
	$fsv=1;
    }

    # skip locked clones
    my $val2;
    my $flock=0;
    if($val2=$self->{'_clone_update_dbm'}->{$disk_id}){
	my($date2,$lock)=split(',',$val2);
	if($lock){
	    $$rnlock++;
	    $flock=1;
	}
    }
    return($flock,$fsv,$facc);
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
    my($self,$id,$t)=@_;
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
		# if can't return an ensembl_id when requested to, fail
		# gracefully
		return 'unk',$id;
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

# write lockfile so timdb knows which clones are being processed
# outside timdb.
sub _write_lockfile{
    my($self,$raclones)=@_;
    my $dir="$UNFIN_ROOT/timdb-lock";
    if(-d $dir){
	my $time=time;
	my $file="$dir/timdb-lock.$time.$$";
	local *OUT;
	if(open(OUT,">$file")){
	    print OUT "# host=".$ENV{'HOST'}." user=".$ENV{'USER'}."\n";
	    foreach my $clone (@$raclones){
		print OUT "$clone\n";
	    }
	    close(OUT);
	    push(@{$self->{'_lockfile'}},$file);
	}else{
	    $self->warn("Could not write lock file $file");
	}
    }else{
	$self->warn("Could not read lock directory $dir");
    }
}

# build the exon/transcript/gene map
sub map_etg{
    my($self)=shift;
    my $p=$self->{'_parser_object'};
    my @clones=(keys %{$self->{'_active_clones'}});
    $self->warn("DEBUG: initialising map of TimDB for ".scalar(@clones)." clones");
    # doing conversion acc->id->acc or id->acc, need it here too
    $p->map_all($self,\@clones);
    $self->warn("DEBUG: TimDB initialisation complete");
    $self->{'_mapped'}=1;
}

# close the dbm clone file, remove lock
sub DESTROY{
    my ($obj) = @_;
    if( $obj->{'_clone_dbm'} ) {
	untie %{$obj->{'_clone_dbm'}};
	$obj->{'_clone_dbm'} = undef;
    }
    if( $obj->{'_accession_dbm'} ) {
	untie %{$obj->{'_accession_dbm'}};
	$obj->{'_accession_dbm'} = undef;
    }
    if( $obj->{'_clone_update_dbm'} ) {
	untie %{$obj->{'_clone_update_dbm'}};
	$obj->{'_clone_update_dbm'} = undef;
    }
    # remove lock file
    if($obj->{'_lockfile'}){
	foreach my $file (@{$obj->{'_lockfile'}}){
	    unlink $file;
	}
    }
}

#END {
#    local *DIR;
#    my $dir="$UNFIN_ROOT/timdb-lock";
#    if(opendir(DIR,$dir)){
#	my @files=readdir(DIR);
#	closedir(DIR);
#	foreach my $file (@files){
#	    if($file=~/\.$$/){
#		unlink "$dir/$file";
#	    }
#	}
#    }
#}
