package XrefMapper::ProcessPaired;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->xref($mapper->xref);
  $self->verbose($mapper->verbose);
  return $self;
}


sub process{
  my ($self) = @_;



  print "Process Pairs\n" if($self->verbose);
  my $object_xref_id;

  my $sth =  $self->xref->dbc->prepare("select MAX(object_xref_id) from object_xref");
  $sth->execute;
  $sth->bind_columns(\$object_xref_id);
  $sth->fetch;
  $object_xref_id++;
  $sth->finish;

  print "Starting at object_xref of $object_xref_id\n" if($self->verbose);

  my $psth = $self->xref->dbc->prepare("select p.accession1, p.accession2 from pairs p");
  my $ox_count_sth =  $self->xref->dbc->prepare('select count(1) from object_xref ox, xref x where ox.xref_id = x.xref_id and ox.ox_status = "DUMP_OUT" and x.accession = ?');
 
  my $ox_transcript_sth =   $self->xref->dbc->prepare('select gtt.transcript_id, ix.query_identity, ix.target_identity  from identity_xref ix, object_xref ox, xref x, gene_transcript_translation gtt where ox.object_xref_id = ix.object_xref_id and ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and gtt.translation_id = ox.ensembl_id and x.accession = ?');

  my $ox_translation_sth =  $self->xref->dbc->prepare('select gtt.translation_id, ix.query_identity, ix.target_identity from identity_xref ix, object_xref ox, xref x, gene_transcript_translation gtt where ox.object_xref_id = ix.object_xref_id and ox.ox_status = "DUMP_OUT" and ox.xref_id = x.xref_id and gtt.transcript_id  = ox.ensembl_id and x.accession = ?');
 
  my $xref_sth =  $self->xref->dbc->prepare("select xref_id from xref where accession = ?");
 
  my $ox_insert_sth = $self->xref->dbc->prepare("insert into object_xref (object_xref_id, xref_id, ensembl_id, ensembl_object_type, linkage_type, ox_status) values(?, ?, ?, ?, 'INFERRED_PAIR', 'DUMP_OUT')");
  local $ox_insert_sth->{RaiseError}; #catch duplicates
  local $ox_insert_sth->{PrintError}; # cut down on error messages



  my $ox_get_id_sth = $self->xref->dbc->prepare("select object_xref_id,ox_status from object_xref where xref_id = ? and ensembl_id = ? and ensembl_object_type = ?");

  my $ox_update_sth =  $self->xref->dbc->prepare('update object_xref set ox_status = "DUMP_OUT", linkage_type = "INFERRED_PAIR" where object_xref_id = ?');
  my $xref_update_sth =  $self->xref->dbc->prepare('update xref set info_type = "INFERRED_PAIR" where xref_id = ?');
  my $ins_dep_ix_sth = $self->xref->dbc->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, ?, ?)");

  $psth->execute() || die "execute failed";
  my ($acc1, $acc2);
  $psth->execute();
  
  my $refseq_count = 0;

  my %change;

  $psth->bind_columns(\$acc1, \$acc2);
  while($psth->fetch()){
    my $count1;
    my $count2;
    $ox_count_sth->execute($acc1);  # translation alignment
    $ox_count_sth->bind_columns(\$count1);
    $ox_count_sth->fetch;

    $ox_count_sth->execute($acc2);  # transcript alignment
    $ox_count_sth->bind_columns(\$count2);
    $ox_count_sth->fetch;

    if(( $count1 and $count2) || (!($count1) and !($count2)) ){
      next; # eithr both matched or neither is.
    }	
    if($count1){ 
      #need xref_id for acc2
      my $xref_id;
      $xref_sth->execute($acc2);
      $xref_sth->bind_columns(\$xref_id);
      if(!$xref_sth->fetch){
#	 print "Could not find xref_id for accession $acc2\n";
	 next;
       }
      next if(!defined($xref_id));
#      print "$acc2\t$xref_id (search using $acc1)\n";
     # insert new object_xref
      # trap error code. if duplicate then just set linkage_type = "INFERRED_PAIR" and ox_status = "DUMP"
      # "maybe" the original failed the cutoff!!! so will have an entry alread but no good.
      $ox_transcript_sth->execute($acc1);
      my $transcript_id=undef;
      my ($q_id,$t_id);
      
      $ox_transcript_sth->bind_columns(\$transcript_id,\$q_id, \$t_id);
      while($ox_transcript_sth->fetch){
	if(defined($transcript_id)){ # remember not all transcripts have translations.

	  $object_xref_id++;
	  $ox_insert_sth->execute($object_xref_id, $xref_id, $transcript_id, "Transcript") ;
	  if($ox_insert_sth->err){
	    my $err = $ox_insert_sth->errstr;
	    if($err =~ /Duplicate/){
	      $change{"UPDATE"}++;
	      # duplicate this can happen as it might have failed the cutoff
	      # find the old object_xref_id and the update the status's
	      my $old_object_xref_id=undef;
	      my $status;
	      $ox_get_id_sth->execute($xref_id, $transcript_id, "Transcript");	      
	      $ox_get_id_sth->bind_columns(\$old_object_xref_id, \$status);
	      $ox_get_id_sth->fetch();
	      if($status eq "DUMP_OUT"){
		print STDERR "Problem status for object_xref_id is DUMP_OUT but this should never happen as it was not found earlier??? (transcript_id = $transcript_id, $xref_id\n";
	      }
	      if(!defined($old_object_xref_id)){
		die "Duplicate but can't find the original?? xref_id = $xref_id, ensembl_id = $transcript_id, type = Transcript\n";
	      }
	      $ox_update_sth->execute($old_object_xref_id)|| die "Could not set update for object_xref_id = $old_object_xref_id";
	      $xref_update_sth->execute($xref_id)|| die "Could not set update for xref_id = $xref_id";
	    }
	    else{
	      die "Problem loading error is $err\n";
	    } 
	  }
	  else{
	    $ins_dep_ix_sth->execute($object_xref_id, $q_id, $t_id);
	    $xref_update_sth->execute($xref_id)|| die "Could not set update for xref_id = $xref_id";
	    $change{"NEW"}++;
	  }
#  	  print "insert $xref_id transcript $transcript_id ........\n";
	  $refseq_count++;
	}

      }
    }
    elsif($count2){
      my $xref_id;
      $xref_sth->execute($acc1);
      $xref_sth->bind_columns(\$xref_id);
      if(!$xref_sth->fetch){
#	print "Could not find xref_id for accession $acc1\n";
	next;
      }
      next if(!defined($xref_id));
#      print "$acc1\t$xref_id (search using $acc2)\n";
      # insert new object_xref
      # trap error code. if duplicate then just set linkage_type = "INFERRED_PAIR" and ox_status = "DUMP"
      # "maybe" the original failed the cutoff!!! so will have an entry alread but no good.
      $ox_translation_sth->execute($acc2);
      my $translation_id = undef;
      my ($q_id, $t_id);
      $ox_translation_sth->bind_columns(\$translation_id, \$q_id, \$t_id);
      while($ox_translation_sth->fetch){
	if(defined($translation_id)){ # remember not all transcripts ahve translations.
	  $object_xref_id++;
	  $ox_insert_sth->execute($object_xref_id, $xref_id, $translation_id, "Translation") ;
	  if($ox_insert_sth->err){
	    $change{"UPDATE"}++;
	    my $err = $ox_insert_sth->errstr;
	    if($err =~ /Duplicate/){
	      # duplicate this can happen as it might have failed the cutoff
	      # find the old object_xref_id and the update the status's
	      my $old_object_xref_id=undef;
	      my $status;
	      $ox_get_id_sth->execute($xref_id, $translation_id, "Translation");	      
	      $ox_get_id_sth->bind_columns(\$old_object_xref_id,\$status);
	      $ox_get_id_sth->fetch();
	      if($status eq "DUMP_OUT"){
		print STDERR "Problem status for object_xref_id is DUMP_OUT but this should never happen as it was not found earlier??? (trasnlation_id = $translation_id, $xref_id\n";
	      }
	      if(!defined($old_object_xref_id)){
		die "Duplicate but can't find the original?? xref_id = $xref_id, ensembl_id = $translation_id, type = Translation\n";
	      }
	      $ox_update_sth->execute($old_object_xref_id)|| die "Could not set update for object_xref_id = $old_object_xref_id";
	      $xref_update_sth->execute($xref_id)|| die "Could not set update for xref_id = $xref_id";
	    }
	    else{
	      die "Problem loading error is $err\n";
	    } 
	  }
	  else{
	    $ins_dep_ix_sth->execute($object_xref_id, $q_id, $t_id);
	    $xref_update_sth->execute($xref_id)|| die "Could not set update for xref_id = $xref_id";
	    $change{"NEW"}++;
	  }
#	  print "insert $xref_id translation $translation_id ........\n";
	  $refseq_count++;
	}
      }
    }
    else{
      print STDERR "HMMM how did i get here. This should be impossible. [logic error]\n";
   }
  }
  $psth->finish;
  $ox_count_sth->finish;
  $ox_transcript_sth->finish;
  $ox_translation_sth->finish;
  $ox_update_sth->finish;
  $xref_update_sth->finish;
  $xref_sth->finish;
  $ins_dep_ix_sth->finish;
  foreach my $key (keys %change){
    print "\t$key\t".$change{$key}."\n" if($self->verbose);
  }
  print "$refseq_count new relationships added\n" if($self->verbose);
  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('processed_pairs',now())");
  $sth_stat->execute();
  $sth_stat->finish;
}

1;
