
#
# BioPerl module for DB::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::Obj - Object representing an instance of an EnsEMBL DB

=head1 SYNOPSIS

    $db = new Bio::EnsEMBL::DB::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );

    $clone = $db->get_clone('X45667');

    $contig = $db->get_Contig("dJ52N12.02793");

    $gene  = $db->get_Gene('HG45501');

    

=head1 DESCRIPTION

This object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). From the object you can pull
out other objects by their stable identifier, such as Clone (accession number),
Exons, Genes and Transcripts. The clone gives you a DB::Clone object, from
which you can pull out associated genes and features. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::Obj;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DBSQL::Contig;
use Bio::EnsEMBL::DBSQL::Clone;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use DBI;

use Bio::EnsEMBL::DBSQL::DummyStatement;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  print "Got",join(',',@args),"\n";
  my ($db,$host,$driver,$user,$password,$debug) = 
      $self->_rearrange([qw(DBNAME
			    HOST
			    DRIVER
			    USER
			    PASS
			    DEBUG
			    )],@args);
  print "Got $db is db and $user as user\n";

  $db || $self->throw("Database object must have a database name");
  $user || $self->throw("Database object must have a user");

  #
  # This needs to be rethought. We are caching sequences
  # here to allow multiple exons to be retrieved fine
  #
  $self->{'_contig_seq_cache'} = {};
  

  $self->{'_lock_table_hash'} = {};

  if( $debug ) {
     $self->_debug($debug);
 } else {
     $self->_debug(0);
 }
  
  if( ! $driver ) {
      $driver = 'mysql';
  }
  if( ! $host ) {
      $host = 'localhost';
  }
  my $dsn = "DBI:$driver:database=$db;host=$host";

  if( $debug && $debug > 10 ) {
      $self->_db_handle("dummy dbh handle in debug mode $debug");
  } else {

      my $dbh = DBI->connect("$dsn","$user",$password);

      $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");
      
      if( $self->_debug > 3 ) {
	  $self->warn("Using connection $dbh");
      }
     
      $self->_db_handle($dbh);
  }

# set stuff in self from @args
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

   $geneid || $self->throw("Attempting to create gene with no id");

   my $gene = Bio::EnsEMBL::Gene->new();

   # check the database is sensible before we yank out this gene. Yes this is
   # paranoid.

   my $sth1 = $self->prepare("select p1.contig from exon as p1, transcript as p2, exon_transcript as p3 where p2.gene = '$geneid' and p2.id = p3.transcript and p3.exon = p1.id");
   $sth1->execute();
   while( my $rowhash = $sth1->fetchrow_hashref) {
       # get a contig object, which checks that this exists. 
       my $contig = $self->get_Contig($rowhash->{'contig'});
       # if there is no exception then it is there. Get rid of it
       $contig = 0;
   }


   # go over each Transcript
   my $sth = $self->prepare("select id,translation from transcript where gene = '$geneid'");

   my $res = $sth->execute();
   my $seen =0;
   while( my $rowhash = $sth->fetchrow_hashref) {
       my $trans = $self->get_Transcript($rowhash->{'id'});
       my $translation = $self->get_Translation($rowhash->{'translation'});
       $trans->translation($translation);
       $gene->add_Transcript($trans);
       $seen = 1;
   }
   
   if( $seen == 0 ) {
       $self->throw("No gene with $geneid as a name! - Sorry!");
   }
   $gene->id($geneid);

   return $gene;
}

=head2 get_last_update

 Title   : get_last_update
 Usage   : $obj->get_last_update; 
 Function: Reads the meta table of the database to get the last_update time
 Example : get_last_update
 Returns : UNIX TIME of last update
 Args    : none


=cut

sub get_last_update{
    my ($self) = @_;
    
    my $sth = $self->prepare("select last_update from meta");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $datetime = $rowhash->{'last_update'};
    $sth = $self->prepare("select UNIX_TIMESTAMP('".$datetime."')");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    return $rowhash->[0];
}

=head2 get_offset_time

 Title   : get_offset_time
 Usage   : $obj->get_offset_time; 
 Function: Reads the meta table of the database to get the offset_time
 Example : get_offset_time
 Returns : UNIX TIME of offset_time
 Args    : none


=cut

sub get_offset_time {
    my ($self) = @_;

    my $sth = $self->prepare("select offset_time from meta");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    return $rowhash->{'offset_time'};
}

=head2 get_now_offset

 Title   : get_now_offset
 Usage   : $obj->get_now_minus_offset; 
 Function: Gets the current time from the point of view of the database, substracts the
           offset time found in the meta table and gives back unix time of now-offset
 Example : get_now_offset
 Returns : UNIX TIME of now - offset_time
 Args    : none


=cut

sub get_now_offset{
    my ($self) = @_;

    #First, get now, i.e. current time from db, in datetime format
    my $sth = $self->prepare("select now()");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $now = $rowhash->{'now()'};
    
    #Now get the offset time from the meta table, which is in time format
    $sth = $self->prepare("select offset_time from meta");
    $sth->execute();
    $rowhash = $sth->fetchrow_hashref();
    my $offset = $rowhash->{'offset_time'};

    #Perform the subtraction in mysql
    $sth = $self->prepare("select DATE_SUB(\"$now\", INTERVAL \"$offset\" HOUR_SECOND)");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    my $datetime = $rowhash->[0];

    #Trasform the result into unix time and return it
    $sth = $self->prepare("select UNIX_TIMESTAMP('".$datetime."')");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    return $rowhash->[0];
}
    
    

=head2 get_Protein_annseq

 Title   : get_Protein_annseq
 Usage   : get_Protein_annseq ($ENSP); 
 Function: Creates an annseq object for a particular peptide, storing the peptide
           sequence in $annseq->seq, and adding all the protein features as generic
           Seqfeatures
 Example : 
 Returns : $annseq
 Args    : $ENSP


=cut

sub get_Protein_annseq{
    my ($self,$ENSP) = @_;
    my $annseq = Bio::EnsEMBL::AnnSeq->new();
    
    my $sth = $self->prepare("select id from transcript where translation = '$ENSP'");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript = $self->get_Transcript($rowhash->{'id'});
    my $translation = $self->get_Translation($ENSP);
    $transcript->translation($translation);
    my $seq = $transcript->translate();
    $annseq->seq($seq);
    $sth = $self->prepare("select * from proteinfeature where translation = '$ENSP'");
    $res = $sth->execute();
    while( my $rowhash = $sth->fetchrow_hashref) {
	my $analysis = $rowhash->{'analysis'};
	my $sth2 = $self->prepare("select * from analysis where id = '$analysis'");
	my $res2 = $sth2->execute();
	my $rowhash2 = $sth2->fetchrow_hashref;
	my $feature = new Bio::SeqFeature::Generic ( -start => $rowhash->{'seq_start'}, 
						     -end => $rowhash->{'seq_end'},
						     -score =>  $rowhash->{'score'},
						     -primary => $rowhash2->{'gff_feature'},
						     -source => $rowhash2->{'gff_source'});
	$annseq->add_SeqFeature($feature);
    }
    
    return $annseq;   
}

=head2 get_Transcript
    
 Title   : get_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub get_Transcript{
    my ($self,$transid) = @_;
    
    my $trans = Bio::EnsEMBL::Transcript->new();
    # go over each Transcript
    my $sth = $self->prepare("select exon from exon_transcript where transcript = '$transid'");

   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {
       my $exon = $self->get_Exon($rowhash->{'exon'});
       $trans->add_Exon($exon);
   }
   $trans->id($transid);

   return $trans;
}

=head2 get_Translation

 Title   : get_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Translation{
   my ($self,$translation_id) = @_;

   my $sth = $self->prepare("select seq_start,start_exon,seq_end,end_exon from translation where id = '$translation_id'");
   my $res = $sth->execute();
   my $rowhash = $sth->fetchrow_hashref;
   my $out = Bio::EnsEMBL::Translation->new();
   $out->start($rowhash->{'seq_start'});
   $out->end($rowhash->{'seq_end'});
   $out->start_exon_id($rowhash->{'start_exon'});
   $out->end_exon_id($rowhash->{'end_exon'});
   $out->id($translation_id);

   return $out;
}

=head2 get_Exon

 Title   : get_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Exon{
   my ($self,$exonid) = @_;

   my $sth = $self->prepare("select id,contig,created,modified,seq_start,seq_end,strand,phase from exon where id = '$exonid'");
   $sth->execute;
   my $rowhash = $sth->fetchrow_hashref;

   my $exon = Bio::EnsEMBL::Exon->new();
   $exon->contig_id($rowhash->{'contig'});

   my $contig_id = $exon->contig_id();

   # we have to make another trip to the database to get out the contig to clone mapping.
   my $sth2 = $self->prepare("select clone from contig where id = '$contig_id'");
   $sth2->execute;
   my $rowhash2 = $sth2->fetchrow_hashref;

   $exon->clone_id($rowhash2->{'clone'});

   # rest of the attributes
   $exon->id($rowhash->{'id'});
   $exon->created($rowhash->{'created'});
   $exon->modified($rowhash->{'modified'});
   $exon->start($rowhash->{'seq_start'});
   $exon->end($rowhash->{'seq_end'});
   $exon->strand($rowhash->{'strand'});
   $exon->phase($rowhash->{'phase'});
   
   # we need to attach this to a sequence. For the moment, do it the stupid
   # way perhaps?

   my $seq;

   if( $self->_contig_seq_cache($exon->contig_id) ) {
       $seq = $self->_contig_seq_cache($exon->contig_id);
   } else {
       my $contig = $self->get_Contig($exon->contig_id());
       $seq = $contig->seq();
       $self->_contig_seq_cache($exon->contig_id,$seq);
   }

   $exon->attach_seq($seq);

   return $exon;
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

   my  $sth = $self->prepare("select id from contig where clone = \"$id\";");
   $sth->execute();
   my  $rv = $sth->rows;
   if( ! $rv ) {
       # make sure we deallocate sth - keeps DBI happy!
       $sth = 0;
       $self->throw("Clone $id does not seem to occur in the database!");
   }

   my $clone = new Bio::EnsEMBL::DBSQL::Clone( -id => $id,
					-dbobj => $self );

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
   my ($self,$id) = @_;

   my $sth = $self->prepare("select p1.id,p2.id from dna as p1,contig as p2 where p2.id = '$id'");

   $sth->execute;
   my $rowa = $sth->fetchrow_arrayref;
   
   if( ! $rowa->[0] || ! $rowa->[1] ) {
       $self->throw("Contig $id does not exist in the database or does not have DNA sequence");
   }

   my $contig = new Bio::EnsEMBL::DBSQL::Contig ( -dbobj => $self,
					       -id => $id );
   return $contig;
}
=head2 get_all_Clone_id

 Title   : get_all_Clone_id
 Usage   : @cloneid = $obj->get_all_Clone_id
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Clone_id{
   my ($self) = @_;
   my $sth = $self->prepare("select id from clone");
   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
}

=head2 get_all_Gene_id

 Title   : get_all_Gene_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Gene_id{
   my ($self) = @_;
   my $sth = $self->prepare("select id from gene");
   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
}

=head2 get_updated_objects
    
 Title   : get_updated_objects
 Usage   : $obj->get_updated_objects ($recipient_last_update, $recipient_now, $recipient_offset)
 Function: Gets all the objects that have been updated (i.e.change in 
	   version number) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_updated_objects (973036800,973090800)
 Returns : all the objects updated within that timespan
 Args    : $recipient_last_update, $recipient_now

=cut

sub get_updated_objects{
    my ($self, $last, $now_offset) = @_;
    
    $last || $self->throw("Attempting to get updated objects without the recipient db last update time");
    $now_offset  || $self->throw("Attempting to get updated objects without the recipient db current time");

    #First, let us convert the unix times now_offset and last into mysql times
    my $sth = $self->prepare("select FROM_UNIXTIME(".$last.")");
    $sth->execute();
    my $rowhash = $sth->fetchrow_arrayref();
    $last = $rowhash->[0];

    $sth = $self->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];
    

    #First, get all clone ids that have been updated between last and now-offset
    my $sth = $self->prepare("select id from clone where modified > '".$last."' and modified <= '".$now_offset."'");
   
    $sth->execute;
    my @out;
    my @clones;
    while( my $rowhash = $sth->fetchrow_hashref) {
	push(@clones,$rowhash->{'id'});
    }
    
    #Get all clone objects for the ids contained in @clones, and push them in @out
    foreach my $cloneid (@clones) {
	push @out, $self->get_Clone ($cloneid);
    }	
    
    #Get all gene ids that have been updated between last and now-offset
    $sth = $self->prepare("select id from gene where modified > '".$last."' and modified <= '".$now_offset."'");
    $sth->execute;

    my @genes;
    while( $rowhash = $sth->fetchrow_hashref) {
	push(@genes,$rowhash->{'id'});
    }
    
    #Get all clone objects for the ids contained in @clones, and push them in @out
    foreach my $geneid (@genes) {
	push @out, $self->get_Gene ($geneid);
    }	
    return @out;
}

=head2 archive_Gene
    
 Title   : archive_Gene
 Usage   : $obj->archive_gene($gene,$clone,$arcdb)
 Function: Deletes a gene and all its transcripts and exons, 
           and archives partial info in the archive db passed on.
 Example : 
 Returns : nothing
 Args    : $gene, $clone, $arcdb (archive database object)


=cut

sub archive_Gene {
   my ($self,$gene,$clone,$arc_db) = @_;
   my $sth;

   # get transcripts for the gene given 

   foreach my $transcript ($gene->each_Transcript) {
       
       #Get out transcript info needed to write into archive db
       my $seq = $transcript->dna_seq;
       $seq->id($transcript->id);
       print STDERR "The transcript sequence id is ".$seq->id."\n";

       #Temporary, since versions not stored yet...
       !$transcript->version && $transcript->version(1);
       !$gene->version && $gene->version(1);
       !$clone->version && $clone->version(1);
       
       #Finally, write all the info to a new entry in the archive database

       $arc_db->write_seq($seq, $transcript->version, 'transcript', $gene->id, $gene->version, $clone->id, $clone->version);
       
       #Get out translation to write protein into archive db
       #Note: version is the one from transcript!

       $seq = $transcript->translate;
       print STDERR "The protein sequence id is ".$seq->id."\n";

       $arc_db->write_seq($seq, $transcript->version, 'protein', $gene->id, $gene->version, $clone->id, $clone->version);

       #Delete transcript rows
       $sth= $self->prepare("delete from transcript where id = '".$transcript->id."'");
       $sth->execute;
       
       #Get out exons for this transcript
       foreach my $exon ($transcript->each_Exon) {
	   
	   #Get out info needed to write into archive db
	   $seq = $exon->seq;
	   $seq->id($exon->id);
	   print STDERR "The exon sequence id is ".$exon->id."\n";
	   #Temporary, since versions not stored yet...
	   !$exon->version && $exon->version(1);
	   
	   #Write into archive db
	   $arc_db->write_seq($seq, $exon->version, 'exon', $gene->id, $gene->version, $clone->id, $clone->version);
	   
           #Delete exon_transcript rows
	   $sth= $self->prepare("delete from exon_transcript where transcript = '".$transcript->id."'");
	   $sth->execute;
	   #Delete exon rows
	   $sth = $self->prepare("delete from exon where id = '".$exon->id."'");
	   $sth->execute;
       }
   }
   
   # delete gene rows
   $sth = $self->prepare("delete from gene where id = '".$gene->id."'");
   $sth->execute;
}   

=head2 delete_Clone

 Title   : delete_Clone
 Usage   : $obj->delete_Clone($clone_id)
 Function: Deletes clone, including contigs, but not its genes
 Example :
 Returns : 
 Args    :


=cut

sub delete_Clone{
   my ($self,$clone_id,$arc_db) = @_;
   
   $clone_id || $self->throw ("Trying to delete clone without a clone_id\n");
   $arc_db || $self->throw ("Not allowed to delete clones without passing on a valid, connected, archive_DB!\n");
   
   my @contigs;
   # get a list of contigs to zap
   my $sth = $self->prepare("select id from contig where clone = '$clone_id'");
   
   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@contigs,$rowhash->{'id'});
   }
   
   
   # Delete from DNA table, Contig table, Clone table
   
   foreach my $contig ( @contigs ) {
       my $sth = $self->prepare("delete from contig where id = '$contig'");
       $sth->execute;
       $sth = $self->prepare("delete from dna where contig = '$contig'");
       $sth->execute;
   }
   
   $sth = $self->prepare("delete from clone where id = '$clone_id'");
   $sth->execute;
}

=head2 delete_Gene

 Title   : delete_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub delete_Gene{
   my ($self,$geneid) = @_;
   my @trans;
   my %exon;

   # get out exons, transcripts for gene. 

   my $sth = $self->prepare("select id from transcript where gene = '$geneid'");
   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@trans,$rowhash->{'id'});
   }

   foreach my $trans ( @trans ) {
       my $sth = $self->prepare("select exon from exon_transcript where transcript = '$trans'");
       $sth->execute;
       while( my $rowhash = $sth->fetchrow_hashref) {
	   $exon{$rowhash->{'id'}} =1;
       }
   }

   # delete exons, transcripts, gene rows

   foreach my $exon ( keys %exon ) {
       my $sth = $self->prepare("delete from exon where id = '$exon'");
       $sth->execute;
   }

   foreach my $trans ( @trans ) {
       my $sth= $self->prepare("delete from transcript where id = '$trans'");
       $sth->execute;
       $sth= $self->prepare("delete from exon_transcript where transcript = '$trans'");
       $sth->execute;
   }

   $sth = $self->prepare("delete from gene where id = '$geneid'");
   $sth->execute;
}   
       

=head2 geneid_to_cloneid

 Title   : geneid_to_cloneid
 Usage   : @cloneid = $db->geneid_to_cloneid($geneid);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub geneid_to_cloneid{
   my ($self,$geneid) = @_;

   my $sth = $self->prepare("select p1.id from contig as p1, transcript as p2, exon_transcript as p3 where p2.gene = '$geneid' and p2.id = p3.transcript and p3 ");
   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

}

=head2 cloneid_to_geneid

 Title   : cloneid_to_geneid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub cloneid_to_geneid{
   my ($self,$cloneid) = @_;

   my $sth = $self->prepare("select p2.gene from contig as p1, transcript as p2, exon_transcript as p3 where p2.rtranscript =  and p2.id = p3.transcript and p3 ");
   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
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

   my %done;

   if( ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   # get out unique contig ids from gene to check against
   # database.

   foreach my $contig_id ( $gene->unique_contig_ids() ) {
       eval {
	   my $contig = $self->get_Contig($contig_id);
	   # if there is no exception then it is there. Get rid of it
	   $contig = 0;
       };
       if( $@ ) {
	   $self->throw("In trying to write gene " . $gene->id(). " into the database, unable to find contig $contig_id. Aborting write\n\nFull Exception\n\n$@\n");
	   # done before locks, so we are ok.
       }
       
   }



  # $self->_lock_tables('gene','exon','transcript','exon_transcript');

   # gene is big daddy object


   foreach my $trans ( $gene->each_Transcript() ) {
       $self->write_Transcript($trans,$gene);
       my $c = 1;
       foreach my $exon ( $trans->each_Exon() ) {
	   my $sth = $self->prepare("insert into exon_transcript (exon,transcript,rank) values ('". $exon->id()."','".$trans->id()."',".$c.")");
	   $sth->execute();
	   if( $done{$exon->id()} ) { next; }
	   $done{$exon->id()} = 1;
	   $self->write_Exon($exon);
	   $c++;
       }
   }

   my $sth2 = $self->prepare("insert into gene (id) values ('". $gene->id(). "')");
   $sth2->execute();

   #$self->_unlock_tables();

}


=head2 write_all_Protein_features

 Title   : write_all_Protein_features
 Usage   : $obj->write_all_Protein_features($ENSP)
 Function: writes all protein features of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_all_Protein_features {
    my ($self,$prot_annseq,$ENSP) = @_;
    
    my $c=0;
    foreach my $feature ($prot_annseq->all_SeqFeatures()) {
	my $sth = $self->prepare("insert into proteinfeature (id,seq_start, seq_end, score, analysis, translation) values (NULL,"
				 .$feature->start().","
				 .$feature->end().","
				 .$feature->score().",'"
				 .$c."','"
				 .$ENSP."')");
	$sth->execute();
	
	my $sth2 = $self->prepare("insert into analysis (id,db,db_version,program,program_version,gff_source,gff_feature) values ('$c','testens',1,'elia_program',1,'"
				  .$feature->source_tag()."','"
				  .$feature->primary_tag()."')");
	 $sth2->execute();
	$c++;
    }
}

=head2 write_Protein_feature

 Title   : write_Protein_feature
 Usage   : $obj->write_Protein_feature($ENSP, $feature)
 Function: writes a protein feature object of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_Protein_feature {
    my ($self,$ENSP,$feature) = @_;
    
    my $sth = $self->prepare("insert into proteinfeature (seq_start, seq_end, score, translation) values ("
			     .$feature->start()." ,"
			     .$feature->end()." ,'"
			     .$feature->score()." ,'"
			     .$ENSP."'
				)");
    $sth->execute();
}

=head2 write_Feature

 Title   : write_Feature
 Usage   : $obj->write_Feature($feature)
 Function: Writes a feature on the genomic sequence of a contig into the database
 Example :
 Returns : nothing
 Args    : Bio::SeqFeature::Generic


=cut

sub write_Feature {
    my ($self,$feature,$analysisid,$contig) = @_;

    $self->throw("Feature is not a Bio::SeqFeature::Generic")        unless $feature->isa("Bio::SeqFeature::Generic");
    $self->throw("$contig is not a Bio::EnsEMBL::DBSQL::Contig")     unless $contig ->isa("Bio::EnsEMBL::DBSQL::Contig");

    my $contigid = $contig->id;

    
    my $sth = $self->prepare("insert into feature(id,contig,start,end,score,strand,name,analysis) values (NULL,\"" . 
			     $contig ->id          . "\"," .
			     $feature->start       . "," . 
			     $feature->end         . "," . 
			     $feature->score       . "," . 
			     $feature->strand      . ",\"" . 
			     $feature->source_tag . "\"," .
			     $analysisid           . ")");
    my $rv = $sth->execute();

    return $rv;
}

=head2 write_Analysis

 Title   : write_Analysis
 Usage   : $obj->write_Analysis($anal)
 Function: Writes analysis details to the database
           Checks first whether this analysis entry already exists
 Example :
 Returns : int
 Args    : Bio::EnsEMBL::Analysis::Analysis


=cut

sub write_Analysis {
    my ($self,$anal) = @_;

    $self->throw("Argument is not a Bio::EnsEMBL::Analysis::Analysis") unless $anal->isa("Bio::EnsEMBL::Analysis::Analysis");


    # First check whether this entry already exists.
    my $query;
    if ($anal->db ne "" && $anal->db_version ne "") {
        $query = "select id from analysis where db = \""      . $anal->db              . "\" and" .
                " db_version = \""      . $anal->db_version      . "\" and " .
                " program =    \""      . $anal->program         . "\" and " .
                " program_version = \"" . $anal->program_version . "\" and " .
                " gff_source = \""      . $anal->gff_source      . "\" and" .
                " gff_feature = \""     . $anal->gff_feature     . "\"";
    } else {
        $query = "select id from analysis where " .
                " program =    \""      . $anal->program         . "\" and " .
                " program_version = \"" . $anal->program_version . "\" and " .
                " gff_source = \""      . $anal->gff_source      . "\" and" .
                " gff_feature = \""     . $anal->gff_feature     . "\"";
    }
    my $sth = $self->prepare($query);

    my $rv = $sth->execute();

    # Ony write if we have no result

    # temporary id here

    if ($sth->rows == 0) {
        if ($anal->db ne "" && $anal->db_version ne "") {
            $query = "insert into analysis(id,db,db_version,program,program_version,gff_source,gff_feature) values (NULL,\"" .
                $anal->db               . "\","   .
                $anal->db_version       . ",\""   .
                $anal->program          . "\",\"" .
                $anal->program_version  . "\",\"" .
                $anal->gff_source       . "\",\"" .
                $anal->gff_feature      . "\")";
        } else {
            $query = "insert into analysis(id,program,program_version,gff_source,gff_feature) values (NULL,\"" .
                $anal->program          . "\",\"" .
                $anal->program_version  . "\",\"" .
                $anal->gff_source       . "\",\"" .
                $anal->gff_feature      . "\")";
        }
        print("Query is $query\n");

        my $sth2 = $self->prepare($query);
        my $rv   = $sth2->execute;


        $sth = $self->prepare("select last_insert_id()");
        $rv  = $sth->execute;

        $sth = $self->prepare("select last_insert_id()");
        $rv  = $sth->execute;

        if ($sth->rows == 1) {
            my $rowhash = $sth->fetchrow_hashref;
            return $rowhash->{'last_insert_id()'};
        } else {
            $self->throw("Wrong number of rows returned : " . $sth->rows . " : should be 1");
        }
    } else {
        my $rowhash = $sth->fetchrow_hashref;
        return $rowhash->{'id'};
    }

}


=head2 write_Homol_Feature

 Title   : write_Homol_Feature
 Usage   : $obj->write_Homol_Feature($feature)
 Function: Writes a homol feature on the genomic sequence of a contig into the database
 Example :
 Returns : nothing
 Args    : Bio::SeqFeature::Homol


=cut

sub write_Homol_Feature {
    my ($self,$feature,$analysisid,$contig) = @_;

    $self->throw("Wrong number of arguments to write_Homol_Feature") unless $contig;
    $self->throw("Feature is not a Bio::SeqFeature::Homol")          unless $feature->isa("Bio::SeqFeature::Homol");
    $self->throw("$contig is not a Bio::EnsEMBL::DBSQL::Contig")     unless $contig ->isa("Bio::EnsEMBL::DBSQL::Contig");
    $self->throw("No analysis id input")                             unless $analysisid;

    my $contigid = $contig->id;
    my $homol    = $feature->homol_SeqFeature;

    my $rv = $self->write_Feature($feature,$analysisid,$contig);

    $self->throw("Writing homol feature to the database failed for contig " . $contigid . "\n") unless $rv;

    # Now write into the homol table - reading the manual last_insert_id should
    # work on a per client basis so I shouldn't have to lock the tables.

    my $sth    = $self->prepare("insert into homol_feature(feature,hstart,hend,hid) values (last_insert_id()," .
                                $homol->start        . "," .
                                $homol->end          . ",\"" .
                                $homol->seqname      . "\")");

    my $rv  = $sth->execute;

    return $rv;
}


=head2 write_Transcript

 Title   : write_Transcript
 Usage   : $obj->write_Transcript($trans,$gene)
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Transcript{
   my ($self,$trans,$gene) = @_;
   

   if( ! $trans->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$trans is not a EnsEMBL transcript - not dumping!");
   }

   if( ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   # ok - now load this line in

   my $tst = $self->prepare("insert into transcript (id,gene,translation) values ('" . $trans->id . "','" . $gene->id . "','" . $trans->translation->id() . "')");
   $tst->execute();
   $self->write_Translation($trans->translation());
   
   return 1;
}

=head2 write_Translation

 Title   : write_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Translation{
   my ($self,$translation) = @_;
   
   if( !$translation->isa('Bio::EnsEMBL::Translation') ) {
       $self->throw("Is not a translation. Cannot write!");
   }

   
   my $tst = $self->prepare("insert into translation (id,seq_start,start_exon,seq_end,end_exon) values ('" 
			    . $translation->id . "',"
			    . $translation->start . ",'"  
			    . $translation->start_exon_id. "',"
			    . $translation->end . ",'"
			    . $translation->end_exon_id . "')");
   $tst->execute();
   
}
   

=head2 write_Exon

 Title   : write_Exon
 Usage   : $obj->write_Exon($exon)
 Function: writes a particular exon into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Exon{
   my ($self,$exon) = @_;
   
   if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
       $self->throw("$exon is not a EnsEMBL exon - not dumping!");
   }

#   my $lockst = $self->prepare("lock exon");
#   $lockst->execute;

   # ok - now load this line in

   # FIXME: better done with placeholders. (perhaps?).

   my $exonst = "insert into exon (id,contig,created,modified,seq_start,seq_end,strand,phase) values ('" .
       $exon->id() . "','" .
	   $exon->contig_id() . "','" .
	       $exon->created(). "','" .
		   $exon->modified . "'," .
		       $exon->start . ",".
			   $exon->end . ",".
			       $exon->strand . ",".
				   $exon->phase . ")";
   
   my $sth = $self->prepare($exonst);
   $sth->execute();
   
#   my $unlockst = $self->prepare("unlock exon");
#   $unlockst->execute;
   
   return 1;
}

=head2 write_Contig

 Title   : write_Contig
 Usage   : $obj->write_Contig($contig,$clone)
 Function: writes a contig and its dna into the database
 Example :
 Returns : 
 Args    :


=cut


sub write_Contig {
   my($self,$contig,$clone)  = @_;

   if( ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI  - can't insert contig for clone $clone");
   }

   my $dna = $contig->seq || $self->throw("No sequence in contig object");
             $dna->id     || $self->throw("No contig id entered.");
             $clone       || $self->throw("No clone entered.");

   my $contigid  = $dna->id;
   my $date      = `date '+%Y-%m-%d'`; chomp $date;
   my $len       = $dna->seq_len;
   my $seqstr    = $dna->seq;
   my $offset    = $contig->offset();
   my $orientation    = $contig->orientation();
   my @sql;

   push(@sql,"lock tables contig write,dna write");
   push(@sql,"insert into dna(contig,sequence,created) values('$contigid','$seqstr','$date')");
   push(@sql,"replace into contig(id,dna,length,clone,offset,orientation) values('$contigid',LAST_INSERT_ID(),$len,'$clone',$offset,$orientation)");
   push(@sql,"unlock tables");   

   foreach my $sql (@sql) {
     my $sth =  $self->prepare($sql);
     my $rv  =  $sth->execute();
     $self->throw("Failed to insert contig $contigid") unless $rv;
   }

   return 1;
}

=head2 write_Clone

 Title   : write_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Clone{
   my ($self,$clone) = @_;

   if( !$clone->isa('Bio::EnsEMBL::DB::CloneI') ) {
       $self->throw("Clone must be a CloneI type, not a $clone");
   }


   my $clone_id = $clone->id();
   my $version = $clone->version();
   my $embl_id = $clone->embl_id();
   my $htg_phase = $clone->htg_phase();
   my @sql;

   push(@sql,"lock tables clone write");
   push(@sql,"insert into clone(id,version,embl_id,htg_phase) values('$clone_id','$version','$embl_id','$htg_phase')");
   push(@sql,"unlock tables");   

   foreach my $sql (@sql) {
     my $sth =  $self->prepare($sql);
     my $rv  =  $sth->execute();
     $self->throw("Failed to insert clone $clone_id") unless $rv;
   }

   foreach my $contig ( $clone->get_all_Contigs() ) {
       $self->write_Contig($contig,$clone_id);
   }

   
}


=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \" \" ");
 Function: prepares a SQL statement on the DBI handle

           If the debug level is greater than 10, provides information into the
           DummyStatement object
 Example :
 Returns : A DBI statement handle object
 Args    : a SQL string


=cut

sub prepare{
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query!");
   }

   if( $self->_debug > 10 ) {
       print STDERR "Prepared statement $string\n";
       my $st = Bio::EnsEMBL::DBSQL::DummyStatement->new();
       $st->_fileh(\*STDERR);
       $st->_statement($string);
       return $st;
   }

   # should we try to verify the string?

   return $self->_db_handle->prepare($string);
}

=head2 _contig_seq_cache

 Title   : _contig_seq_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _contig_seq_cache{
   my ($self,$id,$seq) = @_;

   if( $seq ) {
       
       #
       # Every 100 hits, flush the cache
       #
       if( $self->{'_contig_seq_cnt'} > 100 ) {
	   $self->_flush_seq_cache;
	   $self->{'_contig_seq_cnt'} = 0;
       }

       $self->{'_contig_seq_cnt'}++;
       $self->{'_contig_seq_cache'}->{$id} = $seq;
   }

   return $self->{'_contig_seq_cache'}->{$id};
}

=head2 _flush_seq_cache

 Title   : _flush_seq_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _flush_seq_cache{
   my ($self,@args) = @_;

   $self->{'_contig_seq_cache'} = {};

}

=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_debug'} = $value;
    }
    return $self->{'_debug'};
    
}


=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of _db_handle
 Args    : newvalue (optional)


=cut

sub _db_handle{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}

=head2 _lock_tables

 Title   : _lock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _lock_tables{
   my ($self,@tables) = @_;
   
   my $state;
   foreach my $table ( @tables ) {
       if( $self->{'_lock_table_hash'}->{$table} == 1 ) {
	   $self->warn("$table already locked. Relock request ignored");
       } else {
	   if( $state ) { $state .= ","; } 
	   $state .= "$table write";
	   $self->{'_lock_table_hash'}->{$table} = 1;
       }
   }

   my $sth = $self->prepare("lock tables $state");
   my $rv = $sth->execute();
   $self->throw("Failed to lock tables $state") unless $rv;

}

=head2 _unlock_tables

 Title   : _unlock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _unlock_tables{
   my ($self,@tables) = @_;

   my $sth = $self->prepare("unlock tables");
   my $rv = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY{
   my ($obj) = @_;

   $obj->_unlock_tables();

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}


