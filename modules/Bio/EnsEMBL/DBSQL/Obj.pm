#
#
# BioPerl module for DBSQL::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::Obj - Object representing an instance of an EnsEMBL DB

=head1 SYNOPSIS

    $db = new Bio::EnsEMBL::DBSQL::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );

    $clone  = $db->get_clone('X45667');

    $contig = $db->get_Contig("dJ52N12.02793");

    $gene   = $db->get_Gene('HG45501');

    

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
use Bio::EnsEMBL::Ghost;
use Bio::EnsEMBL::DBSQL::Contig;
use Bio::EnsEMBL::DBSQL::Clone;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Analysis;
use DBI;

use Bio::EnsEMBL::DBSQL::DummyStatement;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($db,$host,$driver,$user,$password,$debug) = 
      $self->_rearrange([qw(DBNAME
			    HOST
			    DRIVER
			    USER
			    PASS
			    DEBUG
			    )],@args);
  $db || $self->throw("Database object must have a database name");
  $user || $self->throw("Database object must have a user");

  #
  # This needs to be rethought. We are caching sequences
  # here to allow multiple exons to be retrieved fine
  # And now more cache's. I think cache's might be a fact of life...
  # 

  $self->{'_contig_seq_cache'} = {};
  $self->_analysis_cache({});
  $self->{'_contig_seq_cnt'} = 0;
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
 Usage   : $obj->get_Gene($geneid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene('ENSG00000009151','evidence')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene{
   my ($self,$geneid, $supporting) = @_;
   my @out;

   if ($supporting && $supporting eq 'evidence') {
       @out = $self->get_Gene_array_supporting('evidence',$geneid);
   }
   else {
       @out = $self->get_Gene_array_supporting('without',$geneid);
   }
   return $out[0];
}


=head2 get_Gene_array

 Title   : get_Gene_array
 Usage   :
 Function: old method, present for historical reasons, points to new method
           get_gene_array_supporting without asking for supp.evidence
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene_array {
    my ($self,@geneid) = @_;

    if( @geneid == 0 ) {
	$self->throw("Attempting to create gene with no id");
    }
    
    my @out=$self->get_Gene_array_supporting('without',@geneid); 

    return @out;
}

=head2 get_Gene_array_supporting

 Title   : get_Gene_array_supporting
 Usage   : $obj->get_Gene_array_supporting($supporting,@geneid)
 Function: Gets an array of genes, with transcripts and exons. If $supporting
           equal to 'evidence' the supporting evidence for each exon is also read
           from the supporting evidence table
 Example : $obj->get_Gene_array_supporting ('evidence',@geneid)
 Returns : an array of gene objects
 Args    : 'evidence' and gene id array


=cut

sub get_Gene_array_supporting {
    my ($self,$supporting,@geneid) = @_;

    $supporting || $self->throw("You need to specify whether to retrieve supporting evidence or not!");

    if( @geneid == 0 ) {
	$self->throw("Attempting to create gene with no id");
    }

    
    my (@out, @sup_exons);
    my $inlist = join(',',map "'$_'", @geneid);
    $inlist = "($inlist)";
				
    # I know this SQL statement is silly.
    #
    
				  
    my $sth = $self->prepare("select p3.gene,p4.id,p3.id,p1.exon,p1.rank,p2.seq_start,p2.seq_end,UNIX_TIMESTAMP(p2.created),UNIX_TIMESTAMP(p2.modified),p2.strand,p2.phase,p5.seq_start,p5.start_exon,p5.seq_end,p5.end_exon,p5.id,p6.version,p3.version,p2.version,p5.version,p4.clone 
                              from gene as p6,contig as p4, transcript as p3, exon_transcript as p1, exon as p2,translation as p5,geneclone_neighbourhood as p7 
                              where p6.id in $inlist and p3.gene = p6.id and p4.clone = p7.clone and p7.gene = p6.id and p2.contig = p4.id and p1.exon = p2.id and p3.id = p1.transcript and p5.id = p3.translation order by p3.gene,p3.id,p1.rank");
    
    $sth->execute();
    
    my $current_gene_id       = '';
    my $current_transcript_id = '';
    
    my ($gene,$trans);
    
    while( (my $arr = $sth->fetchrow_arrayref()) ) {
	my ($geneid,$contigid,$transcriptid,$exonid,$rank,$start,$end,$exoncreated,$exonmodified,$strand,$phase,$trans_start,$trans_exon_start,$trans_end,$trans_exon_end,$translationid,$geneversion,$transcriptversion,$exonversion,$translationversion,$cloneid) = @{$arr};
	
	#print STDERR "Got exon $exonid\n";
	
	if( ! defined $phase ) {
	    $self->throw("Bad internal error! Have not got all the elements in gene array retrieval");
	}
	
	if( $geneid ne $current_gene_id ) {

	    if( $transcriptid eq $current_transcript_id ) {
		$self->throw("Bad internal error. Switching genes without switching transcripts");
	    } 
	    
	    $gene = Bio::EnsEMBL::Gene->new();
	    $gene->id($geneid);
	    
	    #   $sth = $self->_dbobj->prepare("select version from gene where id='".$gene->id."'");
	    #   $sth->execute();
	    #   my $rowhash = $sth->fetchrow_hashref();
	    
	    $gene->version($geneversion);
	    $gene->add_cloneid_neighbourhood($cloneid);
	    $current_gene_id = $geneid;
	    push(@out,$gene);
	    #print STDERR "Made new gene\n";
	}
	
	if( $transcriptid ne $current_transcript_id ) {
	    $trans = Bio::EnsEMBL::Transcript->new();
	    $trans->id($transcriptid);
	    $trans->version($transcriptversion);
	    $current_transcript_id = $transcriptid;
	    
	    my $translation = Bio::EnsEMBL::Translation->new();
	    $translation->start        ($trans_start);
	    $translation->end          ($trans_end);
	    $translation->start_exon_id($trans_exon_start);
	    $translation->end_exon_id  ($trans_exon_end);
	    $translation->id           ($translationid);
	    $translation->version      ($translationversion);
	    $trans->translation        ($translation);
	    $gene ->add_Transcript     ($trans);
	}
	
	my $exon = Bio::EnsEMBL::Exon->new();
	
	$exon->clone_id ($cloneid);
	$exon->contig_id($contigid);
	$exon->id       ($exonid);
	$exon->created  ($exoncreated);
	$exon->modified ($exonmodified);
	$exon->start    ($start);
	$exon->end      ($end);
	$exon->strand   ($strand);
	$exon->phase    ($phase);
	$exon->version  ($exonversion);

	#
	# Attach the sequence, cached if necessary...
	#
	if ($supporting && $supporting eq 'evidence') {
	    push @sup_exons, $exon;
	}
	
	my $seq;
	
	if( $self->_contig_seq_cache($exon->contig_id) ) {
	    $seq = $self->_contig_seq_cache($exon->contig_id);
	} else {
	    my $contig = $self->get_Contig($exon->contig_id());
	    $seq = $contig->seq();
	    $self->_contig_seq_cache($exon->contig_id,$seq);
	}
	
	$exon ->attach_seq($seq);
	$trans->add_Exon($exon);

    }
    
    if ($supporting && $supporting eq 'evidence') {
	#print "get_Gene: getting supporting evidence for array starting with ".@sup_exons[0]->id."\n";
	$self->get_supporting_evidence(@sup_exons);
    }

    return @out;
}


=head2 donor_locator
    
 Title   : get_donor_locator
 Usage   : $obj->get_donor_locator; 
 Function: Reads the meta table of the database to get the donor_database_locator
 Example : get_donor_locator
 Returns : locator string
 Args    : none


=cut

sub get_donor_locator{
    my ($self) = @_;
    
    my $sth = $self->prepare("select donor_database_locator from meta");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $donor = $rowhash->{'donor_database_locator'};
    ($donor eq "") && $self->throw ("No value stored for database locator in meta table!");
    return $rowhash->{'donor_database_locator'};
}

=head2 get_last_update_offset

 Title   : get_last_update_offset
 Usage   : $obj->get_last_update_offset; 
 Function: Reads the meta table of the database to get the last_update time - offset time
 Example : get_last_update_offset
 Returns : UNIX TIME of last update - offset time
 Args    : none

=cut

sub get_last_update_offset{
    my ($self) = @_;
    
    #Get the last update time
    my $sth = $self->prepare("select last_update from meta");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $last = $rowhash->{'last_update'};

    #Now get the offset time from the meta table, which is in time format
    $sth = $self->prepare("select offset_time from meta");
    $sth->execute();
    $rowhash = $sth->fetchrow_hashref();
    my $offset = $rowhash->{'offset_time'};

    #Perform the subtraction in mysql
    $sth = $self->prepare("select DATE_SUB(\"$last\", INTERVAL \"$offset\" HOUR_SECOND)");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    my $datetime = $rowhash->[0];

    $sth = $self->prepare("select UNIX_TIMESTAMP('".$datetime."')");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    my $last_offset = $rowhash->[0];
    ($last_offset eq "") && $self->throw ("No value stored for last_update in meta table!");
    return $last_offset;
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
    
    #Get the last update time
    my $sth = $self->prepare("select last_update from meta");
    my $res = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $last = $rowhash->{'last_update'};

    $sth = $self->prepare("select UNIX_TIMESTAMP('".$last."')");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    $last = $rowhash->[0];
    ($last eq "") && $self->throw ("No value stored for last_update in meta table!");
    return $last;
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
    my $offset = $self->get_offset;

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
    
=head2 get_offset

 Title   : get_offset
 Usage   : $obj->get_offset; 
 Function: Gets the offset time found in the meta table
 Example : get_offset
 Returns : UNIX TIME of offset_time
 Args    : none


=cut

sub get_offset{
    my ($self) = @_;

     #Now get the offset time from the meta table, which is in time format
    my $sth = $self->prepare("select offset_time from meta");
    $sth->execute();
    my $rowhash = $sth->fetchrow_hashref();
    my $offset = $rowhash->{'offset_time'};

    #Trasform the result into unix time and return it
    $sth = $self->prepare("select UNIX_TIMESTAMP('".$offset."')");
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
    my $seen = 0;
    my $trans = Bio::EnsEMBL::Transcript->new();
    # go over each Transcript
    my $sth = $self->prepare("select exon from exon_transcript where transcript = '$transid'");
    
    my $res = $sth->execute();
    while( my $rowhash = $sth->fetchrow_hashref) {
	my $exon = $self->get_Exon($rowhash->{'exon'});
	$trans->add_Exon($exon);
	$seen = 1;
    }
    if( $seen == 0 ) {
	$self->throw("transcript $transid is not present in db");
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

   my $sth = $self->prepare("select version,seq_start,start_exon,seq_end,end_exon from translation where id = '$translation_id'");
   my $res = $sth->execute();
   my $rowhash = $sth->fetchrow_hashref;
   if( !defined $rowhash ) {
       $self->throw("no translation of $translation_id");
   }

   my $out = Bio::EnsEMBL::Translation->new();
   $out->version($rowhash->{'version'});
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

   my $sth = $self->prepare("select id,version,contig,UNIX_TIMESTAMP(created),UNIX_TIMESTAMP(modified),seq_start,seq_end,strand,phase from exon where id = '$exonid'");
   $sth->execute;
   my $rowhash = $sth->fetchrow_hashref;
   if( ! defined $rowhash ) {
       $self->throw("No exon of this id $exonid");
   }

   my $exon = Bio::EnsEMBL::Exon->new();
   $exon->contig_id($rowhash->{'contig'});
   $exon->version($rowhash->{'version'});

   my $contig_id = $exon->contig_id();

   # we have to make another trip to the database to get out the contig to clone mapping.
   my $sth2 = $self->prepare("select clone from contig where id = '$contig_id'");
   $sth2->execute;
   my $rowhash2 = $sth2->fetchrow_hashref;

   $exon->clone_id($rowhash2->{'clone'});


   # rest of the attributes
   $exon->id($rowhash->{'id'});
   $exon->created($rowhash->{'UNIX_TIMESTAMP(created)'});
   $exon->modified($rowhash->{'UNIX_TIMESTAMP(modified)'});
#   print STDERR "Got exon with ",$exon->created," ",$exon->modified,"\n";
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

=head2 get_updated_Objects
    
 Title   : get_updated_Objects
 Usage   : $obj->get_updated_Objects ($recipient_last_update, $recipient_now)
 Function: Gets all the objects that have been updated (i.e.change in 
	   version number) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_updated_Objects (973036800,973090800)
 Returns : database objects (clones and genes)
 Args    : $recipient_last_update, $recipient_now

=cut

sub get_updated_Objects{
    my ($self, $last_offset, $now_offset) = @_;
    
    $last_offset || $self->throw("Attempting to get updated objects without the recipient db last update time");
    $now_offset  || $self->throw("Attempting to get updated objects without the recipient db current time");
    ($last_offset>$now_offset) && $self->throw("Last update more recent than now-offset time, serious trouble");

    #First, let us convert the unix times now_offset and last into mysql times
    my $sth = $self->prepare("select FROM_UNIXTIME(".$last_offset.")");
    $sth->execute();
    my $rowhash = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];
    print STDERR "Last= $last_offset\n";

    $sth = $self->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];
    print STDERR "now: $now_offset\n";

    #First, get all clone ids that have been updated between last and now-offset
    $sth = $self->prepare("select id from clone where stored > '".$last_offset." - 00:30:00' and stored <= '".$now_offset."'");
   
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
    $sth = $self->prepare("select id from gene where stored > '".$last_offset."' and stored <= '".$now_offset."'");
    $sth->execute;

    my @genes;
    while( $rowhash = $sth->fetchrow_hashref) {
	push(@genes,$rowhash->{'id'});
    }
    
    #Get all gene objects for the ids contained in @clones, and push them in @out
    foreach my $geneid (@genes) {
	push @out, $self->get_Gene ($geneid);
    }	
    return @out;
}

=head2 get_updated_Ghosts
    
 Title   : get_updated_Ghosts
 Usage   : $obj->get_updated_Ghosts ($recipient_last_update, $recipient_now_offset)
 Function: Gets all the ghosts for objects that have been deleted (i.e.permanently from 
	   the donor db) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_updated_Ghosts (973036800,973090800)
 Returns : ghost objects
 Args    : $recipient_last_update, $recipient_now_offset

=cut

sub get_updated_Ghosts{
    my ($self, $last_offset, $now_offset) = @_;
    my @out;
    
    $last_offset || $self->throw("Attempting to get updated objects without the recipient db last update time");
    $now_offset  || $self->throw("Attempting to get updated objects without the recipient db current time");
    ($last_offset>$now_offset) && $self->throw("Last update more recent than now-offset time, serious trouble");
    
    #First, let us convert the unix times now_offset and last into mysql times
    my $sth = $self->prepare("select FROM_UNIXTIME(".$last_offset.")");
    $sth->execute();
    my $rowhash = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];
    
    $sth = $self->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];
    
    #Get all ghosts that have deleted times between last and now-offset
    $sth = $self->prepare("select id,version,obj_type,deleted,stored from ghost where stored > '".$last_offset."' and stored <= '".$now_offset."'");
    $sth->execute();
    while(my $rowhash = $sth->fetchrow_hashref()) {
	my $ghost = Bio::EnsEMBL::Ghost->new();
	$ghost->id($rowhash->{'id'});
	$ghost->version($rowhash->{'version'});
	$ghost->obj_type($rowhash->{'obj_type'});
	$ghost->deleted($rowhash->{'deleted'});
	$ghost->_stored($rowhash->{'stored'});
	push @out, $ghost;
    }
    return @out;
}

=head2 get_Ghost
    
 Title   : get_Ghost
 Usage   : $obj->get_Ghost ($ghost_id,$ghost_version,$ghost_obj_type)
 Function: Gets a ghost by id, version,obj_type  
 Example : $obj->get_Ghost ('test','1','transcript')
 Returns : ghost objects
 Args    : ghost id, version and object type

=cut

sub get_Ghost{
    my ($self, $g_id, $g_obj_type) = @_;
    my @out;
    
    $g_id || $self->throw("Attempting to get a ghost object without an id");
    $g_obj_type || $self->throw("Attempting to get a ghost object without an object type");

    my $sth = $self->prepare("select id,version,obj_type,deleted,stored from ghost where id='".$g_id."' and obj_type = '".$g_obj_type."'");
    $sth->execute();
    my  $rv = $sth->rows;
    ! $rv && $self->throw("Ghost not found in database!");
    my $rowhash = $sth->fetchrow_hashref();
    my $ghost = Bio::EnsEMBL::Ghost->new();
    $ghost->id($rowhash->{'id'});
    $ghost->version($rowhash->{'version'});
    $ghost->obj_type($rowhash->{'obj_type'});
    $ghost->deleted($rowhash->{'deleted'});
    $ghost->_stored($rowhash->{'stored'});
    return $ghost;
}

=head2 write_Ghost
    
 Title   : write_Ghost
 Usage   : $obj->write_Ghost ($ghost)
 Function: Writes a ghost to the database  
 Example : $obj->write_Ghost ($ghost)
 Returns : 
 Args    : ghost object

=cut

sub write_Ghost{
    my ($self, $ghost) = @_;
    
    $ghost || $self->throw("Attempting to write a ghost without a ghost object");
    $ghost->isa("Bio::EnsEMBL::Ghost") || $self->throw("$ghost is not an EnsEMBL ghost - not dumping!");
    
    my $sth = $self->prepare("insert into ghost (id, version, obj_type,deleted,stored) values('".$ghost->id."','".$ghost->version."','".$ghost->obj_type."','".$ghost->deleted."',now())");
    $sth->execute();
    return 1;
}

=head2 archive_Gene
    
 Title   : archive_Gene
 Usage   : $obj->archive_gene($gene,$clone,$arcdb)
 Function: Deletes a gene and all its transcripts and exons, 
           and archives partial info in the archive db passed on.
 Example : 
 Returns : nothing
 Args    : $gene, $arcdb (archive database object)


=cut

sub archive_Gene {
   my ($self,$gene,$arc_db) = @_;
   my $sth;

   # get transcripts for the gene given 
   
   foreach my $transcript ($gene->each_Transcript) {
       
       #Get out transcript info needed to write into archive db
       my $seq = $transcript->dna_seq;
       $seq->id($transcript->id);
       
       #Temporary, since versions not stored yet...
       !$transcript->version && $transcript->version(1);
       !$gene->version && $gene->version(1);
       
       #Finally, write all the info to a new entry in the archive database
       
       $arc_db->write_seq($seq, $transcript->version, 'transcript', $gene->id, $gene->version);
       
       #Get out translation to write protein into archive db
       #Note: version is the one from transcript!

       $seq = $transcript->translate;

       $arc_db->write_seq($seq, $transcript->version, 'protein', $gene->id, $gene->version);

       #Delete transcript rows
       $sth= $self->prepare("delete from transcript where id = '".$transcript->id."'");
       $sth->execute;
       
       #First get all exons of the gene using $gene->each_unique_Exon method, to write them
       foreach my $exon ($gene->each_unique_Exon) {
	   #Get out info needed to write into archive db
	   $seq = $exon->seq;
	   $seq->id($exon->id);
	   
	   #Temporary, since versions not stored yet...
	   !$exon->version && $exon->version(1);
	   
	   #Write into archive db
	   $arc_db->write_seq($seq, $exon->version, 'exon', $gene->id, $gene->version);
       }

       #Then use the each_Exon method from transcript to delete exons
       foreach my $exon ($transcript->each_Exon) {
	   
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

=head2 delete_Exon

 Title   : delete_Clone
 Usage   : $obj->delete_Exon($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->delete_Exon(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Exon{
    my ($self,$exon_id) = @_;
    $exon_id || $self->throw ("Trying to delete an exon without an exon_id\n");
    
    #Delete exon_transcript rows
    my $sth= $self->prepare("delete from exon_transcript where transcript = '".$exon_id."'");
    $sth->execute;
    #Delete exon rows
    $sth = $self->prepare("delete from exon where id = '".$exon_id."'");
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
   my ($self,$clone_id) = @_;
   
   $clone_id || $self->throw ("Trying to delete clone without a clone_id\n");
   
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
   my $translation;
   # get out exons, transcripts for gene. 

   my $sth = $self->prepare("select id,translation from transcript where gene = '$geneid'");
   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@trans,$rowhash->{'id'});
       $translation = $rowhash->{'translation'};
   }

   foreach my $trans ( @trans ) {
       my $sth = $self->prepare("select exon from exon_transcript where transcript = '$trans'");
       $sth->execute;
       while( my $rowhash = $sth->fetchrow_hashref) {
	   $exon{$rowhash->{'exon'}} =1;
       }
   }

   my $sth2 = $self->prepare("delete from translation where id = '$translation'");
   $sth2->execute;

   # delete exons, transcripts, gene rows

   foreach my $exon ( keys %exon ) {
       print STDERR "Got exon $exon in delete_Gene\n";
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

   my $sth = $self->prepare("select clone from geneclone_neighbourhood where gene = '$geneid'");

   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'clone'});
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

   my $sth = $self->prepare("select count(*),cont.clone ,ex.contig,tran.gene  from contig as cont, transcript as tran, exon_transcript as et, exon as ex where ex.id = et.exon and tran.id = et.transcript and cont.clone = '$cloneid'  and cont.id = ex.contig group by tran.gene");

   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'clone'});
   }

   return @out;
}

=head2 replace_last_update
    
 Title   : replace_last_update(@$now_offset)
 Usage   : $obj->replace_last_update($now_offset)
 Function: Replaces the time in the last update field of the meta table with the now_offset time of the recipient
 Example : 
 Returns : nothing
 Args    : 

=cut

sub replace_last_update {
    my ($self, $now_offset) = @_;
    
    $now_offset || $self->throw("Trying to replace last update without a now-offset time\n");
    
    my $last_offset= $self->get_last_update;
    
    my $sth = $self->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $sth->execute();
    my $rowhash = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];
    
    $sth = $self->prepare("select FROM_UNIXTIME(".$last_offset.")");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];
    
    my $donor=$self->get_donor_locator;
    my $offset=$self->get_offset;
    
    $sth = $self->prepare("delete from meta where last_update = '".$last_offset."'");  
    $sth->execute;
 
    $sth = $self->prepare("insert into meta (last_update,donor_database_locator,offset_time) values ('".$now_offset."','".$donor."','".$offset."')");
    $sth->execute;
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
   my $old_gene;

   my %done;

   if( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
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

   $old_gene = $self->get_Gene($gene->id);
   
   if ( !defined($old_gene) || ($gene->version > $old_gene->version)) {

       !$gene->created() && $gene->created(0);
       !$gene->modified() && $gene->modified(0);

       my $sth2 = $self->prepare("insert into gene (id,version,created,modified,stored) values ('". 
				 $gene->id       . "','".
				 $gene->version  . "',FROM_UNIXTIME(".
				 $gene->created  . "),FROM_UNIXTIME(".
				 $gene->modified . "),now())");
       $sth2->execute();
   }
   else {
       print "Got this gene with this version already, no need to write in db\n";
   }

   foreach my $cloneid ($gene->each_cloneid_neighbourhood) {

       my $sth = $self->prepare("select gene,clone from geneclone_neighbourhood where gene='".$gene->id."' && clone='$cloneid'");
       $sth->execute();
       my $rowhash =  $sth->fetchrow_arrayref();
       my  $rv = $sth->rows;
       if( ! $rv ) {
	   $sth = $self->prepare("insert into geneclone_neighbourhood (gene,clone) values ('" . 
				 $gene->id . "','". 
				 $cloneid ."')");
	   $sth->execute();
       }
   }
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
 Usage   : $obj->write_Feature($contig,@features)
 Function: Writes a feature on the genomic sequence of a contig into the database
 Example :
 Returns : nothing
 Args    : Bio::EnsEMBL::SeqFeatureI


=cut

sub write_Feature {
    my ($self,$contig,@features) = @_;

    #
    # Yes - we need to rethink how we are writing features into the
    # database. This is a little obtuse and clunky now
    #

    $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI")          unless (defined($contig) && $contig->isa("Bio::EnsEMBL::DB::ContigI"));
    
    my $contigid = $contig->id;
    my $analysis;
    
    my $sth = $self->prepare("insert into feature(id,contig,seq_start,seq_end,score,strand,name,analysis,hstart,hend,hid) values (?,?,?,?,?,?,?,?,?,?,?)");
    
    # Put the repeats in a different table, and also things we need to write
    # as fsets.
    my @repeats;
    my @fset;

    FEATURE :
    foreach my $feature ( @features ) {
	
	if( ! $feature->isa('Bio::EnsEMBL::SeqFeatureI') ) {
	    $self->throw("Feature $feature is not a feature!");
	}

	eval {
	    $feature->validate();
	};

	if ($@) {
	    print(STDERR "Feature invalid. Skipping feature\n");
	    next FEATURE;
	}

	
	if($feature->isa('Bio::EnsEMBL::Repeat')) {
	    push(@repeats,$feature);
	} elsif ( $feature->sub_SeqFeature ) {
	    push(@fset,$feature);
	} else {    
	    if (!defined($feature->analysis)) {
		$self->throw("Feature " . $feature->seqname . " " . $feature->source_tag ." doesn't have analysis. Can't write to database");
	    } else {
		$analysis = $feature->analysis;
	    }
	    
	    my $analysisid = $self->write_Analysis($analysis);


	    if ( $feature->isa('Bio::EnsEMBL::FeaturePair') ) {
		my $homol = $feature->feature2;
	    
		$sth->execute('NULL',
			      $contig->id,
			      $feature->start,
			      $feature->end,
			      $feature->score,
			      $feature->strand,
			      $feature->source_tag,
			      $analysisid,
			      $homol->start,
			      $homol->end,
			      $homol->seqname);
	    } else {
		$sth->execute('NULL',
			      $contig->id,
			      $feature->start,
			      $feature->end,
			      $feature->score,
			      $feature->strand,
			      $feature->source_tag,
			      $analysisid,
			      -1,
			      -1,
			      "__NONE__");
	    }
	}
    }

    my $sth2 = $self->prepare("insert into repeat_feature(id,contig,seq_start,seq_end,score,strand,analysis,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?)");

    foreach my $feature (@repeats) {
	if (!defined($feature->analysis)) {
	    $self->throw("Feature " . $feature->seqname . " " . $feature->source_tag ." doesn't have analysis. Can't write to database");
	} else {
	    $analysis = $feature->analysis;
	}

	my $analysisid = $self->write_Analysis($analysis);
	my $homol = $feature->feature2;

	$sth2->execute('NULL',
		       $contig->id,
		       $feature->start,
		       $feature->end,
		       $feature->score,
		       $feature->strand,
		       $analysisid,
		       $homol->start,
		       $homol->end,
		       $homol->seqname);
    }


    # Now the predictions
    # we can't block do these as we need to get out the id wrt to the features
    foreach my $feature ( @fset ) {
#	print STDERR "Adding in a fset feature ",$feature->gff_string,"\n";

	if (!defined($feature->analysis)) {

	    $self->throw("Feature " . $feature->seqname . " " . 
			              $feature->source_tag . 
			 " doesn't have analysis. Can't write to database");
	} else {
	    $analysis = $feature->analysis;
	}

	my $analysisid = $self->write_Analysis($analysis);
	my $score      = $feature->score();

	if( !defined $score ) { $score = "-1000"; }

	my $sth3 = $self->prepare("insert into fset(id,score) values ('NULL',$score)");
	   $sth3->execute();

	# get out this id. This looks really clunk I know. Any better ideas... ?

	my $sth4 = $self->prepare("select LAST_INSERT_ID()");
	   $sth4->execute();
	my $arr = $sth4->fetchrow_arrayref();
	my $fset_id = $arr->[0];

	# now write each sub feature
	my $rank = 1;

	foreach my $sub ( $feature->sub_SeqFeature ) {
	    my $sth5 = $self->prepare("insert into feature(id,contig,seq_start,seq_end,score,strand,analysis,hstart,hend,hid) values('NULL','".$contig->id."',"
				      .$sub->start   .","
				      .$sub->end     . ","
				      .$sub->score   . ","
				      .$sub->strand  . ","
				      .$analysisid   . ",-1,-1,'__NONE__')");
	    $sth5->execute();
	    my $sth6 = $self->prepare("insert into fset_feature(fset,feature,rank) values ($fset_id,LAST_INSERT_ID(),$rank)");
	    $sth6->execute();
	    $rank++;
	}
    }
    
    return 1;
}

=head2 write_supporting_evidence

 Title   : write_supporting_evidence
 Usage   : $obj->write_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : None


=cut

sub write_supporting_evidence {
    my ($self,$exon) = @_;

    $self->throw("Argument must be Bio::EnsEMBL::Exon. You entered [$exon]\n") unless $exon->isa("Bio::EnsEMBL::Exon");

     my $sth  = $self->prepare("insert into supporting_feature(id,exon,seq_start,seq_end,score,strand,analysis,name,hstart,hend,hid) values(?,?,?,?,?,?,?,?,?,?,?)");
    
    FEATURE: foreach my $f ($exon->each_Supporting_Feature) {

	eval {
	    $f->validate();
	};

	if ($@) {
	    print(STDERR "Supporting feature invalid. Skipping feature\n");
	    next FEATURE;
	}

	my $analysisid = $self->write_Analysis($f->analysis);
	
	if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
	    $sth->execute('NULL',
			  $exon->id,
			  $f->start,
			  $f->end,
			  $f->score,
			  $f->strand,
			  $analysisid,
			  $f->source_tag,
			  $f->hstart,
			  $f->hend,
			  $f->hseqname
			  );
	} else {
	    $self->warn("Feature is not a Bio::EnsEMBL::FeaturePair");
	}
    }
}


=head2 get_supporting_evidence

 Title   : get_supporting_evidence
 Usage   : $obj->get_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : array of exon objects, needed to know which exon to attach the evidence to


=cut

sub get_supporting_evidence {
    my ($self,@exons) = @_;

    my $instring = "'";
    my %exhash;

    if (@exons == 0) {
	$self->throw("No exon objects were passed on!");
    }

    foreach my $exon (@exons) {

	$exhash{$exon->id} = $exon;

	$instring = $instring . $exon->id . "','";
    }
    
    $instring = substr($instring,0,-2);
   
    my $sth = $self->prepare("select * from supporting_feature where exon in (" . $instring . ")");
    $sth->execute;

    my %anahash;

    while (my $rowhash = $sth->fetchrow_hashref) {
	my $f1 = new Bio::EnsEMBL::SeqFeature;
	my $f2 = new Bio::EnsEMBL::SeqFeature;
	
	my $f = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					      -feature2 => $f2);

	my $exon = $rowhash->{exon};

#	$f1->seqname($rowhash->{contig});
	$f1->seqname("Supporting_feature");
	$f1->start  ($rowhash->{seq_start});
	$f1->end    ($rowhash->{seq_end});
	$f1->strand ($rowhash->{strand});
	$f1->source_tag($rowhash->{name});
	$f1->primary_tag('similarity');
	$f1->score  ($rowhash->{score});
	
	$f2->seqname($rowhash->{hid});
	$f2->start  ($rowhash->{hstart});
	$f2->end    ($rowhash->{hend});
	$f2->strand ($rowhash->{strand});
	$f2->source_tag($rowhash->{name});
	$f2->primary_tag('similarity');
	$f2->score  ($rowhash->{score});

	my $analysisid = $rowhash->{analysis};

	if ($anahash{$analysisid}) {
	    $f->analysis($anahash{$analysisid});

	} else {
	    $f->analysis($self->get_Analysis($analysisid));
	    $anahash{$analysisid} = $f->analysis;
	}
	
	$f->validate;

	$exhash{$exon}->add_Supporting_Feature($f);
    }

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
    my $analysisid = $self->exists_Analysis($anal);

    return $analysisid if $analysisid;

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

    my $sth  = $self->prepare($query);
    my $rv   = $sth->execute;
    
    
    $sth = $self->prepare("select last_insert_id()");
    $rv  = $sth->execute;
    
    if ($sth->rows == 1) {
	my $rowhash = $sth->fetchrow_hashref;
	return $rowhash->{'last_insert_id()'};
    } else {
	$self->throw("Wrong number of rows returned : " . $sth->rows . " : should be 1");
    }

}


=head2 exists_Homol_Feature

 Title   : exists_Homol_Feature
 Usage   : $obj->exists_Homol_Feature($feature)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : nothing
 Args    : Bio::SeqFeature::Homol


=cut

sub exists_Homol_Feature {
    my ($self,$feature,$analysisid,$contig) = @_;

    $self->throw("Feature is not a Bio::SeqFeature::Homol") unless $feature->isa("Bio::SeqFeature::Homol");
    
    my $homol = $feature->homol_SeqFeature;

    if (!defined($homol)) {
	$self->throw("Homol feature doesn't exist");
    }
    
    my $query = "select f.id  from feature as f,homol_feature as h where " .
			     "  f.id = h.feature " . 
			     "  and h.hstart = "    . $homol->start . 
			     "  and h.hend   = "    . $homol->end   . 
			     "  and h.hid    = '"   . $homol->seqname . 
			     "' and f.contig = '"   . $contig->id . 
			     "' and f.seq_start = " . $feature->start . 
			     "  and f.seq_end = "   . $feature->end .
			     "  and f.score = "     . $feature->score . 
			     "  and f.name = '"     . $feature->source_tag . 
			     "' and f.analysis = "  . $analysisid;


    my $sth = $self->prepare($query);
    my $rv  = $sth->execute;
    my $rowhash;

    if ($rv && $sth->rows > 0) {
	my $rowhash = $sth->fetchrow_hashref;
	return $rowhash->{'id'};
    } else {
	return 0;
    }
}
    

=head2 get_Analysis

 Title   : get_Analysis
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Analysis {
    my ($self,$id) = @_;

    my $sth = $self->prepare("select * from analysis where id = $id");
    my $rv  = $sth->execute;
    my $rh  = $sth->fetchrow_hashref;

    if ($sth->rows) {

	my $anal = new Bio::EnsEMBL::Analysis::Analysis(-db    => $rh->{db},
					      -db_version      => $rh->{db_version},
					      -program         => $rh->{program},
					      -program_version => $rh->{program_version},
					      -gff_source      => $rh->{gff_source},
					      -gff_feature     => $rh->{gff_feature},
					      -id              => $rh->{id},
					      );

	return $anal;
    }  else {
	$self->throw("Can't fetch analysis id $id\n");
    }
    
}


=head2 exists_Analysis

 Title   : get_Analysis
 Usage   : $obj->exists_Analysis($anal)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : Analysis id if the entry exists
 Args    : Bio::EnsEMBL::Analysis::Analysis


=cut

sub exists_Analysis {
    my ($self,$anal) = @_;

    $self->throw("Object is not a Bio::EnsEMBL::Analysis::Analysis") unless $anal->isa("Bio::EnsEMBL::Analysis::Analysis");

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
    
    if( exists $self->_analysis_cache->{$query} ) {
	return $self->_analysis_cache->{$query};
    }

    my $sth = $self->prepare($query);
    

    my $rv = $sth->execute();


    if ($rv && $sth->rows > 0) {
	my $rowhash = $sth->fetchrow_hashref;
	my $anaid = $rowhash->{'id'}; 
	$self->_analysis_cache->{$query} = $anaid;
	return $anaid;
    } else {
	return 0;
    }
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
   my $old_trans;

   if( ! $trans->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$trans is not a EnsEMBL transcript - not dumping!");
   }

   if( ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   # ok - now load this line in

   
   eval {
       $old_trans=$self->get_Transcript($trans->id);
   };
    

   if ( $@ || ($trans->version > $old_trans->version)) {

       my $tst = $self->prepare("insert into transcript (id,gene,translation,version) values ('" . $trans->id . "','" . $gene->id . "','" . $trans->translation->id() . "',".$trans->version.")");
       $tst->execute();
       $self->write_Translation($trans->translation());
   }
   else {
       print "Transcript already present in the database with the same version number $old_trans [",$old_trans->version,"], no need to write it in\n";
   }
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
    my $old_transl;
    
    if( !$translation->isa('Bio::EnsEMBL::Translation') ) {
	$self->throw("Is not a translation. Cannot write!");
    }
    
    eval {
	$old_transl=$self->get_Translation($translation->id);
    };
    

    if ( $@ || ($translation->version > $old_transl->version)) {
	if( !defined $translation->version  ) {
	    $self->throw("No version number on translation");
	}
	my $tst = $self->prepare("insert into translation (id,version,seq_start,start_exon,seq_end,end_exon) values ('" 
				 . $translation->id . "',"
				 . $translation->version . ","
				 . $translation->start . ",'"  
				 . $translation->start_exon_id. "',"
				 . $translation->end . ",'"
				 . $translation->end_exon_id . "')");
	$tst->execute();
    }
    else {
	print "Translation already present in the database with the same version number, no need to write it in\n";
    }
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
   my $old_exon;

   if( ! $exon->isa('Bio::EnsEMBL::Exon') ) {
       $self->throw("$exon is not a EnsEMBL exon - not dumping!");
   }

#   my $lockst = $self->prepare("lock exon");
#   $lockst->execute;

   # ok - now load this line in

   # FIXME: better done with placeholders. (perhaps?).

   eval {
       $old_exon=$self->get_Exon($exon->id);
   };
   
   if  ( $@ || ($exon->version > $old_exon->version)) {
#       print(STDERR "Inserting " . $exon->created . " " . $exon->modified . "\n");
       my $exonst = "insert into exon (id,version,contig,created,modified,seq_start,seq_end,strand,phase,stored,end_phase) values ('" .
	   $exon->id() . "'," .
	   $exon->version() . ",'".
	   $exon->contig_id() . "', FROM_UNIXTIME(" .
	   $exon->created(). "), FROM_UNIXTIME(" .
	   $exon->modified() . ")," .
			       $exon->start . ",".
				   $exon->end . ",".
				       $exon->strand . ",".
					   $exon->phase . ",now(),".
					       $exon->end_phase . ")";
       
       my $sth = $self->prepare($exonst);
       $sth->execute();

       # Now the supporting evidence

       $self->write_supporting_evidence($exon);
   }
   else {
       print "Exon with the same version number already present, no need to write it in";
   }
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
   my $date      = $contig->seq_date;
   my $len       = $dna->seq_len;
   my $seqstr    = $dna->seq;
   my $offset    = $contig->offset();
   my $orientation    = $contig->orientation();
   my $order = $contig->order();
   my @sql;

   push(@sql,"lock tables contig write,dna write");
   push(@sql,"insert into dna(contig,sequence,created) values('$contigid','$seqstr',FROM_UNIXTIME($date))");
   push(@sql,"replace into contig(id,dna,length,clone,offset,orientation,corder) values('$contigid',LAST_INSERT_ID(),$len,'$clone',$offset,$orientation,$order)");
   push(@sql,"unlock tables");   

   foreach my $sql (@sql) {
     my $sth =  $self->prepare($sql);
     my $rv  =  $sth->execute();
     $self->throw("Failed to insert contig $contigid") unless $rv;
   }

   # write sequence features. We write all of them together as it
   # is more efficient

   $self->write_Feature($contig,$contig->get_all_SeqFeatures);

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

   my $clone_id = $clone->id;
   $clone || $self->throw("Trying to write a clone without a clone object!\n");
   if( !$clone->isa('Bio::EnsEMBL::DB::CloneI') ) {
       $self->throw("Clone must be a CloneI type, not a $clone");
   }
   
   my @sql;

   push(@sql,"lock tables clone write");
   push(@sql,"insert into clone(id,version,embl_id,embl_version,htg_phase,created,modified,stored) values('$clone_id','".$clone->version."','".$clone->embl_id."','".$clone->embl_version."','".$clone->htg_phase."',FROM_UNIXTIME(".$clone->created."),FROM_UNIXTIME(".$clone->modified."),now())");
   push(@sql,"unlock tables");   

   foreach my $sql (@sql) {
     my $sth =  $self->prepare($sql);
     #print STDERR "Executing $sql\n";
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

=head2 _analysis_cache

 Title   : _analysis_cache
 Usage   : $obj->_analysis_cache()
 Function: 
 Returns : reference to a hash
 Args    : newvalue (optional)


=cut

sub _analysis_cache{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_analysis_cache'} = $value;
    }
    return $obj->{'_analysis_cache'};

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
   my $rv  = $sth->execute();
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



