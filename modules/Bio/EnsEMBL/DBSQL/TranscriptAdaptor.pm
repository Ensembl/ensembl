# EnsEMBL Transcript reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# based on 
# Elia Stupkas Gene_Obj
# 
# Date : 20.02.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::TranscriptAdaptor - MySQL Database queries to generate and store transcripts/translations.

=head1 SYNOPSIS

Transcripts and Translations are stored and fetched in this
object. Translations never go alone any more. The database only
accepts them (at the moment) in a transcript.  

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Elia Stupka  : elia@ebi.ac.uk
  Ewan Birney  : 

=head1 APPENDIX

=cut

;

package Bio::EnsEMBL::DBSQL::TranscriptAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );



=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $transcriptobj->fetch_by_dbID( $dbid )
 Function: 
 Example : $obj->get( .. )
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!

=cut


sub fetch_by_dbID {
    my ($self,$transcriptId) = @_;
}




=head2 _store_exons_in_transcript

 Title   : _store_exons_in_transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _store_exons_in_transcript{
   my ($self,$trans,@exons) = @_;

   if( !ref $trans || !$trans->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw(" $trans is not a transcript");
   }
   #print STDERR "Got ",scalar(@exons),"to store...\n";

   my $exon;
   while ( ($exon = shift @exons)) {
       #print STDERR "Handling exon",$exon->id,":",$exon->sticky_rank,"\n";

       if( $#exons >= 0 && $exons[0]->id eq $exon->id ) {
        
	   # sticky exons.
	   my @sticky_exons;
	   push(@sticky_exons,$exon);
	   while( my $newexon = shift @exons ) {
	       if( $newexon->id eq $exon->id ) {
                            
		   push(@sticky_exons,$newexon);
                   
	       } else {
               
		   unshift(@exons,$exon);
		   last;
	       }
	   }
           
	   my $sticky = $self->_make_sticky_exon(@sticky_exons);
	   #print STDERR "Added sticky exon... $sticky\n";
	   $trans->add_Exon($sticky);
           
       } else {
	   $trans->add_Exon($exon);
       }
   }

}


=head2 get_Transcript
    
 Title   : get_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub fetch_by_dbID {
    my ($self,$transcriptId) = @_;

    my $seen = 0;
    my $trans = Bio::EnsEMBL::Transcript->new();
    my $exonAdaptor = $self->db->get_exonAdaptor();

    my $sth = $self->prepare("select exon from exon_transcript where transcript = '$transid'");
    $sth->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $exon = $exonAdaptor->fetch_by_dbID($rowhash->{'exon'});
	$trans->add_Exon($exon);
	$seen = 1;
    }
    $sth = $self->prepare("select version,translation from transcript where id = '$transid'");
    $sth->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	my $translation = $self->get_Translation($rowhash->{'translation'});
	$trans->translation($translation);
	$trans->version($rowhash->{'version'});
    }
    if ($seen == 0 ) {
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

   my $sth     = $self->_db_obj->prepare("select version,seq_start,start_exon,seq_end,end_exon from translation where id = '$translation_id'");
   my $res     = $sth->execute();
   my $rowhash = $sth->fetchrow_hashref;

   if( !defined $rowhash ) {
       $self->throw("no translation of $translation_id");
   }

   my $out = Bio::EnsEMBL::Translation->new();

   $out->version      ($rowhash->{'version'});
   $out->start        ($rowhash->{'seq_start'});
   $out->end          ($rowhash->{'seq_end'});
   $out->start_exon_id($rowhash->{'start_exon'});
   $out->end_exon_id  ($rowhash->{'end_exon'});
   $out->id           ($translation_id);

   return $out;
}



=head2 get_Virtual_Contig
    
 Title   : get_Virtual_Contig
 Usage   : $gene_obj->get_Virtual_Contig($transcript,$max_length)
 Function: Gets a Bio::EnsEMBL::DB::Virtual Contig object which 
           spans the whole sequence on which the given 
           Bio::EnsEMBL::Transcript object lies, as long 
           as its length does not exceed max_length. If max_length
           is exceeded, undef is returned instead.
 Example : $gene_obj->get_Virtual_Contig($transcript,50000)
 Returns : VirtualContig Object (or undef)
 Args    : Bio::EnsEMBL::Transcript object and max_length int variable


=cut


# coord change to context Chromosome
# get_context ...    
sub get_Virtual_Contig{
}




=head2 get_Transcript_in_VC_coordinates
    
 Title   : get_Transcript_in_VC_coordinates
 Usage   : $gene_obj->get_Transcript_in_VC_coordinates($transcript_id)
 Function: Gets a Bio::EnsEMBL::Transcript object in vc coordinates
 Example : $gene_obj->get_Virtual_Contig($transcript_id)
 Returns : Bio::EnsEMBL::Transcript
 Args    : transcript id


=cut




sub get_Transcript_in_VC_coordinates
{
    
    
}



=head2 store

 Title   : store
 Usage   : $transcriptAdaptor->store( $transcript )
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    : needs a gene ...


=cut

sub store {
   my ($self,$trans,$gene) = @_;
   my $old_trans;

   if( ! $trans->isa('Bio::EnsEMBL::Transcript') ) {
       $self->throw("$trans is not a EnsEMBL transcript - not dumping!");
   }

   if( ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not dumping!");
   }

   # ok - now load this line in
   my $tst = $self->prepare("
        insert into transcript (id, gene, translation, version) 
        values (?, ?, ?, ?)
        ");
                
   $tst->execute(
        $trans->id,
        $gene->id, 
        $trans->translation->id,
        $trans->version   
        );

   #print STDERR "Going to look at gene links\n";

   foreach my $dbl ( $trans->each_DBLink ) {
       #print STDERR "Going to insert for",$trans->id," ",$dbl->primary_id," ",$dbl->database,"\n";
       my $sth3 = $self->_db_obj->prepare("insert into transcriptdblink (transcript_id,external_id,external_db) values ('". 
					  $trans->id        . "','".
					  $dbl->primary_id . "','".
					  $dbl->database   . "')");
       $sth3->execute();
       
   }
  # write the translation

   if( !$translation->isa('Bio::EnsEMBL::Translation') ) {
     $self->throw("Is not a translation. Cannot write!");
   }
   
   if ( !defined $translation->version  ) {
     $self->throw("No version number on translation");
   }
    
   $tst = $self->prepare("insert into translation (id,version,seq_start,start_exon,seq_end,end_exon) values ('" 
			    . $translation->id . "',"
			    . $translation->version . ","
			    . $translation->start . ",'"  
			    . $translation->start_exon_id. "',"
			    . $translation->end . ",'"
			    . $translation->end_exon_id . "')");
   $tst->execute();
   return 1;
}


sub get_New_external_id {
    my ($self,$table,$stub,$number) = @_;

    $table .= "_external";
    if( !defined $number ) {
	$number = 1;
    }

    my @out;


    my $lsth   = $self->_db_obj->prepare("lock table $table write");
    $lsth->execute;

    # wrap critical region in an eval so we can catch errors and release table

    eval {

	my $query = "select max(external_id) as id from $table where external_id like '$stub%'";
	
	my $sth   = $self->_db_obj->prepare($query);
	my $res   = $sth->execute;
	my $row   = $sth->fetchrow_hashref;
	my $id    = $row->{id};
	
	if (!defined($id) || $id eq "") {
	    $id = $stub . "00000000000";
	}
	
	if ($id =~ /\D+(\d+)$/) {
	    
	    my $newid  = $1;
	    my $i;
	    
	    foreach $i ( 1..$number ) {

		$newid++;
		
		
		if (length($newid) > 11) {
		    if ($newid =~ /^0/) {
			$newid =~ s/^0//;
		    } else {
			$self->throw("Can't truncate number string to generate new id [$newid]");
		    }
		}
		my $c = $stub . $newid;
		my $query = "insert into $table (internal_id,external_id) values (NULL,'$c')";
		my $sth   = $self->_db_obj->prepare($query);
		my $res   = $sth->execute;
		
		push(@out,$c);
	    }
	    
	    
	} else {
	    $self->throw("[$id] does not look like an object id (e.g. ENST00000019784)");
	}
    };

    my $error = undef;

    if( $@ ) {
	$error = $@;
    }


    my $usth   = $self->_db_obj->prepare("unlock tables");
    $usth->execute;


    if( defined $error ) {
	$self->throw("Problem in making IDs. Unlocked tables. \n\n Error $@");
    }

    return @out;
    
}

sub get_NewId {
    my ($self,$table,$stub) = @_;

    my $query = "select max(id) as id from $table where id like '$stub%'";
    my $sth   = $self->_db_obj->prepare($query);
    my $res   = $sth->execute;
    my $row   = $sth->fetchrow_hashref;
    my $id    = $row->{id};

    if (!defined($id) || $id eq "") {
	$id = $stub . "00000000000";
    }

    if ($id =~ /\D+(\d+)$/) {

	my $newid  = $1;

	$newid++;
	

	if (length($newid) > 11) {
	    if ($newid =~ /^0/) {
		$newid =~ s/^0//;
	    } else {
		$self->throw("Can't truncate number string to generate new id [$newid]");
	    }
	}
	$newid = $stub . $newid;

	return $newid;
    } else {
	$self->throw("[$id] does not look like an object id (e.g. ENST00000019784)");
    }
    

}



=head2 get_new_GeneID

 Title   : get_new_GeneID
 Usage   : my $id = $geneobj->get_new_GeneID
 Function: 
 Example : 
 Returns : Gets the next unused gene id from the database
 Args    : none


=cut

sub get_new_GeneID {
    my ($self,$stub) = @_;

    $stub = "ENSG" unless defined($stub);

    return $self->get_NewId("gene",$stub);
    
}

=head2 get_new_TranscriptID

 Title   : get_new_TranscriptID
 Usage   : my $id = $geneobj->get_new_TranscriptID
 Function: 
 Example : 
 Returns : Gets the next unused transcript id from the database
 Args    : none


=cut

sub get_new_TranscriptID {
    my ($self,$stub) = @_;

    $stub = "ENST" unless defined($stub);

    return $self->get_NewId("transcript",$stub);

}

=head2 get_new_ExonID

 Title   : get_new_ExonID
 Usage   : my $id = $geneobj->get_new_ExonID
 Function: 
 Example : 
 Returns : Gets the next unused exon id from the database
 Args    : none


=cut

sub get_new_ExonID {
    my ($self,$stub) = @_;

    $stub = "ENSE" unless defined($stub);

    return $self->get_NewId("exon",$stub);
}


sub get_new_TranslationID {
    my ($self,$stub) = @_;

    $stub  = "ENSP" unless defined($stub);

    return $self->get_NewId("translation",$stub);

}

