# EnsEMBL Gene reading writing adaptor for mySQL
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

Bio::EnsEMBL::DBSQL::GeneAdaptor - MySQL Database queries to generate and store gens.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Elia Stupka  : elia@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut


package Bio::EnsEMBL::DBSQL::GeneAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;


use vars '@ISA';


@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


=head2 list_geneIds

 Title   : list_geneIds
 Usage   : $geneAdaptor->list_geneIds
 Function: Gets an array of internal ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_geneIds {
   my ($self) = @_;

   my @out;
   my $sth = $self->prepare("SELECT gene_id FROM gene");
   $sth->execute;

   while (my ($id) = $sth->fetchrow) {
       push(@out, $id);
   }

   return @out;
}

=head2 list_stable_geneIds

 Title   : list_stable_geneIds
 Usage   : $geneAdaptor->list_stable_geneIds
 Function: Gets an array of stable ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_stable_geneIds {
   my ($self) = @_;

   my @out;
   my $sth = $self->prepare("SELECT stable_id FROM gene_stable_id");
   $sth->execute;

   while (my ($id) = $sth->fetchrow) {
       push(@out, $id);
   }

   return @out;
}


=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $geneobj->fetch_by_dbID( $geneid)
 Function: gets one gene out of the db
 Example : $obj->get($dbID)
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag

=cut


sub fetch_by_dbID {
  my ( $self, $geneId ) = @_;
  my $exonAdaptor = $self->db->get_ExonAdaptor();
  my @exons = $exonAdaptor->fetch_by_geneId( $geneId );

  my %exonIds;

  #
  # We assumme that if we get no exons this is a throwing situation
  #
  if( scalar( @exons ) == 0 ) {
    $self->throw("No exons for gene $geneId, assumming no gene");
  }


  foreach my $exon ( @exons ) {
    $exonIds{$exon->dbID} = $exon;
  }

  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  
  # fetching all exons by gene
  # fetching all transcripts
  # adding the exons
  # adding the transcripts
  my %transcriptExons;
  my %transcripts;

  my $query = qq{
    SELECT tscript.gene_id
      , tscript.transcript_id
      , e_t.exon_id, e_t.rank
      , gene.analysis_id
      , gene.type
      , tscript.translation_id
    FROM gene
      , transcript tscript
      , exon_transcript e_t
    WHERE gene.gene_id = tscript.gene_id
      AND tscript.transcript_id = e_t.transcript_id
      AND gene.gene_id = $geneId

    ORDER BY tscript.gene_id
      , tscript.transcript_id
      , e_t.rank
    };

  my $sth = $self->prepare( $query );
  $sth->execute();

  my $first = 1;
  while( my @arr = $sth->fetchrow_array() ) {
    # building a gene
    if( $first ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->adaptor($self);
      $gene->dbID( $geneId );
      my $ana = $self->db->get_AnalysisAdaptor->fetch_by_dbID($arr[4]);
      $gene->analysis($ana);
      $gene->type($arr[5]);
      $first = 0;
    }

    # store an array of exon ids for each transcript
    if( !exists $transcriptExons{$arr[1]} ) {
      $transcriptExons{$arr[1]} = [];
    }

    push( @{$transcriptExons{$arr[1]}}, $arr[2] );
    $transcripts{$arr[1]} = $arr[6];
  }

  if( $first ) {
    return undef;
  }
  
  foreach my $transcriptId ( keys %transcripts ) {
    # should be fetch_by_geneId ..
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->dbID( $transcriptId );
    $transcript->adaptor( $self->db->get_TranscriptAdaptor() );
    $transcript->_translation_id($transcripts{$transcriptId} );

    foreach my $exonId ( @{$transcriptExons{$transcriptId}} ) {
      $transcript->add_Exon( $exonIds{$exonId} );
    }
    $gene->add_Transcript( $transcript );
  }
  
  return $gene;
}


=head2 fetch_by_stable_id

 Title   : fetch_by_stable_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_stable_id{
   my ($self,$id) = @_;

   my $sth = $self->prepare("select gene_id from gene_stable_id where stable_id = '$id'");
   $sth->execute;

   my ($dbID) = $sth->fetchrow_array();

   if( !defined $dbID ) {
       $self->throw("No stable id with $id, cannot fetch");
   }

   return $self->fetch_by_dbID($dbID);
}

=head2 fetch_by_contig_list

 Title   : fetch_by_contig_list
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_contig_list{
   my ($self,@list) = @_;

   my $str;

   foreach my $id ( @list ) {
       $str .= "'$id',";
   }
   $str =~ s/\,$//g;
   $str = "($str)";

   # 
   # this is non-optimised, because we are going to make multiple
   # trips to the database. should fix here
   #

   my $sth = $self->prepare("select distinct(t.gene_id) from transcript t,exon_transcript et,exon e,contig c where c.name in $str and c.contig_id = e.contig_id and et.exon_id = e.exon_id and et.transcript_id = t.transcript_id");
   $sth->execute;

   my @out;

   while( my ($gid) = $sth->fetchrow_array() ) {
       push(@out,$self->fetch_by_dbID($gid));
   }

   return @out;
	    
}


=head2 fetch_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              The slice we want genes on
  Function  : retrieve all the genes on this slice. Uses contig list from
              Assembly Mapper and Gene->transform. Uses fetch_by_dbID to find
              Genes initially and then makes transformed versions.
  Returntype: list of Bio::EnsEMBL::Gene
  Exceptions: none
  Caller    : Bio::EnsEMBL::Slice

=cut

sub fetch_by_Slice {
  my ( $self, $slice ) = @_;
  my @out;

  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type
    ( $slice->assembly_type() );
  
  $mapper->register_region( $slice->chr_name(),
			    $slice->chr_start(),
			    $slice->chr_end());
  
  my @cids = $mapper->list_contig_ids( $slice->chr_name(),
				       $slice->chr_start(),
				       $slice->chr_end());

  # no genes found so return
  if ( scalar (@cids) == 0 ) {
    return undef;
  }

  my $str = "(".join( ",",@cids ).")";

  my $sth = $self->prepare("
     SELECT distinct(t.gene_id) 
     FROM   transcript t,exon_transcript et,exon e 
     WHERE  e.contig_id in $str 
     AND    et.exon_id = e.exon_id 
     AND    et.transcript_id = t.transcript_id");

  $sth->execute;
  
  while( my ($geneid) = $sth->fetchrow ) {
    my $gene = $self->fetch_by_dbID( $geneid );
    my $newgene = $gene->transform( $slice );    
    push( @out, $newgene );
  }

  return @out;
}



=head2 fetch_by_Transcript_id

 Title   : fetch_by_Transcript_id
 Usage   : $gene_obj->get_Gene_by_Transcript_id($transid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_Transcript_id {
    my $self = shift;
    my $transid = shift;

    # this is a cheap SQL call
    #my $sth = $self->prepare("select gene_id from transcript where transcript_id = '$transid'");
    my $sth = $self->prepare("	SELECT	tr.gene_id 
				FROM	transcript as tr 
				WHERE	tr.transcript_id = $transid");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    my $gene = $self->fetch_by_dbID( $geneid );
    return $gene;
}



=head2 fetch_by_Peptide_id 

 Title   : fetch_by_Peptide_id
 Usage   : $geneAdaptor->fetch_by_Peptide_id($peptideid)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : peptide id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_Peptide_id {
    my $self = shift;
    my $peptideid = shift;

    # this is a cheap SQL call
    my $sth = $self->prepare("	SELECT	tr.gene_id 
				FROM	transcript as tr, 
					translation_stable_id as trs 
			    	WHERE	trs.stable_id = '$peptideid' 
				AND	trs.translation_id = tr.translation_id
			    ");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    return $self->fetch_by_dbID($geneid);
}


=head2 fetch_by_maximum_DBLink

 Title   : fetch_by_maximum_DBLink
 Usage   : $geneAdaptor->fetch_by_maximum_DBLink($ext_id)
 Function: gets one gene out of the db with 
 Returns : gene object (with transcripts, exons)
 Args    : 


=cut

sub fetch_by_maximum_DBLink {
    my $self = shift;
    my $external_id = shift;

    my @genes=$self->fetch_by_DBEntry( $external_id );

    my $biggest;
    my $max=0;
    my $size=scalar(@genes);
    if ($size > 0) {
	foreach my $gene (@genes) {
	    my $size = (scalar($gene->each_unique_Exon));
	    if ($size > $max) {
		$biggest = $gene;
		$max=$size;
	    }
	}
	return $biggest;
    }
    return;
}


=head2 get_description

 Title   : get_description
 Usage   : $geneAdptor->get_description($dbID)
 Function: gets gene description line 
 Returns : a string
 Args    : 


=cut

sub get_description {
  my ($self,$dbID) = @_;

  if( !defined $dbID ) {
      $self->throw("must call with dbID");
  }

  my $sth = $self->prepare("select description from gene_description where gene_id = $dbID");
  $sth->execute;
  my @array = $sth->fetchrow_array();
  return $array[0];
}

=head2 get_stable_entry_info

 Title   : get_stable_entry_info
 Usage   : $geneAdptor->get_stable_entry_info($gene)
 Function: gets stable info for gene and places it into the hash
 Returns : 
 Args    : 


=cut

sub get_stable_entry_info {
  my ($self,$gene) = @_;

  if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     $self->throw("Needs a gene object, not a $gene");
  }

  my $sth = $self->prepare("select stable_id,UNIX_TIMESTAMP(created),UNIX_TIMESTAMP(modified),version from gene_stable_id where gene_id = ".$gene->dbID);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $gene->{'_stable_id'} = $array[0];
  $gene->{'_created'}   = $array[1];
  $gene->{'_modified'}  = $array[2];
  $gene->{'_version'}   = $array[3];
  

  return 1;
}

=head2 fetch_by_DBEntry

 Title   : fetch_by_DBLink
 Usage   : $geneAdptor->fetch_by_DBLink($ext_id)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_DBEntry {
  my $self = shift;
  my $external_id = shift;
  my @genes = ();

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();


  my @ids = $entryAdaptor->geneids_by_extids($external_id);
  foreach my $gene_id ( @ids ) {
    my $gene = $self->fetch_by_dbID( $gene_id );
    if( $gene ) {
      push( @genes, $gene );
    }
  }
  return @genes;
}




=head2 store

 Title   : store
 Usage   : $geneAdaptor->store($gene)
 Function: writes a particular gene into the database. Assumes that everything 
           has dbIDs ....
 Example :
 Returns : nothing
 Args    : $gene object


=cut

sub store {
   my ($self,$gene) = @_;
   my %done;
   #print STDERR "storing gene\n";
   my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
   #print STDERR "have transcript adaptor\n";
   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Must store a gene object, not a $gene");
   }

   if( !defined $gene->analysis ) {
       $self->throw("Genes must have an analysis object!");
   }
   #print STDERR "storing analysis_object\n";
   my $analysisId = $self->db->get_AnalysisAdaptor()->store( $gene->analysis );


   if ( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
   }
 
   my $trans_count = scalar($gene->each_Transcript);
   print $trans_count."\n";
   print $gene->type()."\n";
   my $type = $gene->type;
   #print STDERR "inserting into genetable\n";
   my $sth2 = $self->prepare("insert into gene(type, analysis_id, transcript_count) values('$type', $analysisId, $trans_count)" );
   $sth2->execute();
   
   $gene->adaptor( $self );
   $gene->dbID( $sth2->{'mysql_insertid'} );
   #print STDERR "have gene dbID\n";
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();
   #print STDERR "have dbEntryAdaptor\n";
  
   foreach my $dbl ( $gene->each_DBLink ) {
     #print STDERR $dbl."\n";
     $dbEntryAdaptor->store( $dbl, $gene->dbID, "Gene" );
   }
   #print "have stored all dbLinks\n";
   # write exons at this level to avoid duplicates
   my $exonAdaptor = $self->db->get_ExonAdaptor();
   my @ex = $gene->get_all_Exons;

   foreach my $exon($gene->get_all_Exons){
     $exonAdaptor->store( $exon );
   }


   # write exons transcripts and exon_transcript table
   foreach my $trans ( $gene->each_Transcript() ) {
     $transcriptAdaptor->store($trans,$gene);
   }
   
   return $gene->dbID;
}


sub remove {
  my $self = shift;
  my $gene = shift;

  if( ! defined $gene->dbID() ) {
    return;
  }

  my $sth= $self->prepare( "delete from gene where gene_id = ? " );
  $sth->execute( $gene->dbID );
  $sth= $self->prepare( "delete from gene_stable_id where gene_id = ? " );
  $sth->execute( $gene->dbID );
  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  # my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();
  foreach my $trans ( $gene->each_Transcript() ) {
    $transcriptAdaptor->remove($trans,$gene);
  }

  $gene->{'_dbID'} = undef;
}


=head2 get_Interpro_by_geneid

 Title   : get_Interpro_by_geneid
 Usage   : @interproid = $geneAdaptor->get_Interpro_by_geneid($gene->id);
 Function: gets interpro accession numbers by geneid. A hack really -
           we should have a much more structured system than this
 Example :
 Returns : 
 Args    :


=cut

sub get_Interpro_by_geneid {
   my ($self,$gene) = @_;
    my $sql="
	SELECT	i.interpro_ac, 
		x.description 
        FROM	transcript t, 
		protein_feature pf, 
		interpro i, 
                xref x,
		gene_stable_id gsi
	WHERE	gsi.stable_id = '$gene' 
	    AND	t.gene_id = gsi.gene_id
	    AND	t.translation_id = pf.translation 
	    AND	i.id = pf.hit_id 
	    AND	i.interpro_ac = x.dbprimary_acc";
   
   my $sth = $self->prepare($sql);
   $sth->execute;

   my @out;
   my %h;
   while( (my $arr = $sth->fetchrow_arrayref()) ) {
       if( $h{$arr->[0]} ) { next; }
       $h{$arr->[0]}=1;
       my $string = $arr->[0] .":".$arr->[1];
       
       push(@out,$string);
   }


   return @out;
}





sub create_tables {
# read sql from geneAdaptor.sql
}



1;
__END__

