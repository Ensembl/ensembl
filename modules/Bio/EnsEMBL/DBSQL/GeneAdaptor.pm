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

$gene_adaptor = $db_adaptor->get_GeneAdaptor();

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


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 list_geneIds

  Arg [1]    : none
  Example    : @gene_ids = $gene_adaptor->list_geneIds();
  Description: Gets an array of internal ids for all genes in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

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

  Arg [1]    : list_stable_gene_ids
  Example    : @stable_ids = $gene_adaptor->list_stable_gene_ids();
  Description: Returns a list stable ids for all genes in the current db
  Returntype : list of strings
  Exceptions : none
  Caller     : ?

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

  Arg [1]    : int $geneId 
               the unique internal database id of the Gene to be retrieved
  Example    : $gene = $gene_adaptor->fetch_by_dbID
  Description: Retrieves a gene object from the database using its unique
               internal identifier.
  Returntype : Bio::EnsEMBL::Gene in contig coordinates
  Exceptions : thrown if no exons exist for the gene with dbID $geneId
  Caller     : general

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
      $gene->source($self->db->source());
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

  Arg [1]    : string $id 
               The stable id of the gene to retrieve
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944');
  Description: Retrieves a gene object from the database via its stable id
  Returntype : Bio::EnsEMBL::Gene in contig coordinates
  Exceptions : thrown if no gene of stable_id $id exists in the database
  Caller     : general

=cut

sub fetch_by_stable_id{
   my ($self,$id) = @_;

   my $sth = $self->prepare("SELECT gene_id 
                             FROM gene_stable_id 
                             WHERE stable_id = '$id'");
   $sth->execute;

   my ($dbID) = $sth->fetchrow_array();

   if( !defined $dbID ) {
       $self->throw("No stable id with $id, cannot fetch");
   }

   return $self->fetch_by_dbID($dbID);
}


=head2 fetch_by_contig_list

  Arg [1]    : list of ints @list
               the contigs to retrieve genes from
  Example    : @genes = $gene_adaptor->fetch_by_contig_list(1, 2, 3, 4);
  Description: Retrieves all genes which are present on list of contigs
               denoted by their unique database ids
  Returntype : list of Bio::EnsEMBL::Genes in contig coordinates
  Exceptions : none
  Caller     : general

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

   my $sth = $self->prepare("SELECT distinct(t.gene_id) 
                             FROM transcript t,exon_transcript et,
                                  exon e,contig c 
                             WHERE c.name IN $str 
                             AND c.contig_id = e.contig_id 
                             AND et.exon_id = e.exon_id 
                             AND et.transcript_id = t.transcript_id");
   $sth->execute;

   my @out;

   while( my ($gid) = $sth->fetchrow_array() ) {
       push(@out,$self->fetch_by_dbID($gid));
   }

   return @out;
}


=head2 fetch_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice to fetch genes from
  Example    : $genes = $gene_adaptor->fetch_by_slice($slice);
  Description: Retrieves all genes which are present on a slice
  Returntype : list of Bio::EnsEMBL::Genes in slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_by_Slice {
  my ( $self, $slice, $type ) = @_;
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



=head2 fetch_by_transcript_id

  Arg [1]    : int $transid 
               unique database identifier for the transcript whose gene should
               be retrieved.
  Example    : $gene = $gene_adaptor->fetch_by_transcript_id($transcript);
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : ?

=cut

sub fetch_by_transcript_id {
    my $self = shift;
    my $transid = shift;

    # this is a cheap SQL call
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

  Arg [1]    : string $peptideid
               the stable id of a translation of the gene that should
               be obtained
  Example    : $gene = $gene_adaptor->fetch_by_Peptide_id('ENSP00000278194');
  Description: retrieves a gene via the stable id of one of its translations
               this method should be renamed, but for now it stays
  Returntype : Bio::EnsEMBL::Gene in contig coordinates
  Exceptions : none
  Caller     : geneview

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

  Arg [1]    : int $external_id
               the unique identifier of the DBEntry of the gene to retrieve
  Example    : my $gene = $gene_adaptor->fetch_by_maximum_DBLink($ext_id);
  Description: retrieves the gene with the  most most exons with the external 
               database link identified by $external_id
  Returntype : Bio::EnsEMBL::Gene in contig coordinates
  Exceptions : none
  Caller     : ?

=cut


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
	    my $size = (scalar($gene->get_all_Exons));
	    if ($size > $max) {
		$biggest = $gene;
		$max=$size;
	    }
	}
	return $biggest;
    }
    return;
}


=head2 get_stable_entry_info

  Arg [1]    : Bio::EnsEMBL::Gene $gene
  Example    : $gene_adaptor->get_stable_entry_info($gene);
  Description: gets stable info for a gene. this is not usually done at
               creation time for speed purposes, and can be lazy-loaded later
               if it is needed..
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::Gene

=cut

sub get_stable_entry_info {
  my ($self,$gene) = @_;

  if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     $self->throw("Needs a gene object, not a $gene");
  }

  my $sth = $self->prepare("SELECT stable_id, UNIX_TIMESTAMP(created),
                                   UNIX_TIMESTAMP(modified), version 
                            FROM gene_stable_id 
                            WHERE gene_id = ".$gene->dbID);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $gene->{'_stable_id'} = $array[0];
  $gene->{'_created'}   = $array[1];
  $gene->{'_modified'}  = $array[2];
  $gene->{'_version'}   = $array[3];
  

  return 1;
}


=head2 fetch_by_DBEntry

  Arg [1]    : in $external_id
               the external identifier for the gene to be obtained
  Example    : @genes = $gene_adaptor->fetch_by_DBEntry($ext_id)
  Description: retrieves a list of genes with an external database 
               idenitifier $external_id
  Returntype : list of Bio::EnsEMBL::DBSQL::Gene in contig coordinates
  Exceptions : none
  Caller     : ?

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

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $gene_adaptor->store($gene);
  Description: Stores a gene in the database
  Returntype : the database identifier of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene or if 
               $gene does not have an analysis object
  Caller     : general

=cut

sub store {
   my ($self,$gene) = @_;
   my %done;

   my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();

   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Must store a gene object, not a $gene");
   }
   if( !defined $gene->analysis ) {
       $self->throw("Genes must have an analysis object!");
   }

   my $aA = $self->db->get_AnalysisAdaptor();
   my $analysisId = $aA->exists( $gene->analysis() );

   if( defined $analysisId ) {
     $gene->analysis()->dbID( $analysisId );
   } else {
     $analysisId = $self->db->get_AnalysisAdaptor()->store( $gene->analysis );
   }

   if ( !defined $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("$gene is not a EnsEMBL gene - not writing!");
   }
 
   my $trans_count = scalar($gene->get_all_Transcripts);
   my $type = $gene->type;
   my $sth2 = $self->prepare("INSERT INTO gene(type, analysis_id, 
                                               transcript_count) 
                              VALUES('$type', $analysisId, $trans_count)" );
   $sth2->execute();
   
   $gene->adaptor( $self );
   $gene->dbID( $sth2->{'mysql_insertid'} );

   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbl ( $gene->each_DBLink ) {
     $dbEntryAdaptor->store( $dbl, $gene->dbID, "Gene" );
   }

   # write exons at this level to avoid duplicates
   my $exonAdaptor = $self->db->get_ExonAdaptor();
   my @ex = $gene->get_all_Exons;

   foreach my $exon($gene->get_all_Exons){
     $exonAdaptor->store( $exon );
   }

   # write exons transcripts and exon_transcript table
   foreach my $trans ( $gene->get_all_Transcripts() ) {
     $transcriptAdaptor->store($trans,$gene);
   }
   
   return $gene->dbID;
}


=head2 remove

  Arg [1]    : Bio::EnsEMBL::Gene $gene 
               the gene to remove from the database 
  Example    : $gene_adaptor->remove($gene);
  Description: Removes a gene from the database
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

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
  foreach my $trans ( $gene->get_all_Transcripts() ) {
    $transcriptAdaptor->remove($trans,$gene);
  }

  $gene->{'_dbID'} = undef;
}



=head2 get_Interpro_by_geneid

  Arg [1]    : string $gene
               the stable if of the gene to obtain
  Example    : @i = $gene_adaptor->get_Interpro_by_geneid($gene->stable_id()); 
  Description: gets interpro accession numbers by gene stable id.
               A hack really - we should have a much more structured 
               system than this
  Returntype : list of strings 
  Exceptions : none 
  Caller     : domainview?

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
	    AND	t.translation_id = pf.translation_id 
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



=head2 get_description

  Arg [1]    : int $dbID
               the database identifier of the gene whose description should be
               obtained
  Example    : $description = $gene_adaptor->get_description(1234);
  Description: Retrieves a description for a gene which is created during
               the gene build and stored in the the gene_description table.
               In the future this method should check the family database
               for a description if one does not exist in the core, however,
               this has not been implemented yet.  
  Returntype : string
  Exceptions : thrown if $dbId arg is not defined
  Caller     : geneview

=cut

sub get_description {
  my ($self, $dbID) = @_;

  if( !defined $dbID ) {
      $self->throw("must call with dbID");
  }

  my $sth = $self->prepare("select description from gene_description where gene_id = $dbID");
  $sth->execute;
  my @array = $sth->fetchrow_array();
  return $array[0];
}


1;
__END__

