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

use strict;

use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);




=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal

=cut

sub _tables {
  my $self = shift;
  
  return (['gene', 'g'],[ 'gene_description', 'gd'], ['gene_stable_id', 'gsi'] );
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

sub _columns {
  my $self = shift;

  return qw( g.seq_region_id g.seq_region_start g.seq_region_end g.seq_region_strand
             g.analysis_id g.type g.display_xref_id gd.description gsi.stable_id 
             gsi.version  );
}




=head2 list_dbIDs

  Arg [1]    : none
  Example    : @gene_ids = @{$gene_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all genes in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("gene");
}

=head2 list_stable_ids

  Arg [1]    : none
  Example    : @stable_gene_ids = @{$gene_adaptor->list_stable_dbIDs()};
  Description: Gets an array of stable ids for all genes in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("gene_stable_id", "stable_id");
}



=head2 fetch_by_dbID

  Arg [1]    : int $geneId 
               the unique internal database id of the Gene to be retrieved
  Arg [2]    : int $chromosomal_coordinates (optional)
               if defined, try to return chromosomal coordinates.
  Example    : $gene = $gene_adaptor->fetch_by_dbID
  Description: Retrieves a gene object from the database using its unique
               internal identifier.
  Returntype : Bio::EnsEMBL::Gene in contig coordinates
  Exceptions : thrown if no exons exist for the gene with dbID $geneId
  Caller     : general

=cut

sub fetch_by_dbID {
  my ( $self, $geneId, $chr_coordinates ) = @_;

  my $exonAdaptor = $self->db->get_ExonAdaptor();
  my @exons = @{$exonAdaptor->fetch_all_by_gene_id( $geneId )};
  my %exonIds;

  #
  # We assumme that if we get no exons this is a throwing situation
  #
  if( scalar( @exons ) == 0 ) {
    $self->throw("No exons for gene $geneId, assumming no gene");
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
  my $ana;
  my $first = 1;
  my $gene;

  while( my @arr = $sth->fetchrow_array() ) {
    # building a gene
    if( $first ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->adaptor($self);
      $gene->dbID( $geneId );
      $ana = $self->db->get_AnalysisAdaptor->fetch_by_dbID($arr[4]);
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
  
  #discard duplicate exons, add analysis object to exons
  foreach my $exon ( @exons ) {
    $exon->analysis($ana);
    $exonIds{$exon->dbID} = $exon;
  }
  
  foreach my $transcriptId ( keys %transcripts ) {
    # should be fetch_by_geneId ..
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->dbID( $transcriptId );
    $transcript->adaptor( $self->db->get_TranscriptAdaptor() );
    
    # Test for truth because translation_id will be 0 if not set
    # because column is NOT NULL in schema.
    if (my $translation_id = $transcripts{$transcriptId}) {
        $transcript->_translation_id($translation_id);
    }

    foreach my $exonId ( @{$transcriptExons{$transcriptId}} ) {
      
      $transcript->add_Exon( $exonIds{$exonId} );
    }
    # store also a type for the transcript
    $transcript->type($gene->type);
    $gene->add_Transcript( $transcript );
  }
  
  # if chromosomal coordinates are needed, transform with empty slice

  if( defined $chr_coordinates ) {
    my $sa = $self->db->get_SliceAdaptor();
    my $empty_slice = Bio::EnsEMBL::Slice->new
      ( 
       -empty => 1,
       -adaptor => $sa 
      );
    $gene->transform( $empty_slice );
  }

  return $gene;
}


=head2 fetch_by_stable_id

  Arg  1     : string $id 
               The stable id of the gene to retrieve
  Arg [2]    : string $coordinate_system_name
  Arg [3]    : string $coordinate_system_version

  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944', 'chromosome');
  Description: Retrieves a gene object from the database via its stable id
  Returntype : Bio::EnsEMBL::Gene in given coordinate system
  Exceptions : if we cant get the gene in given coord system
  Caller     : general

=cut

sub fetch_by_stable_id{
   my ($self,$id, $cs_name, $cs_version) = @_;

   my $constraint = "gsi.stable_id = $id";

   # should be only one :-)
   my $genes = $self->SUPER::generic_fetch( $constraint );

   if( ! @$genes ) { return undef }

   my @new_genes = map { $_->transform( $cs_name, $cs_version ) } @$genes;


   return $new_genes->[0];
 }


=head2 fetch_by_exon_stable_id

  Arg [1]    : string $id
               The stable id of the gene to retrieve
  Arg [2]    : (optional) boolean $chr_coords
               flag indicating genes should be returned in chromosomal
               coords instead of contig coords.
  Example    : $gene = $gene_adaptor->fetch_by_exon_stable_id('ENSE00000148944');
  Description: Retrieves a gene object from the database via an exon stable id
  Returntype : Bio::EnsEMBL::Gene in contig coordinates by default or in
               chromosomal coords if the $chr_coords arg is set to 1.
  Exceptions : thrown if no gene of stable_id $id exists in the database
  Caller     : general

=cut

sub fetch_by_exon_stable_id{
   my ($self,$id, $chr_coords) = @_;

   my $sth = $self->prepare(
     "SELECT t.gene_id
        FROM transcript as t,
             exon_transcript as et,
             exon_stable_id as esi
       WHERE t.transcript_id = et.transcript_id 
         AND et.exon_id = esi.exon_id 
         AND esi.stable_id = '$id'");
   $sth->execute;

   my ($dbID) = $sth->fetchrow_array();

   if( !defined $dbID ) {
       $self->throw("No stable id with $id, cannot fetch");
   }

   my $gene = $self->fetch_by_dbID($dbID);

   if($chr_coords) {
     #transform gene to chromosomal coords
     my $slice_adaptor = $self->db->get_SliceAdaptor;
     my $slice = new Bio::EnsEMBL::Slice(-empty => 1,
                                         -adaptor => $slice_adaptor);
     $gene->transform($slice);
   }
  
   return $gene;
 }




=head2 fetch_all_by_domain

  Arg [1]    : string $domain
               the domain to fetch genes from
  Arg [2]    : (optional) boolean $empty_flag
               true if lightweight genes are desired (for speed purposes)
  Example    : my @genes = $gene_adaptor->fetch_all_by_domain($domain);
  Description: retrieves a listref of genes whose translation contain interpro
               domain $domain.
  Returntype : list of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : domainview

=cut

sub fetch_all_by_domain {
  my ($self, $domain, $empty_flag) = @_;

  unless($domain) {
    $self->throw("domain argument is required");
  }

  my $sth = $self->prepare("SELECT tr.gene_id
                            FROM interpro i,
                                 protein_feature pf,
                                 transcript tr
                            WHERE i.interpro_ac = ?
                            AND   i.id = pf.hit_id
                            AND   pf.translation_id = tr.translation_id
                            GROUP BY tr.gene_id");
 
  $sth->execute($domain);
  my @gene_ids = ();
  my $gene_id;

  $sth->bind_columns(\$gene_id);
  while($sth->fetch()) {
    push @gene_ids, $gene_id;
  }

  #may want to use proxy...
  return $self->db->get_GeneAdaptor->fetch_all_by_gene_id_list(\@gene_ids, 
							   $empty_flag);
}


  
=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice to fetch genes from
  Arg [2]
  Example    : $genes = $gene_adaptor->fetch_all_by_Slice($slice);
  Description: Retrieves all genes which are present on a slice
  Returntype : listref of Bio::EnsEMBL::Genes in slice coordinates
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice {
  my ( $self, $slice, $logic_name ) = @_;

  my @out = ();

  $logic_name ||= '';

  my $key = uc($slice->name .":" . $logic_name);

  #check the cache which uses the slice name as it key
  if($self->{'_slice_gene_cache'}{$key}) {
    return $self->{'_slice_gene_cache'}{$key};
  }

  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type
    ( $slice->assembly_type() );
  
  my @cids = $mapper->list_contig_ids( $slice->chr_name(),
				       $slice->chr_start(),
				       $slice->chr_end());
  
  # no genes found so return
  if ( scalar (@cids) == 0 ) {
    return [];
  }

  my $str = "(".join( ",",@cids ).")";

  my $where = "WHERE e.contig_id in $str 
               AND   et.exon_id = e.exon_id 
               AND   et.transcript_id = t.transcript_id
               AND   g.gene_id = t.gene_id";

  if($logic_name) {
    #determine analysis id via logic_name
    my $analysis = 
      $self->db->get_AnalysisAdaptor()->fetch_by_logic_name($logic_name);
    unless(defined($analysis) && $analysis->dbID()) {
      $self->warn("No analysis for logic name $logic_name exists");
      return [];
    }
    
    my $analysis_id = $analysis->dbID;
    $where .= " AND g.analysis_id = $analysis_id";
  }
    
  my $sql = "
    SELECT distinct(t.gene_id)
    FROM   transcript t,exon_transcript et,exon e, gene g 
    $where";
    
  my $sth = $self->db->prepare($sql);
  $sth->execute;
  
  while( my ($geneid) = $sth->fetchrow ) {
    my $gene = $self->fetch_by_dbID( $geneid );
    my $newgene = $gene->transform( $slice );

    if( $newgene->start() <= $slice->length() &&
	$newgene->end() >= 1 ) {
      # only take the gene if its really overlapping the Slice
      push( @out, $newgene );
    }
 }

  #place the results in an LRU cache
  $self->{'_slice_gene_cache'}{$key} = \@out;

  return \@out;
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

=head2 fetch_by_transcript_stable_id

  Arg [1]    : string $transid 
               unique database identifier for the transcript whose gene should
               be retrieved.
  Example    : none
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : ?

=cut

sub fetch_by_transcript_stable_id {
    my $self = shift;
    my $transid = shift;

    # this is a cheap SQL call
     my $sth = $self->prepare("	SELECT	tr.gene_id 
				FROM	transcript as tr, transcript_stable_id tcl 
				WHERE	tcl.stable_id = \"$transid\"
                                  AND   tr.transcript_id = tcl.transcript_id");
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

 Title   : fetch_by_maximum_DBLink
 Usage   : $geneAdaptor->fetch_by_maximum_DBLink($ext_id)
 Function: gets one gene out of the db with 
 Returns : gene object (with transcripts, exons)
 Args    : 


=cut

sub fetch_by_maximum_DBLink {
  my $self = shift;
  my $external_id = shift;

    my $genes=$self->fetch_all_by_DBEntry( $external_id );

    my $biggest;
    my $max=0;
    my $size=scalar(@$genes);
    if ($size > 0) {
	foreach my $gene (@$genes) {
	    my $size = scalar(@{$gene->get_all_Exons});
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


=head2 fetch_all_by_DBEntry

  Arg [1]    : in $external_id
               the external identifier for the gene to be obtained
  Example    : @genes = @{$gene_adaptor->fetch_all_by_DBEntry($ext_id)}
  Description: retrieves a list of genes with an external database 
               idenitifier $external_id
  Returntype : listref of Bio::EnsEMBL::DBSQL::Gene in contig coordinates
  Exceptions : none
  Caller     : ?

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;
  my $external_id = shift;
  my @genes = ();

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();


  my @ids = $entryAdaptor->list_gene_ids_by_extids($external_id);
  foreach my $gene_id ( @ids ) {
    my $gene = $self->fetch_by_dbID( $gene_id );
    if( $gene ) {
      push( @genes, $gene );
    }
  }
  return \@genes;
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
 
   my $trans_count = scalar(@{$gene->get_all_Transcripts});
   my $type = "";

   if (defined($gene->type)) {
       $type = $gene->type;
   }

   # assuming that the store is used during the Genebuil process, set
   # the display_xref_id to 0.  This ought to get re-set during the protein
   # pipeline run.  This probably update to the gene table has yet to be
   # implemented.
   my $xref_id = 0;

   my $sth2 = $self->prepare("INSERT INTO gene(type, analysis_id, 
                                               transcript_count, 
                                               display_xref_id) 
                              VALUES(?,?,?,?)");
   $sth2->execute("$type", $analysisId, $trans_count, $xref_id);

   my $gene_dbID = $sth2->{'mysql_insertid'};

   if (defined($gene->stable_id)) {
     if (!defined($gene->created) || 
         !defined($gene->modified) ||
         !defined($gene->version)) {
       $self->throw("Trying to store incomplete stable id information for gene");
     }

     my $statement = "INSERT INTO gene_stable_id(gene_id," .
                                   "version, stable_id, created, modified)".
                      " VALUES(" . $gene_dbID . "," .
                               $gene->version . "," .
                               "'" . $gene->stable_id . "'," .
                               "FROM_UNIXTIME(".$gene->created."),".
                               "FROM_UNIXTIME(".$gene->modified."))";
     my $sth = $self->prepare($statement);
     $sth->execute();
   }


   #
   # store the dbentries associated with this gene
   #
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbe ( @{$gene->get_all_DBEntries} ) {
     $dbEntryAdaptor->store( $dbe, $gene_dbID, "Gene" );
   }

   #
   # Update the genes display xref if it is set
   #
   my $display_xref = $gene->display_xref;
   if($display_xref) {
     if(my $dbe_id = $dbEntryAdaptor->exists($display_xref)) {
       my $dispxref_sth = $self->prepare('UPDATE gene SET display_xref_id = ? 
                                        WHERE gene_id = ?');
       $dispxref_sth->execute($dbe_id, $gene_dbID);
     }
   }
       
   # write exons transcripts and exon_transcript table
   my $trans = $gene->get_all_Transcripts;

   #force lazy loading of translations before new exon dbIDs are set
   map {$_->translation} @$trans;

   foreach my $t ( @$trans ) {
     $transcriptAdaptor->store($t,$gene_dbID );
   }

   $gene->adaptor( $self );
   $gene->dbID( $gene_dbID );
   
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
  foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
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
  Returntype : listref of strings 
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


   return \@out;
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

  my $sth = $self->prepare("SELECT description 
                            FROM   gene_description 
                            WHERE  gene_id = ?");
  $sth->execute($dbID);
  my @array = $sth->fetchrow_array();
  return $array[0];
}



=head2 deleteObj

  Arg [1]    : none
  Example    : none 
  Description: Responsible for cleaning up this objects references to other
               objects so that proper garbage collection can occur. 
  Returntype : none
  Exceptions : none
  Caller     : DBConnection::DeleteObj

=cut

sub deleteObj {
  my $self = shift;

  #print STDERR "\t\tGeneAdaptor::deleteObj\n";
  
  #call superclass destructor
  $self->SUPER::deleteObj();

  #flush the cache
  %{$self->{'_slice_gene_cache'}} = ();
}


							

=head2 get_display_xref

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               retrieve display_xref for this gene
  Example    : none
  Description: Retrieves the display_xref for a gene.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : thrown if $dbId arg is not defined
  Caller     : general

=cut

sub get_display_xref {
  my ($self, $gene ) = @_;

  if( !defined $gene ) {
      $self->throw("Must call with a Gene object");
  }

  my $sth = $self->prepare("SELECT e.db_name,
                                   x.display_label,
                                   x.xref_id
                            FROM   gene g, 
                                   xref x, 
                                   external_db e
                            WHERE  g.gene_id = ?
                              AND  g.display_xref_id = x.xref_id
                              AND  x.external_db_id = e.external_db_id
                           ");
  $sth->execute( $gene->dbID );


  my ($db_name, $display_label, $xref_id ) = $sth->fetchrow_array();
  if( !defined $xref_id ) {
    return undef;
  }
  my $db_entry = Bio::EnsEMBL::DBEntry->new
    (
     -dbid => $xref_id,
     -adaptor => $self->db->get_DBEntryAdaptor(),
     -dbname => $db_name,
     -display_id => $display_label
    );

  return $db_entry;
}



=head2 update

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $gene_adaptor->update($gene);
  Description: Updates  a gene in the database
  Returntype : None
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene
               warn if trying to update the number of attached transcripts.  This
               is a far more complex process and is not yet implemented.
               warn if the method is called on a gene that does not exist in the
               database.
  Caller     : general

=cut

sub update {
   my ($self,$gene) = @_;
   my $update = 0;

   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     $self->throw("Must update a gene object, not a $gene");
   }

   my $xref_id;
   if($gene->display_xref && $gene->display_xref->dbID) {
     $xref_id = $gene->display_xref->dbID;
   } else {
     $xref_id = 0;
   }

   my $tcount = scalar @{$gene->get_all_Transcripts};

   my $sth = $self->prepare("UPDATE gene
                          SET    type = ?,
                                 analysis_id = ?,
                                 display_xref_id = ?,
                                 transcript_count = ?
                          WHERE  gene_id = ?");

   $sth->execute($gene->type, 
		 $gene->analysis->dbID, 
		 $xref_id, 
		 $tcount, 
		 $gene->dbID);
}


1;
__END__

