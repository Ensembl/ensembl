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

Bio::EnsEMBL::DBSQL::GeneAdaptor - A database aware adaptor responsible for the
retrieval and storage of Gene objects.

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::DBAdaptor;

  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(...);

  $gene_adaptor = $db->get_GeneAdaptor();

  $gene = $gene_adaptor->fetch_by_dbID(1234);

  $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000184129');

  @genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA2')};

  $slice_adaptor = $db->get_SliceAdaptor;
  $slice = $slice_adaptor->fetch_by_region('chromosome', '1', 1, 1000000);
  @genes = @{$gene_adaptor->fetch_all_by_Slice($slice)};


=head1 CONTACT

Contact the EnsEMBL development list for questions or information:
ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::GeneAdaptor;

use strict;


use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );

use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns the names, aliases of the tables to use for queries
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _tables {
  my $self = shift;

  return ([ 'gene', 'g' ],
          [ 'gene_stable_id', 'gsi' ],
          [ 'xref', 'x' ],
          [ 'external_db' , 'exdb' ]);
}


# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method
#               returns a list of columns to use for queries
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _columns {
  my $self = shift;

  return qw( g.gene_id g.seq_region_id g.seq_region_start g.seq_region_end
	     g.seq_region_strand g.analysis_id g.biotype g.display_xref_id
	     g.description g.status g.source 
	     gsi.stable_id gsi.version UNIX_TIMESTAMP(gsi.created_date)
	     UNIX_TIMESTAMP(gsi.modified_date) 
	     x.display_label x.dbprimary_acc x.description x.version 
	     exdb.db_name exdb.status exdb.release exdb.db_display_name );
}


sub _left_join {
  return ( [ 'gene_stable_id', "gsi.gene_id = g.gene_id" ],
	   [ 'xref', "x.xref_id = g.display_xref_id" ],
	   [ 'external_db', "exdb.external_db_id = x.external_db_id" ] );
}



=head2 list_dbIDs

  Arg [1]    : none
  Example    : @gene_ids = @{$gene_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all genes in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("gene");
}

=head2 list_stable_ids

  Arg [1]    : none
  Example    : @stable_gene_ids = @{$gene_adaptor->list_stable_dbIDs()};
  Description: Gets an listref of stable ids for all genes in the current db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("gene_stable_id", "stable_id");
}


=head2 fetch_by_display_label

  Arg [1]    : string $label
  Example    : my $gene = $geneAdaptor->fetch_by_display_label( "BRCA2" );
  Description: returns the gene which has the given display label or undef if
               there is none. If there are more than 1, only the first is reported.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_display_label {
  my $self = shift;
  my $label = shift;

  my ( $gene ) = @{$self->generic_fetch( "x.display_label = \"$label\"" )};
  return $gene;
}



=head2 fetch_by_stable_id

  Arg  1     : string $id 
               The stable id of the gene to retrieve
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944');
  Description: Retrieves a gene object from the database via its stable id.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database).  It may
               be converted to a different coordinate system through a call to
               transform() or transfer().  If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene in given coordinate system
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
   my ($self,$id) = @_;

   #because of the way this query is constructed (with a left join to the
   # gene_stable_id table), it is faster to do 2 queries, getting the gene_id
   # in the first query

   my $sth = $self->prepare("SELECT gene_id from gene_stable_id " . 
                            "WHERE  stable_id = ?");
   $sth->execute($id);

   my ($dbID) = $sth->fetchrow_array();
   $sth->finish;

   return undef if(!$dbID);

   return $self->fetch_by_dbID($dbID);
 }


=head2 fetch_by_exon_stable_id

  Arg [1]    : string $id
               The stable id of an exon of the gene to retrieve
  Example    : $gene = $gene_adptr->fetch_by_exon_stable_id('ENSE00000148944');
  Description: Retrieves a gene object from the database via an exon stable id.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database).  It may
               be converted to a different coordinate system through a call to
               transform() or transfer().  If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene (or undef)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_exon_stable_id{
   my ( $self, $id ) = @_;

   my $sth = $self->prepare(
     "SELECT t.gene_id
        FROM transcript as t,
             exon_transcript as et,
             exon_stable_id as esi
       WHERE t.transcript_id = et.transcript_id 
         AND et.exon_id = esi.exon_id 
         AND esi.stable_id = ?");
   $sth->execute("$id");

   my ($dbID) = $sth->fetchrow_array();

   return undef if(!defined($dbID));

   my $gene = $self->fetch_by_dbID( $dbID );

   return $gene;
}




=head2 fetch_all_by_domain

  Arg [1]    : string $domain
               the domain to fetch genes from
  Example    : my @genes = $gene_adaptor->fetch_all_by_domain($domain);
  Description: retrieves a listref of genes whose translation contain interpro
               domain $domain.  The genes are returned in their native coord
               system (i.e. the coord_system they are stored in). If the coord
               system needs to be changed, then tranform
               or transfer should be called on the individual objects returned.
  Returntype : list of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : domainview
  Status     : Stable

=cut

sub fetch_all_by_domain {
  my ($self, $domain) = @_;

  unless($domain) {
    throw("domain argument is required");
  }

  my $sth = $self->prepare("SELECT tr.gene_id " .
                           "FROM interpro i, " .
                               " protein_feature pf, " .
                               " transcript tr, " .
                               " translation tl " .
                           "WHERE i.interpro_ac = ? " .
                           "AND   i.id = pf.hit_id " .
                           "AND   pf.translation_id = tl.translation_id ".
                           "AND   tr.transcript_id = tl.transcript_id " .
                           "GROUP BY tr.gene_id");

  $sth->execute($domain);

  my @array = @{$sth->fetchall_arrayref()};
  $sth->finish();
  
  my @gene_ids = map {$_->[0]} @array;

  return $self->fetch_all_by_dbID_list(\@gene_ids);
}



=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on.
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately rather than
               lazy loaded later.
  Example    : @genes = @{$gene_adaptor->fetch_all_by_Slice()};
  Description: Overrides superclass method to optionally load transcripts
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of genes whose
               transcripts are going to be used.
  Returntype : reference to list of transcripts
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Transcripts
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self  = shift;
  my $slice = shift;
  my $logic_name = shift;
  my $load_exons = shift;

  my $genes = $self->SUPER::fetch_all_by_Slice($slice, $logic_name);

  # if there are 0 or 1 genes still do lazy-loading
  if(!$load_exons || @$genes < 2) {
    return $genes;
  }

  # preload all of the transcripts now, instead of lazy loading later
  # faster than 1 query per transcript

  # first check if transcripts are already preloaded
  # coorectly we should check all of them ..
  return $genes if( exists $genes->[0]->{'_transcript_array'} );

  # get extent of region spanned by transcripts
  my ($min_start, $max_end);
  foreach my $g (@$genes) {
    if(!defined($min_start) || $g->seq_region_start() < $min_start) {
      $min_start = $g->seq_region_start();
    }
    if(!defined($max_end) || $g->seq_region_end() > $max_end) {
      $max_end   = $g->seq_region_end();
    }
  }

  my $ext_slice;

  if($min_start >= $slice->start() && $max_end <= $slice->end()) {
    $ext_slice = $slice;
  } else {
    my $sa = $self->db()->get_SliceAdaptor();
    $ext_slice = $sa->fetch_by_region
      ($slice->coord_system->name(), $slice->seq_region_name(),
       $min_start,$max_end, $slice->strand(), $slice->coord_system->version());
  }

  # associate transcript identifiers with genes

  my %g_hash = map {$_->dbID => $_} @$genes;

  my $g_id_str = '(' . join(',', keys %g_hash) . ')';

  my $sth = $self->prepare("SELECT gene_id, transcript_id " .
                           "FROM   transcript " .
                           "WHERE  gene_id IN $g_id_str");

  $sth->execute();

  my ($g_id, $tr_id);
  $sth->bind_columns(\$g_id, \$tr_id);

  my %tr_g_hash;

  while($sth->fetch()) {
    $tr_g_hash{$tr_id} = $g_hash{$g_id};
  }

  $sth->finish();

  my $ta = $self->db()->get_TranscriptAdaptor();
  my $transcripts = $ta->fetch_all_by_Slice($ext_slice,1);

  # move transcripts onto gene slice, and add them to genes
  foreach my $tr (@$transcripts) {
    if( !exists $tr_g_hash{$tr->dbID()} ) {
      next;
    }

    my $new_tr;
    if($slice != $ext_slice) {
      $new_tr = $tr->transfer($slice) if($slice != $ext_slice);
      if(!$new_tr) {
	throw("Unexpected. Transcript could not be transfered onto Gene slice.");
      }
    } else {
      $new_tr = $tr;
    }


    $tr_g_hash{$tr->dbID()}->add_Transcript($new_tr);
  }

  return $genes;
}






=head2 fetch_by_transcript_id

  Arg [1]    : int $transid
               unique database identifier for the transcript whose gene should
               be retrieved. The gene is returned in its native coord
               system (i.e. the coord_system it is stored in). If the coord
               system needs to be changed, then tranform or transfer should
               be called on the returned object.  undef is returned if the
               gene or transcript is not found in the database.
  Example    : $gene = $gene_adaptor->fetch_by_transcript_id( 1241 );
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub fetch_by_transcript_id {
    my ( $self, $trans_id ) = @_;

    # this is a cheap SQL call
     my $sth = $self->prepare("SELECT	tr.gene_id " .
				"FROM	transcript as tr " .
				"WHERE	tr.transcript_id = ?");
    $sth->execute($trans_id);

    my ($geneid) = $sth->fetchrow_array();

    $sth->finish();

    return undef if( !defined $geneid );

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
  Status     : Stable

=cut

sub fetch_by_transcript_stable_id {
    my ( $self, $trans_stable_id) = @_;

    # this is a cheap SQL call
    my $sth = $self->prepare(
        "SELECT	tr.gene_id " .
				"FROM	  transcript as tr, transcript_stable_id tcl " .
        "WHERE	tcl.stable_id = ? " .
        "AND    tr.transcript_id = tcl.transcript_id");
    $sth->execute("$trans_stable_id");

    my ($geneid) = $sth->fetchrow_array();
    $sth->finish;
    if( !defined $geneid ) {
        return undef;
    }
    my $gene = $self->fetch_by_dbID( $geneid );
    return $gene;
}





=head2 fetch_by_translation_stable_id

  Arg [1]    : string $translation_stable_id
               the stable id of a translation of the gene that should
               be obtained
  Example    : $gene = $gene_adaptor->fetch_by_translation_stable_id
                 ( 'ENSP00000278194' );
  Description: retrieves a gene via the stable id of one of its translations.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : geneview
  Status     : Stable

=cut

sub fetch_by_translation_stable_id {
    my ( $self, $translation_stable_id ) = @_;

    # this is a cheap SQL call
    my $sth = $self->prepare
      ("SELECT	tr.gene_id " .
			 "FROM	transcript as tr, " .
       "      translation as tl, " .
			 "		  translation_stable_id as trs " .
			 "WHERE	trs.stable_id = ? " .
			 "AND	  trs.translation_id = tl.translation_id " .
       "AND   tr.transcript_id = tl.transcript_id");

    $sth->execute("$translation_stable_id");

    my ($geneid) = $sth->fetchrow_array();
    $sth->finish;
    if( !defined $geneid ) {
        return undef;
    }
    return $self->fetch_by_dbID($geneid);
}



=head2 fetch_all_by_external_name

  Arg [1]    : string $external_id
               the external identifier for the gene to be obtained
  Example    : @genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA2')}
  Description: retrieves a list of genes with an external database 
               idenitifier $external_id.  The genes returned are in their
               native coordinate system.  I.e. in the coordinate system they
               are stored in the database in.  If another coordinate system
               is required then the Gene::transfer or Gene::transform method 
               can be used.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : goview
  Status     : Stable

=cut

sub fetch_all_by_external_name {
  my ( $self, $external_id) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my ( @ids, @result );
  @ids = $entryAdaptor->list_gene_ids_by_extids( $external_id );

  my $genes = $self->fetch_all_by_dbID_list(\@ids);

  my %genes_by_dbIDs = map { $_->dbID(),$_ } @$genes;
  @result = map { $genes_by_dbIDs{ $_ } } @ids;

  return \@result;
}



=head2 fetch_all_alt_alleles

  Arg [1]    : Bio::EnsEMBL::Gene $gene
  Example    : my @alt_genes = @{$gene_adaptor->fetch_all_alt_alleles($gene);}
               foreach my $alt_gene (@alt_genes) {
                 print "Alternate allele: " . $alt_gene->stable_id() . "\n";
               }
  Description: Retrieves genes which are alternate alleles to a provided gene.
               Alternate alleles in Ensembl are genes which are similar and are
               on an alternative haplotype of the same region.  There are not 
               currently very many of these.  This method will return a 
               reference to an empty list if no alternative alleles are found.
  Returntype : reference to a list of genes
  Exceptions : throw if incorrect arg provided
               warning if gene arg does not have dbID
  Caller     : Gene::get_all_alt_alleles
  Status     : Stable

=cut

sub fetch_all_alt_alleles {
  my $self = shift;
  my $gene = shift;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('Bio::EnsEMBL::Gene argument is required');
  }

  my $gene_id = $gene->dbID();

  if(!$gene_id) {
    warning('Cannot retrieve alternate alleles for gene without dbID');
    return [];
  }

  my $sth = $self->prepare("SELECT aa1.gene_id " .
                           "FROM   alt_allele aa1, alt_allele aa2 " .
                           "WHERE  aa1.alt_allele_id = aa2.alt_allele_id " .
                           "AND    aa2.gene_id = ? " .
                           "AND    aa1.gene_id <> ?");

  $sth->execute($gene_id, $gene_id);

  my @alt_ids;
  my $row;
  while($row = $sth->fetchrow_arrayref()) {
    push @alt_ids, $row->[0];
  } 
  $sth->finish();
  
  if(@alt_ids) {
    return $self->fetch_all_by_dbID_list(\@alt_ids);
  }
  
  return [];
  
}




=head2 store_alt_alleles


  Arg [1]    : reference to list of Bio::EnsEMBL::Genes $genes
  Example    : $gene_adaptor->store_alt_allele([$gene1, $gene2, $gene3]);
  Description: This method creates a group of alternative aleles (i.e. locus)
               from a set of genes.  The genes should be genes from alternate
               haplotypes which are similar.  The genes must already be stored
               in this database. At least 2 genes must be in the list reference
               provided.
  Returntype : none
  Exceptions : throw on incorrect arguments
               throw on sql error (e.g. duplicate unique id)
  Caller     : ?
  Status     : Stable

=cut

sub store_alt_alleles {
  my $self = shift;
  my $genes = shift;

  if(!ref($genes) eq 'ARRAY') {
    throw('List reference of Bio::EnsEMBL::Gene argument expected.');
  }

  my $num_genes = scalar(@$genes);

  if($num_genes < 2) {
    throw("At least 2 genes must be provided to construct alternate alleles.");
  }

  return if(!@$genes);
  
  #
  #insert the first gene seperately in order to get a unique identifier for
  #the set of alleles
  #
  my $gene = $genes->[0];

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('List reference of Bio::EnsEMBL::Gene argument expected.');
  }

  my $gene_id = $gene->dbID();

  if(!$gene_id) {
    throw("Genes must have dbIDs in order to construct alternate alleles.");
  }

  my $sth = $self->prepare("INSERT INTO alt_allele (gene_id) VALUES (?)");
  $sth->execute($gene->dbID());
  
  my $alt_allele_id = $sth->{'mysql_insertid'};
  $sth->finish();

  #
  # Insert all subsequent alt alleles using the alt_allele identifier
  # from the first insert
  #

  $sth = $self->prepare("INSERT INTO alt_allele (alt_allele_id, gene_id) " .
                        "VALUES (?,?)");
  
  for(my $i = 1; $i < $num_genes; $i++) {
    my $gene = $genes->[$i];

    if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
      throw("List reference of Bio::EnsEMBL::Gene argument expected"); 
    }
    
    $gene_id = $gene->dbID();
    
    if(!$gene_id) {
      #This is an error but we have already inserted into the database
      #delete the already inserted entries to restore the state of the
      #database
      $sth->finish();
      $sth->prepare("DELETE FROM alt_allele WHERE alt_allele_id = ?");
      $sth->execute($alt_allele_id);
      $sth->finish();
      throw('Genes must have dbIDs in order to construct alternate alleles.');
    }

    eval {
      $sth->execute($alt_allele_id, $gene_id);
    };

    if($@) {
      #an error occured, revert the db to the previous state
      $sth = $self->prepare("DELETE FROM alt_allele WHERE alt_allele_id = ?");
      $sth->execute($alt_allele_id);
      $sth->finish();
      throw("An SQL error occured inserting alternate alleles:\n$@");
    }
  }
  
  $sth->finish();

  return;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $gene_adaptor->store($gene);
  Description: Stores a gene in the database
  Returntype : the database identifier of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene or if 
               $gene does not have an analysis object
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ( $self, $gene ) = @_;

  if(!ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
    throw("Must store a gene object, not a $gene");
  }

  my $db = $self->db();

  if($gene->is_stored($db)) {
    return $gene->dbID();
  }

  # ensure coords are correct before storing
  $gene->recalculate_coordinates();

  my $analysis = $gene->analysis();
  throw("Genes must have an analysis object.") if(!defined($analysis));

  my $analysis_id;
  if($analysis->is_stored($db)) {
    $analysis_id = $analysis->dbID();
  } else {
    $analysis_id = $db->get_AnalysisAdaptor->store($analysis);
  }

  my $type = $gene->biotype || "";

  my $original = $gene;
  my $original_transcripts = $gene->get_all_Transcripts();
  my $seq_region_id;
  ($gene, $seq_region_id) = $self->_pre_store($gene);

  my $store_gene_sql =
        "INSERT INTO gene " .
           "SET biotype = ?, " .
               "analysis_id = ?, " .
               "seq_region_id = ?, " .
               "seq_region_start = ?, " .
               "seq_region_end = ?, " .
               "seq_region_strand = ?, ".
	       "description = ?, " .
               "source = ?, ".
               "status = ? ";

  my $sth = $self->prepare( $store_gene_sql );
   $sth->execute(
		 $type,
		 $analysis_id,
		 $seq_region_id,
		 $gene->start(),
		 $gene->end(),
		 $gene->strand(),
		 $gene->description(),
		 $gene->source(),
		 $gene->status()
		);
  $sth->finish();

   my $gene_dbID = $sth->{'mysql_insertid'};

   # store stable ids if they are available
   if (defined($gene->stable_id)) {

     my $statement = "INSERT INTO gene_stable_id
                         SET gene_id = ?,
                             stable_id = ?,
                             version = ?, ";
     if( $gene->created_date() ) {
       $statement .= "created_date = from_unixtime( ".$gene->created_date()."),";
     } else {
       $statement .= "created_date = \"0000-00-00 00:00:00\",";
     }

     if( $gene->modified_date() ) {
       $statement .= "modified_date = from_unixtime( ".$gene->modified_date().")";
     } else {
       $statement .= "modified_date = \"0000-00-00 00:00:00\"";
     }

     $sth = $self->prepare($statement);
     $sth->execute( $gene_dbID, $gene->stable_id(), $gene->version());
     $sth->finish();
   }



   # store the dbentries associated with this gene
   my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

   foreach my $dbe ( @{$gene->get_all_DBEntries} ) {
     $dbEntryAdaptor->store( $dbe, $gene_dbID, "Gene" );
   }


   # we allow transcripts not to share equal exons and instead have copies
   # For the database we still want sharing though, to have easier time with
   # stable ids. So we need to have a step to merge exons together before store
   my %exons;

   foreach my $trans ( @{$gene->get_all_Transcripts} ) {
     foreach my $e ( @{$trans->get_all_Exons} ) {
       my $key = $e->hashkey();

       if( exists $exons{ $key } ) {
         $trans->swap_exons( $e, $exons{$key} );
       } else {
         $exons{$key} = $e;
       }
     }
   }

  my $transcript_adaptor = $db->get_TranscriptAdaptor();

  my $transcripts = $gene->get_all_Transcripts();

  for(my $i = 0; $i < @$transcripts; $i++) {
    my $new = $transcripts->[$i];
    my $old = $original_transcripts->[$i];

    $transcript_adaptor->store($new, $gene_dbID );

    # update the original transcripts since we may have made copies of
    # them by transforming the gene
    $old->dbID($new->dbID());
    $old->adaptor($new->adaptor());
    if($new->translation) {
      $old->translation->dbID($new->translation()->dbID);
      $old->translation->adaptor($new->translation()->adaptor);
    }
  }

  # update gene to point to display xref if it is set
  if(my $display_xref = $gene->display_xref) {
    my $dxref_id;

    if($display_xref->is_stored($db)) {
      $dxref_id = $display_xref->dbID();
    } else {
      $dxref_id = $dbEntryAdaptor->exists($display_xref);
    }

    if(defined($dxref_id)) {
      $sth = $self->prepare
        ("UPDATE gene SET display_xref_id = ? WHERE gene_id = ?");
      $sth->execute($dxref_id, $gene_dbID);
      $sth->finish();
      $display_xref->dbID($dxref_id);
      $display_xref->adaptor($dbEntryAdaptor);
      $display_xref->dbID($dxref_id);
      $display_xref->adaptor($dbEntryAdaptor);
    } else {
      warning("Display_xref ".$display_xref->dbname().":".
              $display_xref->display_id() . " is not stored in database.\n".
              "Not storing relationship to this gene.");
      $display_xref->dbID(undef);
      $display_xref->adaptor(undef);
    }
  }

  # set the adaptor and dbID on the original passed in gene not the
  # transfered copy
  $original->adaptor( $self );
  $original->dbID( $gene_dbID );

  # store gene attributes if there are any
  my $attr_adaptor = $db->get_AttributeAdaptor();
  $attr_adaptor->store_on_Gene($gene,
                               $gene->get_all_Attributes);

  return $gene_dbID;
}




=head2 remove

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               the gene to remove from the database
  Example    : $gene_adaptor->remove($gene);
  Description: Removes a gene completely from the database.  All associated
               transcripts, exons, stable_identifiers, descriptions, etc.
               are removed as well. Use with caution!
  Returntype : none
  Exceptions : throw on incorrect arguments 
               warning if gene is not stored in this database
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $gene = shift;

  if(!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  if( !$gene->is_stored($self->db()) ) {
    warning("Cannot remove gene " . $gene->dbID() . ". Is not stored in " .
            "this database.");
    return;
  }

  # remove all object xrefs associated with this gene

  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
  foreach my $dbe (@{$gene->get_all_DBEntries()}) {
    $dbe_adaptor->remove_from_object($dbe, $gene, 'Gene');
  }

  # remove all alternative allele entries associated with this gene
  my $sth = $self->prepare("delete from alt_allele where gene_id = ?");
  $sth->execute($gene->dbID());
  $sth->finish();

  # remove the attributes associated with this transcript
  my $attrib_adaptor = $self->db->get_AttributeAdaptor;  
  $attrib_adaptor->remove_from_Gene($gene);

  # remove all of the transcripts associated with this gene

  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
    $transcriptAdaptor->remove($trans);
  }

  # remove the gene stable identifier

  $sth = $self->prepare( "delete from gene_stable_id where gene_id = ? " );
  $sth->execute( $gene->dbID );
  $sth->finish();

  # remove this gene from the database

  $sth = $self->prepare( "delete from gene where gene_id = ? " );
  $sth->execute( $gene->dbID );
  $sth->finish();

  # unset the gene identifier and adaptor thereby flagging it as unstored

  $gene->dbID(undef);
  $gene->adaptor(undef);

  return;
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
  Status     : Stable

=cut

sub get_Interpro_by_geneid {
   my ($self,$gene) = @_;
    my $sql="
	SELECT	i.interpro_ac, 
		x.description 
  FROM	transcript t,
    translation tl, 
		protein_feature pf, 
		interpro i, 
    xref x,
		gene_stable_id gsi
	WHERE	gsi.stable_id = '$gene' 
	    AND	t.gene_id = gsi.gene_id
      AND tl.transcript_id = t.transcript_id
	    AND	tl.translation_id = pf.translation_id 
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




=head2 update

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $gene_adaptor->update($gene);
  Description: Updates the type, analysis and display_xref of a gene in the
               database.
  Returntype : None
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene
  Caller     : general
  Status     : Stable

=cut

sub update {
   my ($self,$gene) = @_;
   my $update = 0;

   if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     throw("Must update a gene object, not a $gene");
   }

   my $update_gene_sql = "
        UPDATE gene
           SET biotype = ?,
               analysis_id = ?,
               display_xref_id = ?,
               status = ?
         WHERE gene_id = ?";

   my $display_xref = $gene->display_xref();
   my $display_xref_id;

   if( $display_xref && $display_xref->dbID() ) {
     $display_xref_id = $display_xref->dbID();
   } else {
     $display_xref_id = undef;
   }

   my $sth = $self->prepare( $update_gene_sql );
   $sth->execute($gene->type(), 
		 $gene->analysis->dbID(),
		 $display_xref_id,
		 $gene->status(),
		 $gene->dbID()
		);

   # maybe should update stable id or gene description ???
}


# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Example    : none 
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Genes
#  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();
  my $dbEntryAdaptor = $self->db()->get_DBEntryAdaptor();

  my @genes;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my ( $gene_id, $seq_region_id, $seq_region_start, $seq_region_end, 
       $seq_region_strand, $analysis_id, $biotype, $display_xref_id, 
       $gene_description, $stable_id, $version, $created_date, 
       $modified_date, $xref_display_id, $status, $source, 
       $xref_primary_acc, $xref_desc, $xref_version, $external_name, 
       $external_db, $external_status, $external_release, $external_db_name );

  $sth->bind_columns( \$gene_id, \$seq_region_id, \$seq_region_start,
          \$seq_region_end, \$seq_region_strand, \$analysis_id, \$biotype,
          \$display_xref_id, \$gene_description, \$status, \$source, 
	  \$stable_id, \$version,
          \$created_date, \$modified_date, 
	  \$xref_display_id, \$xref_primary_acc, \$xref_desc, \$xref_version,
          \$external_db, \$external_status,
          \$external_release, \$external_db_name );

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;

  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_sr_name;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
  }

  FEATURE: while($sth->fetch()) {
    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

    #
    # remap the feature coordinates to another coord system 
    # if a mapper was provided
    #
    if($mapper) {

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
			 $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
          $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
                               $cmp_cs_vers);
      } else {
        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
                               $asm_cs_vers);
      }
    }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice) {
      if($dest_slice_start != 1 || $dest_slice_strand != 1) {
        if($dest_slice_strand == 1) {
          $seq_region_start = $seq_region_start - $dest_slice_start + 1;
          $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start = $dest_slice_end - $seq_region_end + 1;
          $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
          $seq_region_strand *= -1;
        }

      }

      #throw away features off the end of the requested slice or on different seq_region
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	 ( $dest_slice_sr_name ne $sr_name )) {
	next FEATURE;
      }

      $slice = $dest_slice;
    }

    my $display_xref;

    if( $display_xref_id ) {
      $display_xref = Bio::EnsEMBL::DBEntry->new_fast
        ({ 'dbID' => $display_xref_id,
           'adaptor' => $dbEntryAdaptor,
           'display_id' => $xref_display_id,
           'primary_id' => $xref_primary_acc,
           'version'    => $xref_version,
           'description' => $xref_desc,
           'release' => $external_release,
           'dbname' => $external_db,
           'db_display_name' => $external_db_name
         });
      $display_xref->status( $external_status );
    }				

    #finally, create the new gene
    push @genes, Bio::EnsEMBL::Gene->new
      ( '-analysis'      =>  $analysis,
        '-biotype'       =>  $biotype,
        '-start'         =>  $seq_region_start,
        '-end'           =>  $seq_region_end,
        '-strand'        =>  $seq_region_strand,
        '-adaptor'       =>  $self,
        '-slice'         =>  $slice,
        '-dbID'          =>  $gene_id,
        '-stable_id'     =>  $stable_id,
        '-version'       =>  $version,
	'-created_date'  =>  $created_date || undef,
	'-modified_date' =>  $modified_date || undef,
        '-description'   =>  $gene_description,
        '-external_name' =>  $external_name,
        '-external_db'   =>  $external_db,
        '-external_status' => $external_status,
        '-display_xref' => $display_xref,
	'-status'   => $status,
        '-source'       => $source );
  }

  return \@genes;
}


=head2 cache_gene_seq_mappings

  Args       : none
  Example    : $gene_adaptor->cache_gene_seq_mappings();
  Description: caches all the assembly mappings needed for genes
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : At Risk
             : New experimental code.

=cut

sub cache_gene_seq_mappings{
  my ($self) = @_;

  #get the sequence level to map too

  my $sql = (<<SSQL);
  SELECT	name 
   FROM	coord_system 
    WHERE attrib like "%sequence_level%"
SSQL
  my $sth = $self->prepare($sql);
  $sth->execute();
  
  my $sequence_level = $sth->fetchrow_array();
  
  $sth->finish();

  my $csa = $self->db->get_CoordSystemAdaptor();
  my $ama = $self->db->get_AssemblyMapperAdaptor();

  my $cs1 = $csa->fetch_by_name($sequence_level);
  

  #get level to map to two

  my $mcc =  $self->db->get_MetaCoordContainerAdaptor();
  my $csnew = $mcc->fetch_all_CoordSystems_by_feature_type('gene');

  foreach my $cs2 (@$csnew){
    my $am = $ama->fetch_by_CoordSystems($cs1,$cs2);
    $am->register_all();    
  };
}
  
##########################
#                        #
#  DEPRECATED METHODS    #
#                        #
##########################


=head2 fetch_by_maximum_DBLink

 Description: DEPRECATED - use fetch_all_by_external_name instead


=cut

sub fetch_by_maximum_DBLink {
  my ( $self, $external_id ) = @_;
  
  deprecate( "use fetch_all_by_external_name instead" );

  my $genes=$self->fetch_all_by_external_name( $external_id);
  
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



=head2 get_display_xref

  Description: DEPRECATED use $gene->display_xref

=cut

sub get_display_xref {
  my ($self, $gene ) = @_;

  deprecate( "display xref should retrieved from Gene object directly" );

  if( !defined $gene ) {
      throw("Must call with a Gene object");
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


=head2 get_description

  Description: DEPRECATED, use gene->get_description

=cut

sub get_description {
  my ($self, $dbID) = @_;

  deprecate( "Gene description should be loaded on gene retrieval. Use gene->get_description()" );

  if( !defined $dbID ) {
      throw("must call with dbID");
  }

  my $sth = $self->prepare("SELECT description 
                            FROM   gene_description 
                            WHERE  gene_id = ?");
  $sth->execute($dbID);
  my @array = $sth->fetchrow_array();
  return $array[0];
}




=head2 fetch_by_Peptide_id

  Description: DEPRECATED, use fetch_by_translation_stable_id()

=cut

sub fetch_by_Peptide_id {
    my ( $self, $translation_stable_id) = @_;

    deprecate( "Please use better named fetch_by_translation_stable_id \n".caller(2) );

    $self->fetch_by_translation_stable_id($translation_stable_id);
}


=head2 get_stable_entry_info

  Description: DEPRECATED use $gene->stable_id instead

=cut

sub get_stable_entry_info {
  my ($self,$gene) = @_;

  deprecated( "stable id info is loaded on default, no lazy loading necessary" );

  if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     throw("Needs a gene object, not a $gene");
  }

  my $sth = $self->prepare("SELECT stable_id, UNIX_TIMESTAMP(created),
                                   UNIX_TIMESTAMP(modified), version 
                            FROM gene_stable_id 
                            WHERE gene_id = ".$gene->dbID);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $gene->{'stable_id'} = $array[0];
  $gene->{'created'}   = $array[1];
  $gene->{'modified'}  = $array[2];
  $gene->{'version'}   = $array[3];

  return 1;
}



=head2 fetch_all_by_DBEntry

  Description: DEPRECATED - Use fetch_all_by_external_name instead

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;
  deprecate('This method has been deprecated because there was another.' .
            "Method which did exactly the same thing.\n" .
            'Use fetch_all_by_external_name instead.');
  return $self->fetch_all_by_external_name(@_);
}



1;

