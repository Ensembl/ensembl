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
  
  return ([ 'gene', 'g' ],[ 'gene_description', 'gd' ], [ 'gene_stable_id', 'gsi' ],
	  [ 'xref', 'x' ], [ 'external_db' , 'exdb' ]);
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

  return qw( g.gene_id g.seq_region_id g.seq_region_start g.seq_region_end 
	     g.seq_region_strand g.analysis_id g.type g.display_xref_id 
	     gd.description gsi.stable_id gsi.version x.display_label 
	     exdb.db_name exdb.status );
}


sub _left_join {
  return ( [ 'gene_description', "gd.gene_id = g.gene_id" ],
	   [ 'gene_stable_id', "gsi.gene_id = g.gene_id" ],
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

sub fetch_by_stable_id {
   my ($self,$id, $cs_name, $cs_version) = @_;

   my $constraint = "gsi.stable_id = \"$id\"";

   # should be only one :-)
   my $genes = $self->SUPER::generic_fetch( $constraint );

   if( ! @$genes ) { return undef }

   my @new_genes = map { $_->transform( $cs_name, $cs_version ) } @$genes;


   return $new_genes[0];
 }


=head2 fetch_by_exon_stable_id

  Arg [1]    : string $id
               The stable id of the gene to retrieve
  Arg [2]    : (optional) string $coordinate_system_name
  Arg [3]    : (optional) string $coordinate_system_version
  Example    : $gene = $gene_adaptor->fetch_by_exon_stable_id('ENSE00000148944', "chromosome", "NCBI_33" );
  Description: Retrieves a gene object from the database via an exon stable id
  Returntype : Bio::EnsEMBL::Gene in given coordinate system
  Exceptions : throw ?
  Caller     : general

=cut

sub fetch_by_exon_stable_id{
   my ( $self, $id, $cs_name, $cs_version ) = @_;

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

   my $gene = $self->fetch_by_dbID( $dbID, $cs_name, $cs_version);
  
   return $gene;
 }




=head2 fetch_all_by_domain

  Arg [1]    : string $domain
               the domain to fetch genes from
  Example    : my @genes = $gene_adaptor->fetch_all_by_domain($domain);
  Description: retrieves a listref of genes whose translation contain interpro
               domain $domain.
  Returntype : list of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : domainview

=cut

sub fetch_all_by_domain {
  my ($self, $domain, $cs_name, $cs_version) = @_;

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

  my $constraint = "g.gene_id in (".join( ",", @gene_ids ).")";
  my $genes = $self->SUPER::generic_fetch( $constraint );
  my @new_genes = map { $_->transform( $cs_name, $cs_version ) } @$genes;

  return \@new_genes;
}


  

=head2 fetch_by_transcript_id

  Arg [1]    : int $transid 
               unique database identifier for the transcript whose gene should
               be retrieved.
  Arg [2]    : string $coordinate_system_name
  Arg [3]    : string $coordinate_system_version
  Example    : $gene = $gene_adaptor->fetch_by_transcript_id( $transcript, "chromosome" );
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : ?

=cut

sub fetch_by_transcript_id {
    my ( $self, $trans_id, $cs_name, $cs_version ) = @_;

    
    # this is a cheap SQL call
     my $sth = $self->prepare("	SELECT	tr.gene_id 
				FROM	transcript as tr 
				WHERE	tr.transcript_id = $trans_id");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    my $gene = $self->fetch_by_dbID( $geneid, $cs_name, $cs_version );
    return $gene;
}

=head2 fetch_by_transcript_stable_id

  Arg [1]    : string $transid 
               unique database identifier for the transcript whose gene should
               be retrieved.
  Arg [2]    : string $coordinate_system_name
  Arg [3]    : string $coordinate_system_version
  Example    : none
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : ?

=cut

sub fetch_by_transcript_stable_id {
    my ( $self, $trans_stable_id, $cs_name, $cs_version ) = @_;

    # this is a cheap SQL call
     my $sth = $self->prepare("	SELECT	tr.gene_id 
				FROM	transcript as tr, transcript_stable_id tcl 
				WHERE	tcl.stable_id = \"$trans_stable_id\"
                                  AND   tr.transcript_id = tcl.transcript_id");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    my $gene = $self->fetch_by_dbID( $geneid, $cs_name, $cs_version );
    return $gene;
}



=head2 fetch_by_Peptide_id

  Arg [1]    : 
  Example    : 
  Description: DEPRECATED, see below 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub fetch_by_Peptide_id {
    my ( $self, $translation_stable_id, $cs_name, $cs_version ) = @_;

    deprecate( "Please use better named fetch_by_translation_stable_id" );

    $self->fetch_by_translation_stable_id
      ( $translation_stable_id, 
	$cs_name, $cs_version 
      );
}


=head2 fetch_by_translation_stable_id

  Arg [1]    : string $translation_stable_id
               the stable id of a translation of the gene that should
               be obtained
  Example    : $gene = $gene_adaptor->fetch_by_translation_stable_id
                 ( 'ENSP00000278194', "toplevel" );
  Description: retrieves a gene via the stable id of one of its translations.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : geneview

=cut

sub fetch_by_translation_stable_id {
    my ( $self, $translation_stable_id, $cs_name, $cs_version ) = @_;

    # this is a cheap SQL call
    my $sth = $self->prepare("	SELECT	tr.gene_id 
				FROM	transcript as tr, 
					translation_stable_id as trs 
			    	WHERE	trs.stable_id = '$translation_stable_id' 
				AND	trs.translation_id = tr.translation_id
			    ");
    $sth->execute;

    my ($geneid) = $sth->fetchrow_array();
    if( !defined $geneid ) {
        return undef;
    }
    return $self->fetch_by_dbID($geneid);
}  



=head2 get_stable_entry_info

  Arg [1]    : Bio::EnsEMBL::Gene $gene
  Example    : $gene_adaptor->get_stable_entry_info($gene);
  Description: DEPRECATED
               gets stable info for a gene. this is not usually done at
               creation time for speed purposes, and can be lazy-loaded later
               if it is needed..
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::Gene

=cut

sub get_stable_entry_info {
  my ($self,$gene) = @_;

  deprecated( "stable id info is loaded on default, no lazy loading necessary" );

  if( !defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene') ) {
     $self->throw("Needs a gene object, not a $gene");
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
  my ( $self, $external_id, $cs_name, $cs_version ) = @_;

  my @genes = ();

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();


  my @ids = $entryAdaptor->list_gene_ids_by_extids( $external_id );
  foreach my $gene_id ( @ids ) {
    my $gene = $self->fetch_by_dbID( $gene_id, $cs_name, $cs_version );
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
   my ( $self, $gene ) = @_;
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
 
   my $type = "";

   if (defined($gene->type)) {
       $type = $gene->type;
   }

   my $store_gene_sql = "
        INSERT INTO gene
           SET type = ?,
               analysis_id = ?,
               display_xref_id = ?,
               seq_region_id = ?,
               seq_region_start = ?,
               seq_region_end = ?,
               seq_region_strand = ? ";

   my $sth2 = $self->prepare( $store_gene_sql );
   $sth2->execute(
		  "$type", 
		  $analysisId, 
		  $gene->display_xref()->dbID(),
		  $gene->slice->seq_region_id(),
		  $gene->start(),
		  $gene->end(),
		  $gene->strand()
		  );

   my $gene_dbID = $sth2->{'mysql_insertid'};

   #
   # store stable ids if they are available
   #
   if (defined($gene->stable_id)) {

     my $statement = "INSERT INTO gene_stable_id
                         SET gene_id = ?,
                             stable_id = ?,
                             version = ?,
                             created = NOW(),
                             modified = NOW()";
     my $sth = $self->prepare($statement);
     $sth->execute( $gene_dbID, $gene->stable_id(), $gene->version() );
   }


   #
   # store the dbentries associated with this gene
   #
   my $dbEntryAdaptor = $self->db->get_DBEntryAdaptor();

   foreach my $dbe ( @{$gene->get_all_DBEntries} ) {
     $dbEntryAdaptor->store( $dbe, $gene_dbID, "Gene" );
   }

   # we allow transcripts not to share equal exons and instead have copies
   # For the database we still want sharing though, to have easier time with
   # stable ids. So we need to have a step to merge exons together before store.
   my %exons;

   foreach my $trans ( @{$self->get_all_Transcripts} ) {
     foreach my $e ( @{$trans->get_all_Exons} ) {
       my $key = $e->start()."-".$e->end()."-".$e->strand()."-".$e->phase();

       if( exists $exons{ $key } ) {
	 $trans->swap_exons( $e, $exons{$key} );
       } else {
	 $exons{$key} = $e;
       }
     }
   }

   foreach my $t ( @{$self->get_all_Transcripts()} ) {
     $transcriptAdaptor->store($t,$gene_dbID );
   }

   $gene->adaptor( $self );
   $gene->dbID( $gene_dbID );

   # should we store gene description? Not strictly necessary, as its later 
   # added to the database. But maybe just for completeness and it should definatly 
   # go into update.
   
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
  $sth = $self->prepare( "delete from gene_description where gene_id = ?" );
  $sth->execute( $gene->dbID() );

  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
    $transcriptAdaptor->remove($trans,$gene);
  }

  $gene->dbID(undef);
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

   my $update_gene_sql = "
        UPDATE gene
           SET type = ?,
               analysis_id = ?,
               display_xref_id = ?,
               seq_region_id = ?,
               seq_region_start = ?,
               seq_region_end = ?,
               seq_region_strand = ?
         WHERE gene_id = ?";

   my $display_xref = $gene->display_xref();
   my $display_xref_id;
 
   if( defined $display_xref && $display_xref->dbID() ) {
     $display_xref_id = $display_xref->dbID();
   } else {
     $display_xref_id = undef;
   }

   my $sth = $self->prepare( $update_gene_sql );
   $sth->execute(
	         $gene->type(), 
		 $gene->analysis->dbID(),
		 $display_xref_id,
		 $gene->slice->seq_region_id(),
		 $gene->start(),
		 $gene->end(),
		 $gene->strand(),
		 
		 $gene->dbID()
		);
   
   # maybe should update stable id or gene description ???
}


=head2 _obj_from_hashref

  Arg [1]    : Hashreference $hashref
  Example    : none 
  Description: PROTECTED implementation of abstract superclass method.
               responsible for the creation of Genes 
  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
  Exceptions : none
  Caller     : internal

=cut

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
  my %rc_hash;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my ( $gene_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand,
       $analysis_id, $type, $display_xref_id, $gene_description, $stable_id, 
       $version, $external_name, $external_db, $external_status );  


  $sth->bind_columns( \$gene_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, 
		      \$seq_region_strand, \$analysis_id, \$type, \$display_xref_id, 
		      \$gene_description, \$stable_id, \$version, \$external_name, 
		      \$external_db, \$external_status );



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
    $asm_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }

  FEATURE: while($sth->fetch()) {
    #get the analysis object
    my $analysis = $aa->fetch_by_dbID($analysis_id);

    my $slice = $slice_hash{$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #
    # remap the feature coordinates to another coord system 
    # if a mapper was provided
    #
    if($mapper) {
      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
			 $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($asm_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
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
    if($dest_slice && ($dest_slice_start != 1 || $dest_slice_strand != 1)) {
      if($dest_slice_strand == 1) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
      } else {
        my $tmp_seq_region_start = $seq_region_start;
        $seq_region_start = $dest_slice_end - $seq_region_end + 1;
        $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
        $seq_region_strand *= -1;
      }

      $slice = $dest_slice;

      #throw away features off the end of the requested slice
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
        next FEATURE;
      }
    }

    my $display_xref;

    if( $display_xref_id ) {
      $display_xref = bless 
	{ 'dbID' => $display_xref_id,
	  'adaptor' => $dbEntryAdaptor
	}, "Bio::EnsEMBL::DBEntry";
    }
				

    #finally, create the new repeat feature
    push @genes, Bio::EnsEMBL::Gene->new
      ( '-analysis'      =>  $analysis,
	'-start'         =>  $seq_region_start,
	'-end'           =>  $seq_region_end,
	'-strand'        =>  $seq_region_strand,
	'-adaptor'       =>  $self,
	'-slice'         =>  $slice,
	'-dbID'          =>  $gene_id,
        '-stable_id'     =>  $stable_id,
        '-version'       =>  $version,
        '-description'   =>  $gene_description,
	'-external_name' =>  $external_name,
        '-external_db'   =>  $external_db,
        '-external_status' => $external_status,
	'-display_xref' => $display_xref );
  }

  return \@genes;
}

=head2 fetch_by_maximum_DBLink

 Title   : fetch_by_maximum_DBLink
 Usage   : $geneAdaptor->fetch_by_maximum_DBLink($ext_id)
 Function: gets one gene out of the db with 
 Returns : gene object (with transcripts, exons)
 Args    : 


=cut

sub fetch_by_maximum_DBLink {
  my ( $self, $external_id, $cs_name, $cs_version ) = @_;
  
  my $genes=$self->fetch_all_by_DBEntry( $external_id, $cs_name, $cs_version );
  
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

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               retrieve display_xref for this gene
  Example    : none
  Description: DEPRECATED
               Retrieves the display_xref for a gene.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : thrown if $dbId arg is not defined
  Caller     : general

=cut

sub get_display_xref {
  my ($self, $gene ) = @_;

  deprecate( "display xref should retrieved from Gene object directly" );

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

=head2 get_description

  Arg [1]    : int $dbID
               the database identifier of the gene whose description should be
               obtained
  Example    : $description = $gene_adaptor->get_description(1234);
  Description: DEPRECATED, use gene->get_description
               Retrieves a description for a gene which is created during
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

  deprecate( "Gene description should be loaded on gene retrieval. Use gene->get_description()" );

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




1;
__END__

