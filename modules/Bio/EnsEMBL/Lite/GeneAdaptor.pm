# EnsEMB Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 15.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneLiteAdaptor - MySQL Database queries to retrieve 
genes quickly from denormalized tables.

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::Lite::GeneAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Utils::Cache; #CPAN LRU Cache

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

my $SLICE_GENE_CACHE_SIZE = 3;
my $MAX_TRANSCRIPT_LENGTH = 3000000;


=head2 new

  Arg [1]    : list of arguments @args
  Example    : $gene_adaptor = new Bio::EnsEMBL::DBSQL::GeneAdaptor($db_adptr);
  Description: Craetes a new GeneAdaptor object
  Returntype : Bio::EnsEMBL::DBSQL::GeneAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my($class, @args) = @_;

  #call superclass constructor
  my $self = $class->SUPER::new(@args);

  #initialize tie hash cache
  tie (%{$self->{'_slice_gene_cache'}}, 
       'Bio::EnsEMBL::Utils::Cache', 
       $SLICE_GENE_CACHE_SIZE);

  #use another cache for 'empty' genes
  tie (%{$self->{'_slice_empty_gene_cache'}},
       'Bio::EnsEMBL::Utils::Cache',
       $SLICE_GENE_CACHE_SIZE);

  return $self;
}



=head2 fetch_all_by_gene_id_list

  Arg [1]    : arrayref $gene_ids  
  Example    : 
  Description: 
  Returntype : listref of Gene objects
  Exceptions : 
  Caller     : 

=cut

sub fetch_all_by_gene_id_list {
  my ($self, $gene_ids, $empty_flag ) = @_;

  my $db = 'core';

  my @genes = ();
  my $core_gene_adaptor = $self->db->get_db_adaptor('core')->get_GeneAdaptor;

  if($empty_flag && @$gene_ids) {
    my $gene_list = join(', ', @$gene_ids); 
    my $sth = $self->prepare("SELECT gene_id, gene_name, chr_name, chr_start, 
                                     chr_end, chr_strand
                              FROM   gene
                              WHERE  db = '$db' AND gene_id IN ($gene_list)");
    
    $sth->execute();
      
    my ($gene_id, $gene_name, $chr_name, $chr_start, $chr_end, $chr_strand);

    $sth->bind_columns(\$gene_id, \$gene_name, \$chr_name, \$chr_start, 
		       \$chr_end, \$chr_strand);

    while($sth->fetch()) {
      my $gene = new Bio::EnsEMBL::Gene();
      $gene->stable_id($gene_name);
      $gene->start($chr_start);
      $gene->end($chr_end);
      $gene->strand($chr_strand);
      $gene->chr_name($chr_name);
      $gene->adaptor($core_gene_adaptor);
      $gene->dbID($gene_id);
      push @genes, $gene;
    }

    return \@genes;
  }

  $self->warn("fetch_by_gene_id_list for non-empty genes not yet implemented");

  return [];
}
  


=head2 fetch_by_gene_id_list

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_gene_id_list instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_gene_id_list {
  my ($self, @args) = @_;

  $self->warn("fetch_by_gene_id_list has been renamed fetch_all_by_gene_id_list\n" . caller);

  return $self->fetch_all_by_gene_id_list(@args);
}




=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice we want genes on
  Arg [2]    : boolean $empty_flag
               A flag as to whether or not empty gene objects should be 
               obtained.  Empty gene objects are light weight and only 
               contain the start, end, source, and name of the gene object
  Function   : retrieve all the genes on this slice. 
               uses www_transcript to get info
  Returntype : listref of Bio::EnsEMBL::Gene objects
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice

=cut

sub fetch_all_by_Slice {
    my ( $self, $slice, $logic_name, $empty_flag ) = @_;

    my $key = $slice->name . (defined ($logic_name) ? ":$logic_name" : "");

    if($empty_flag) {
    # return from cache or the _get_empty_Genes fn while caching results....
        return $self->{'_slice_empty_gene_cache'}{$key} ||= 
                   $self->_get_empty_Genes($slice, $logic_name);
    }

  #check the cache which uses the slice name as it key
    if($self->{'_slice_gene_cache'}{$slice->name()}) {
        return $self->{'_slice_gene_cache'}{$slice->name()};
    }

    my $where = 'WHERE t.chr_name = ? and t.chr_start <= ? and 
                       t.chr_start >= ? and t.chr_end >= ? 
                 AND g.gene_id = t.gene_id
                 AND g.db = t.db';

    my @bind_vals = ( $slice->chr_name, 
		      $slice->chr_end, 
		      $slice->chr_start - $MAX_TRANSCRIPT_LENGTH, 
		      $slice->chr_start);

    if($logic_name) {
      $where .= " and g.analysis = ?";
      push @bind_vals, $logic_name;
    }

    my $sth = $self->prepare(       
      "SELECT t.id, t.transcript_id, t.chr_name, t.chr_start, t.chr_end, 
              t.chr_strand, t.transcript_name, t.translation_id, 
              t.translation_name, t.gene_id, t.type, t.gene_name, t.db, 
              t.exon_structure, t.external_name, t.external_status, 
              t.exon_ids, t.external_db, t.coding_start, t.coding_end, 
              g.chr_start as gene_start, g.chr_end as gene_end,
              g.external_name as gene_external_name, 
              g.external_db as gene_external_db, 
              g.external_status as gene_external_status, 
              g.type as gene_type,
              g.analysis as analysis
        FROM  transcript t, gene g 
        $where");

    $sth->execute( @bind_vals );
 
    return $self->{'_slice_gene_cache'}{$key} =
                $self->_objects_from_sth( $sth, $slice );
}



sub fetch_by_DBEntry {
    my ($self, $db, $dbentry, $chr_coords ) = @_;
    my $sth = $self->prepare( "select gene_name from gene_xref where external_name = ? and db = ?" );
    $sth->execute( $dbentry, $db );
    my( $stable_id ) = $sth->fetchrow;
    return $stable_id ? $self->fetch_by_stable_id( $stable_id, $chr_coords ) : undef;
}

=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
  Arg [2]    : (optional) boolean $chr_coords
               flag set to 1 if genes should be retrieves in chromosomal coords
               by default they are retrieved in contig coords
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG000001432');
  Description: Retrieves a gene via its stable id.  By default genes are 
               retrieved in contig coords, but if the chr_coords flag is set
               then they may be retrieved in chromosomal corods. 
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : geneview, familyview, general

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id, $chr_coords) = @_;
  warn( "~~~~~~~~~~~~~~" );
  my $core_db_adaptor = $self->db->get_db_adaptor('core');

  my $sth = $self->prepare
    ( "SELECT t.id, t.transcript_id, t.chr_name, t.chr_start, t.chr_end, 
              t.chr_strand, t.transcript_name, t.translation_id, 
              t.translation_name, t.gene_id, t.type, t.gene_name, 
              t.db, t.exon_structure, t.external_name, t.external_status, t.exon_ids,
              t.external_db, t.coding_start, t.coding_end, 
              g.chr_start as gene_start, g.chr_end as gene_end,
              g.external_name as gene_external_name, 
              g.external_status as gene_external_status, 
              g.external_db as gene_external_db, g.type as gene_type,
              g.analysis as analysis 
        FROM  transcript t, gene g 
        WHERE g.gene_id = t.gene_id
          AND g.db = t.db 
          AND t.gene_name = ? " );

  my $slice = Bio::EnsEMBL::Slice->new
    ( '-empty' => 1,
      '-adaptor' => $core_db_adaptor->get_SliceAdaptor());
  
  $sth->execute( $stable_id );

  my ( $gene ) = @{$self->_objects_from_sth( $sth, $slice )};

  unless($chr_coords) {
    #transform gene to RawContig coords:
    $gene->transform;
  }

  return $gene;
}


sub fetch_by_transcript_stable_id {
  my ($self, $stable_id ) = @_;
  my $core_db_adaptor = $self->db->get_db_adaptor('core');

  my $sth = $self->prepare
    ( "SELECT t.id, t.transcript_id, t.chr_name, t.chr_start, t.chr_end, 
              t.chr_strand, t.transcript_name, t.translation_id, 
              t.translation_name, t.gene_id, t.type, t.gene_name, t.db, 
              t.exon_structure, t.external_name, t.external_status, t.exon_ids,
              t.external_db, t.coding_start, t.coding_end, 
              g.chr_start as gene_start, g.chr_end as gene_end,
              g.external_name as gene_external_name, 
              g.external_status as gene_external_status, 
              g.external_db as gene_external_db, g.type as gene_type,
              g.analysis as analysis 
        FROM  transcript t, gene g 
        WHERE g.gene_id = t.gene_id
          AND g.db = t.db 
          AND t.transcript_name = ?" );

  my $slice = Bio::EnsEMBL::Slice->new
    ( '-empty' => 1,
      '-adaptor' => $core_db_adaptor->get_SliceAdaptor());
  
  $sth->execute( $stable_id );
  
  my ( $gene ) = @{$self->_objects_from_sth( $sth, $slice )};

  #transform gene to rawcontig coords
  $gene->transform;

  return $gene;
}



sub _objects_from_sth {
  my ( $self, $sth, $slice ) = @_;

  # have to make gene, transcripts, translation, db_link for gene and exons

  my %exon_cache;
  my %gene_cache;
  my %analysis_cache;
  my $core_db_adaptor = $self->db->get_db_adaptor('core');


  my ( $gene, $transcript, $translation ); 
  my ( $exon_id );

  while( my $hr = $sth->fetchrow_hashref() ) {
    my $db = $hr->{'db'};
    unless( $slice && $slice->chr_name ) {
      #retrieve a slice for the entire chromosome
      my $chr = $hr->{'chr_name'};
      $slice = $core_db_adaptor->get_SliceAdaptor->fetch_by_chr_name($chr);
    }
    my $slice_start = $slice->chr_start;

    if( !exists $gene_cache{ $hr->{'db'}."-".$hr->{gene_id} } ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->stable_id( $hr->{'gene_name'} );
      $gene->dbID( $hr->{'gene_id'} );
      $gene->adaptor( $core_db_adaptor->get_GeneAdaptor() );
      #$gene->source( $hr->{'db'} );
      $gene->strand( $hr->{'chr_strand'} );

      #set the gene's analysis
      $analysis_cache{$hr->{'analysis'}} ||= 
	Bio::EnsEMBL::Analysis->new(-logic_name => $hr->{'analysis'}); 
      $gene->analysis($analysis_cache{$hr->{'analysis'}});

      $gene->chr_name( $hr->{'chr_name'} );
      $gene->start( $hr->{'gene_start'} );
      $gene->end( $hr->{'gene_end'} );
      $gene->strand( $hr->{'chr_strand'} );
      if( defined $hr->{'gene_type' } ) {
        $gene->external_status( $hr->{'gene_external_status'} );
	$gene->external_name( $hr->{'gene_external_name'} );
	$gene->external_db( $hr->{'gene_external_db'} );
	$gene->type( $hr->{'gene_type'} );
      }
      $gene_cache{ $hr->{'db'}."-".$hr->{gene_id} } = $gene;
    } else {
      $gene = $gene_cache{ $hr->{'db'}."-".$hr->{gene_id} };
    }

    my $strand = $hr->{'chr_strand'} * $slice->strand;

    # create exons from exon_structure entry
    my @exons = ();
    my @lengths = split( ":", $hr->{'exon_structure'} );
    my @exon_ids = split( ":", $hr->{'exon_ids'} );
    my ( $start, $end );

    my ($start_exon, $end_exon, $transl_coding_start, $transl_coding_end);
    my $coding_start = 0;
    if($hr->{'coding_start'}) { 
      $coding_start = $hr->{'coding_start'} - $slice_start + 1;
    }
    my $coding_end = 0;
    if($hr->{'coding_end'}) {
      $coding_end = $hr->{'coding_end'} - $slice_start + 1;
    }
    
    # lowest chr coord  exon first
    $start = $hr->{'chr_start'} - $slice_start + 1;
    $end = $start + $lengths[0] - 1;
    shift( @lengths );
    $exon_id = shift( @exon_ids );

    my $exon;
    if( ! exists $exon_cache{ "$db:$exon_id" } ) {
      $exon = Bio::EnsEMBL::Exon->new_fast( $slice, $start, $end, $strand);
      #  we need dbIDs for Exons !!!
      #   $exon->dbID( );
      # this is not right for source != core ...
      $exon->adaptor( $core_db_adaptor->get_ExonAdaptor() );
      $exon_cache{"$db:$exon_id"} = $exon;
      $exon->dbID( $exon_id );
      $exon->phase( 0 );
      $exon->end_phase( 0 );      
    } else {
      $exon = $exon_cache{"$db:$exon_id"};
    }

    if($coding_start || $coding_end) { #ignore if both = 0
      #check if this is the start (end on - strand) exon of the translation
      if($coding_start >= $start && $coding_start <= $end) {
	if($strand == 1) {
	  $transl_coding_start = $coding_start - $start + 1;
	  $start_exon = $exon;
	} else {
	  $transl_coding_end = $end - $coding_start + 1;
	  $end_exon = $exon;
	}
      }

      #check if this is the end (start on + strand) exon of the translation
      if($coding_end >= $start && $coding_end <= $end) {
	if($strand == 1) {
	  $transl_coding_end = $coding_end - $start + 1;
	  $end_exon = $exon;
	} else {
	  $transl_coding_start = $end - $coding_end + 1;
	  $start_exon = $exon;
	}
      }
    }
    
    $exon->contig( $slice );
    push( @exons, $exon );

    # now the rest of the exons
    while( @lengths ) {
      my $intron_length = shift( @lengths );
      my $exon_length = shift( @lengths );
      $exon_id = shift( @exon_ids );

      $start = $end + $intron_length + 1;
      $end = $start + $exon_length - 1;

      if( ! exists $exon_cache{ "$db:$exon_id" } ) {
	$exon = Bio::EnsEMBL::Exon->new_fast( $slice, $start, $end, 
					      $hr->{'chr_strand'});
	$exon->adaptor( $core_db_adaptor->get_ExonAdaptor() );
	$exon_cache{"$db:$exon_id"} = $exon;
	$exon->dbID( $exon_id );
	# no phase information stored, 
	# putting something to avoid N padding
	$exon->phase( 0 );
	$exon->end_phase( 0 );
	if( ! $exon_id ) {
	  print STDERR "Exon without dbID: $exon\n";
	}
      } else {
	$exon = $exon_cache{"$db:$exon_id"};
      }

      if($coding_start || $coding_end) {
	#check if this is the start (end on - strand) exon of the translation
	if($coding_start >= $start && $coding_start <= $end) {
	  if($strand == 1) {
	    $transl_coding_start = $coding_start - $start + 1;
	    $start_exon = $exon;
	  } else {
	    $transl_coding_end = $end - $coding_start + 1;
	    $end_exon = $exon;
	  }
	}
	
	#check if this is the end (start on + strand) exon of the translation
	if($coding_end >= $start && $coding_end <= $end) {
	  if($strand == 1) {
	    $transl_coding_end = $coding_end - $start + 1;
	    $end_exon = $exon;
	  } else {
	    $transl_coding_start = $end - $coding_end + 1;
	    $start_exon = $exon;
	  }
	}
      }

      push( @exons, $exon );
    }
    
    # create the transcript 
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->adaptor( $core_db_adaptor->get_TranscriptAdaptor() );
    $transcript->dbID( $hr->{'transcript_id'});
    $transcript->coding_region_start( $coding_start );
    $transcript->coding_region_end( $coding_end );
    $transcript->stable_id( $hr->{ 'transcript_name' });
    $transcript->type( $hr->{ 'type' } );
    $transcript->external_status( $hr->{'external_status'} );
    $transcript->external_name( $hr->{'external_name'} );
    $transcript->external_db( $hr->{'external_db' } );

    # Add the exons
    if( $hr->{'chr_strand'} != 1 ) {
      @exons = reverse( @exons );
    }

    for my $exon ( @exons ) {
      $transcript->add_Exon( $exon );
    }

    # translation ..
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->adaptor( $core_db_adaptor->get_TranslationAdaptor() );
    $translation->stable_id( $hr->{'translation_name'} );
    $translation->dbID( $hr->{'translation_id'} );
    
    #some sanger transcripts don't have proper translations defined
    $start_exon && $translation->start_Exon($start_exon);
    $end_exon && $translation->end_Exon($end_exon);
    if($start_exon && $end_exon) {
      $translation->start($transl_coding_start);
      $translation->end($transl_coding_end);
    }
    
    $transcript->translation( $translation );

    $gene->add_Transcript($transcript);
  }

  my @out = ();
  push @out, values( %gene_cache );
  
#  for my $gene ( @out ) {
#    for my $exon ( $gene->get_all_Exons() ) {
#      if( ! $exon->dbID() ) {
#	print STDERR "Exon $exon has no dbID.\n";
#      }
#    }
#  }

  return \@out;
}


=head2 _get_empty_Genes

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : none
  Description: PRIVATE retrieves a list of 'empty' (partially populated) gene
               objects from the lite database.  Designed for swift retrieval
               for the web.
  Returntype : listref of Bio::EnsEMBL::Gene objects
  Exceptions : none
  Caller     : fetch_by_Slice

=cut

sub _get_empty_Genes {
  my ($self, $slice, $logic_name) = @_;

  my $chr_start = $slice->chr_start();
  my $chr_end = $slice->chr_end();
  my $chr_name = $slice->chr_name();


  my $where = "WHERE  g.chr_name = ? AND g.chr_start <= ? AND 
              g.chr_start >= ? AND g.chr_end >= ?";

  
  my @bind_vals = ( $chr_name, 
		    $chr_end,
		    $chr_start - $MAX_TRANSCRIPT_LENGTH,
		    $chr_start );


  if($logic_name) {
    $where .= " and g.analysis = ?";
    push @bind_vals, $logic_name;
  } 

  my $sth = $self->prepare
    ( "SELECT g.db, g.gene_id, g.gene_name, g.chr_name, g.chr_start, 
              g.chr_end, g.chr_strand, g.type, g.external_name, g.external_db,
              g.external_status, g.analysis
       FROM   gene g
      $where"
    );

  $sth->execute(@bind_vals);
  
  my @out = ();

  my $core_gene_adaptor = $self->db->get_db_adaptor('core')->get_GeneAdaptor;

  my $hashref;

  my %analysis_cache;

  while($hashref = $sth->fetchrow_hashref()) {
    my $gene = new Bio::EnsEMBL::Gene();
    $gene->start($hashref->{'chr_start'} - $chr_start);
    $gene->end($hashref->{'chr_end'} - $chr_start);
    $gene->stable_id( $hashref->{'gene_name'} );
    $gene->dbID( $hashref->{'gene_id'} );
    $gene->adaptor( $core_gene_adaptor );
    #$gene->source( $hashref->{'db'} );
    $gene->strand( $hashref->{'chr_strand'} );

    my $analysis = $analysis_cache{$hashref->{'analysis'}} ||= 
      Bio::EnsEMBL::Analysis->new(-logic_name => $hashref->{'analysis'});

    $gene->analysis($analysis);
				  
    if( defined $hashref->{'type' } ) {
      $gene->external_status( $hashref->{'external_status'} );
      $gene->external_name( $hashref->{'external_name'} );
      $gene->external_db( $hashref->{'external_db'} );
      $gene->type( $hashref->{'type'} );
    }
    push @out, $gene;
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

  #print STDERR "\t\tLite::GeneAdaptor::deleteObj\n";

  #call superclass destructor 
  $self->SUPER::deleteObj;

  #flush the cache
  %{$self->{'_slice_gene_cache'}} = ();
  %{$self->{'_external_gene_cache'}} = ();
}


sub fetch_all_by_external_name {
    my ( $self, $name, $db, $empty_flag ) = @_;

    $db ||= 'core';

    $self->throw('name argument is required') if (!$name);

    if($empty_flag) {
    # return from cache or the _get_empty_Genes fn while caching results....
        return $self->{'_external_empty_gene_cache'}{"$db:$name"} ||= $self->_get_empty_Genes_by_external_name( $name, $db );
    }

  #check the cache which uses the slice name as it key
    if($self->{'_external_gene_cache'}{"$db:$name"}) {
        return $self->{'_external_gene_cache'}{"$db:$name"};
    }

    my $sth = $self->prepare( "SELECT t.id, t.transcript_id, t.chr_name, t.chr_start, t.chr_end,
              t.chr_strand, t.transcript_name, t.translation_id,
              t.translation_name, t.gene_id, t.type, t.gene_name, t.db,
              t.exon_structure, t.external_name, t.exon_ids, t.external_db, t.external_status,
              t.coding_start, t.coding_end,
              g.external_name as gene_external_name,
              g.external_db as gene_external_db, g.external_status as gene_external_status,
              g.type as gene_type,
              g.analysis as analysis
         FROM transcript t, gene_xref as gx, gene as g
        where g.gene_id = t.gene_id AND g.db = t.db and t.db = ? and gx.external_name = ? and gx.gene_id = g.gene_id
        order by g.gene_name, t.gene_name"
    );

    $sth->execute( $db, $name );

    return $self->{'_external_gene_cache'}{"$db:$name"} = $self->_objects_from_sth( $sth, "" );
}

sub count_by_external_name { 
  my ($self, $name, $db) = @_;
  my $sth = $self->prepare( "SELECT count(distinct g.gene_name) FROM   gene g, gene_xref as gx WHERE  g.gene_id = gx.gene_id and g.db = ? and gx.external_name = ?" );
  $sth->execute( $db, $name );
  my ($count) = $sth->fetchrow_array();
  return $count;
}

sub _get_empty_Genes_by_external_name {
  my ($self, $name, $db) = @_;

  my $sth = $self->prepare(
	 "SELECT distinct g.db, g.gene_id, g.gene_name, g.chr_name, g.chr_start,
              g.chr_end, g.chr_strand, g.type, g.external_name, g.external_db, g.external_status
       FROM   gene g, gene_xref as gx
       WHERE  g.gene_id = gx.gene_id and g.db = ? and gx.external_name = ?" );

  $sth->execute( $db, $name );
  my @out = ();
  my $core_gene_adaptor = $self->db->get_db_adaptor('core')->get_GeneAdaptor;

  my $hashref;

  while($hashref = $sth->fetchrow_hashref()) {
    my $gene = new Bio::EnsEMBL::Gene();
    $gene->chr_name($hashref->{'chr_name'});
    $gene->start($hashref->{'chr_start'});
    $gene->end($hashref->{'chr_end'});
    $gene->stable_id( $hashref->{'gene_name'} );
    $gene->dbID( $hashref->{'gene_id'} );
    $gene->adaptor( $core_gene_adaptor );
    #$gene->source( $hashref->{'db'} );
    $gene->strand( $hashref->{'chr_strand'} );

    if( defined $hashref->{'type' } ) {
      $gene->external_status( $hashref->{'external_status'} );
      $gene->external_name( $hashref->{'external_name'} );
      $gene->external_db( $hashref->{'external_db'} );
      $gene->type( $hashref->{'type'} );
    }
    push @out, $gene;
  }

  return \@out;
}


=head2 store

  Arg [1]    : none, string, int, Bio::EnsEMBL::Example $formal_parameter_name
    Additional description lines
    list, listref, hashref
  Example    :  ( optional )
  Description: testable description
  Returntype : none, txt, int, float, Bio::EnsEMBL::Example
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut

sub store {
  my ( $self, $gene ) = @_;

  # first make sure that gene is in chromsomal coordinates
  my $sa = $gene->adaptor()->db()->get_SliceAdaptor;
  my $chr_slice = Bio::EnsEMBL::Slice->new( -empty => 1,
					    -adaptor => $sa );

  $gene->transform( $chr_slice );

  my $sth = $self->prepare
    ( " INSERT INTO gene
        SET
          db = ?,
          analysis = ?,
          type = ?,
          gene_id = ?,
          gene_name = ?,
          chr_name = ?,
          chr_start = ?,
          chr_end = ?,
          chr_strand = ?,
          description = ?,
          external_db = ?,
          external_status = ?,
          external_name = ? " 
    );

  $sth->execute( $gene->source, $gene->analysis->logic_name(),
		 $gene->type(), $gene->dbID(), $gene->stable_id(),
		 $gene->chr_name(), $gene->start(), $gene->end(),
		 $gene->strand(), $gene->description(), $gene->external_db(),
                 $gene->external_status(), 
		 $gene->external_name() );

 # now store transcripts

  $sth = $self->prepare
    ( " INSERT INTO transcript
           SET
               db = ?,
         analysis = ?,
             type = ?,
    transcript_id = ?,
  transcript_name = ?,
         chr_name = ?,
        chr_start = ?,
          chr_end = ?,
       chr_strand = ?,
     coding_start = ?,
       coding_end = ?,
   translation_id = ?,
 translation_name = ?,
          gene_id = ?,
        gene_name = ?,
   exon_structure = ?,
         exon_ids = ?,
     external_db = ?,
  external_status = ?,
    external_name = ? " 
    );

  for my $tr ( @{$gene->get_all_Transcripts() } ) {
    
    my $exon_structure = "";
    my $exon_ids = "";
    my $last = undef;

    my @all_exons = @{$tr->get_all_Exons()};
    if( $gene->strand() ==-1 ) {
      @all_exons = reverse( @all_exons );
    }
    for my $ex ( @all_exons ) {
      if( defined $last ) {
	$exon_structure .= ":";
	$exon_structure .= $ex->start() - $last->end() - 1;
	$exon_structure .= ":";

	$exon_ids .= ",";
      }

      $exon_structure .= $ex->length;
      $exon_ids .= $ex->dbID();

      $last = $ex;
    }

    $sth->execute( $gene->source(), $gene->analysis()->logic_name(),
		   $gene->type, $tr->dbID(), $tr->stable_id, $gene->chr_name,
		   $tr->start(), $tr->end(), $gene->strand(), $tr->coding_region_start(),
		   $tr->coding_region_end(), $tr->translation->dbID(), 
		   $tr->translation->stable_id(), $gene->dbID(), $gene->stable_id(),
		   $exon_structure, $exon_ids, $tr->external_db(), $tr->external_status(),
		   $tr->external_name() );
  }
}




1;
__END__

