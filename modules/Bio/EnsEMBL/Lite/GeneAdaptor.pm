# EnsEMBL Gene reading writing adaptor for mySQL
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
    my ( $self, $slice, $empty_flag ) = @_;

    if($empty_flag) {
    # return from cache or the _get_empty_Genes fn while caching results....
        return $self->{'_slice_empty_gene_cache'}{$slice->name()} ||= 
                   $self->_get_empty_Genes($slice);
    }

  #check the cache which uses the slice name as it key
    if($self->{'_slice_gene_cache'}{$slice->name()}) {
        return $self->{'_slice_gene_cache'}{$slice->name()};
    }

    my $sth = $self->prepare( "SELECT t.id, t.transcript_id, t.chr_name, t.chr_start, t.chr_end, 
              t.chr_strand, t.transcript_name, t.translation_id, 
              t.translation_name, t.gene_id, t.type, t.gene_name, t.db, 
              t.exon_structure, t.external_name, t.exon_ids, t.external_db, 
              t.coding_start, t.coding_end, 
              g.external_name as gene_external_name, 
              g.external_db as gene_external_db, g.type as gene_type 
        FROM  transcript t 
    LEFT JOIN gene g 
           ON g.gene_id = t.gene_id
          AND g.db = t.db 
        WHERE t.chr_name = ? and t.chr_start <= ? and t.chr_start >= ? and
              t.chr_end >= ?"
    );

    $sth->execute( $slice->chr_name, $slice->chr_end, 
		   $slice->chr_start - $MAX_TRANSCRIPT_LENGTH, 
		   $slice->chr_start );
 
    return $self->{'_slice_gene_cache'}{$slice->name} =
                $self->_objects_from_sth( $sth, $slice );

}


=head2 fetch_by_Slice

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Slice instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Slice {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Slice has been renamed fetch_all_by_Slice\n" . caller);

  return $self->fetch_all_by_Slice(@args);
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
  my $core_db_adaptor = $self->db->get_db_adaptor('core');

  my $sth = $self->prepare
    ( "SELECT t.id, t.transcript_id, t.chr_name, t.chr_start, t.chr_end, 
              t.chr_strand, t.transcript_name, t.translation_id, 
              t.translation_name, t.gene_id, t.type, t.gene_name, 
              t.db, t.exon_structure, t.external_name, t.exon_ids,
              t.external_db, t.coding_start, t.coding_end, 
              g.external_name as gene_external_name, 
              g.external_db as gene_external_db, g.type as gene_type 
        FROM  transcript t 
    LEFT JOIN gene g 
           ON g.gene_id = t.gene_id
          AND g.db = t.db 
        WHERE t.gene_name = ? " );

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
              t.exon_structure, t.external_name, t.exon_ids,
              t.external_db, t.coding_start, t.coding_end, 
              g.external_name as gene_external_name, 
              g.external_db as gene_external_db, g.type as gene_type 
        FROM  transcript t 
    LEFT JOIN gene g 
           ON g.gene_id = t.gene_id
          AND g.db = t.db 
        WHERE t.transcript_name = ? " );

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
  my $core_db_adaptor = $self->db->get_db_adaptor('core');


  my ( $gene, $transcript, $translation ); 
  my ( $exon_id );

  while( my $hr = $sth->fetchrow_hashref() ) {
    unless( $slice->chr_name ) {
      #retrieve a slice for the entire chromosome
      my $chr = $hr->{'chr_name'};
      %$slice = %{$core_db_adaptor->get_SliceAdaptor->fetch_by_chr_name($chr)};
    }

    if( !exists $gene_cache{ $hr->{'db'}."-".$hr->{gene_id} } ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->stable_id( $hr->{'gene_name'} );
      $gene->dbID( $hr->{'gene_id'} );
      $gene->adaptor( $core_db_adaptor->get_GeneAdaptor() );
      $gene->source( $hr->{'db'} );
      $gene->strand( $hr->{'chr_strand'} );

      if( defined $hr->{'gene_type' } ) {
	$gene->external_name( $hr->{'gene_external_name'} );
	$gene->external_db( $hr->{'gene_external_db'} );
	$gene->type( $hr->{'gene_type'} );
      }
      $gene_cache{ $hr->{'db'}."-".$hr->{gene_id} } = $gene;
    } else {
      $gene = $gene_cache{ $hr->{'db'}."-".$hr->{gene_id} };
    }

    # create exons from exon_structure entry
    my @exons = ();
    my @lengths = split( ":", $hr->{'exon_structure'} );
    my @exon_ids = split( ":", $hr->{'exon_ids'} );
    my ( $start, $end );

      
    # lowest chr coord  exon first
    $start = $hr->{'chr_start'} - $slice->chr_start + 1;
    $end = $start + $lengths[0] - 1;
    shift( @lengths );
    $exon_id = shift( @exon_ids );

    my $exon;
    if( ! exists $exon_cache{ "$exon_id" } ) {
      $exon = Bio::EnsEMBL::Exon->new_fast( $slice, $start, $end, 
					 $hr->{'chr_strand'}*$slice->strand());
      #  we need dbIDs for Exons !!!
      #   $exon->dbID( );
      # this is not right for source != core ...
      $exon->adaptor( $core_db_adaptor->get_ExonAdaptor() );
      $exon_cache{"$exon_id"} = $exon;
      $exon->dbID( $exon_id );
    } else {
      $exon = $exon_cache{"$exon_id"};
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

      if( ! exists $exon_cache{ "$exon_id" } ) {
	$exon = Bio::EnsEMBL::Exon->new_fast( $slice, $start, $end, 
					      $hr->{'chr_strand'});
	$exon->adaptor( $core_db_adaptor->get_ExonAdaptor() );
	$exon_cache{"$exon_id"} = $exon;
	$exon->dbID( $exon_id );
	# no phase information stored, 
	# putting something to avoid N padding
	$exon->phase( 0 );
	$exon->end_phase( 0 );
	if( ! $exon_id ) {
	  print STDERR "Exon without dbID: $exon\n";
	}
      } else {
	$exon = $exon_cache{"$exon_id"};
      }

      push( @exons, $exon );
    }
    
    # create the transcript 
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->adaptor( $core_db_adaptor->get_TranscriptAdaptor() );
    $transcript->dbID( $hr->{'transcript_id'});
    $transcript->coding_start( $hr->{'coding_start'} -$slice->chr_start() + 1);
    $transcript->coding_end( $hr->{'coding_end'} -$slice->chr_start() + 1);
    $transcript->stable_id( $hr->{ 'transcript_name' });
    $transcript->type( $hr->{ 'type' } );
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
  my ($self, $slice) = @_;

  my $chr_start = $slice->chr_start();
  my $chr_end = $slice->chr_end();
  my $chr_name = $slice->chr_name();
  
  my $sth = $self->prepare
    ( "SELECT g.db, g.gene_id, g.gene_name, g.chr_name, g.chr_start, 
              g.chr_end, g.chr_strand, g.type, g.external_name, g.external_db
       FROM   gene g
       WHERE  g.chr_name = ? AND g.chr_start <= ? AND 
              g.chr_start >= ? AND g.chr_end >= ?" );

  $sth->execute( $chr_name, $chr_end,
		 $chr_start - $MAX_TRANSCRIPT_LENGTH,
		 $chr_start );
  
  my @out = ();

  my $core_gene_adaptor = $self->db->get_db_adaptor('core')->get_GeneAdaptor;

  my $hashref;

  while($hashref = $sth->fetchrow_hashref()) {
    my $gene = new Bio::EnsEMBL::Gene();
    $gene->start($hashref->{'chr_start'} - $chr_start);
    $gene->end($hashref->{'chr_end'} - $chr_start);
    $gene->stable_id( $hashref->{'gene_name'} );
    $gene->dbID( $hashref->{'gene_id'} );
    $gene->adaptor( $core_gene_adaptor );
    $gene->source( $hashref->{'db'} );
    $gene->strand( $hashref->{'chr_strand'} );

    if( defined $hashref->{'type' } ) {
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
}


1;
__END__

