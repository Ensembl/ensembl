# EnsEMBL Gene reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 15.07.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneLiteAdaptor - MySQL Database queries to retrieve genes quickly from denormalized tables.

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
use Bio::EnsEMBL::Gene;


use vars '@ISA';


@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

my $MAX_TRANSCRIPT_LENGTH=3000000;


=head2 fetch_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              The slice we want genes on
  Function  : retrieve all the genes on this slice. 
              uses www_transcript to get info
  Returntype: list of Bio::EnsEMBL::Gene
  Exceptions: none
  Caller    : Bio::EnsEMBL::Slice

=cut

sub fetch_by_Slice {
  my ( $self, $slice ) = @_;
  my @out;
  my $core_DBAdaptor = $self->db->core_DBAdaptor();

  my $sth = $self->prepare
    ( "SELECT id, transcript_id, chr_name, chr_start, chr_end, chr_strand,
              transcript_name, translation_id, translation_name, gene_id,
              type, gene_name, db, exon_structure, external_name,
              external_db, coding_start, coding_end 
        FROM  transcript
        WHERE chr_name = ? and chr_start <= ? and chr_start >= ? and
              chr_end >= ?"
    );
    
  eval {
    $sth->execute( $slice->chr_name, $slice->chr_end, $slice->chr_start-$MAX_TRANSCRIPT_LENGTH, $slice->chr_start );
  };

  return () if($@);

  # have to make gene, transcripts, translation, db_link for gene and exons

  my %exon_cache = ();
  my %gene_cache = ();

  my ( $gene, $transcript, $translation ); 

  while( my $hr = $sth->fetchrow_hashref() ) {

    if( !exists $gene_cache{ $hr->{gene_id} } ) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->stable_id( $hr->{'gene_name'} );
      $gene->dbID( $hr->{'gene_id'} );
      $gene->source( $hr->{'db'} );
      $gene->adaptor( $core_DBAdaptor->get_GeneAdaptor() );
      $gene_cache{ $hr->{gene_id} } = $gene;

      print STDERR "Created Gene\n";

    } else {
      $gene = $gene_cache{ $hr->{gene_id} };
    }

    # create exons from exon_structure entry
    my @exons = ();
    my @lengths = split( ":", $hr->{'exon_structure'} );
    my ( $start, $end );

    # lowest chr coord  exon first
    $start = $hr->{'chr_start'} - $slice->chr_start + 1;
    $end = $start + $lengths[0] - 1;
    shift( @lengths );

    my $exon;
    if( ! exists $exon_cache{ "$start-$end" } ) {
      $exon = Bio::EnsEMBL::Exon->new
	( $start, $end, $hr->{'chr_strand'}*$slice->strand());
      #  we need dbIDs for Exons !!!
      #   $exon->dbID( );
      $exon->contig( $slice );
      $exon->adaptor( $core_DBAdaptor->get_ExonAdaptor() );
      $exon_cache{"$start-$end"} = $exon;
      print STDERR "\tCreated Exon\n";
    } else {
      $exon = $exon_cache{"$start-$end"};
    }
    $exon->contig( $slice );
    push( @exons, $exon );

    # now the rest of the exons
    while( @lengths ) {
      my $intron_length = shift( @lengths );
      my $exon_length = shift( @lengths );
      
      $start = $end + $intron_length + 1;
      $end = $start + $exon_length - 1;

      if( ! exists $exon_cache{ "$start-$end" } ) {
	$exon = Bio::EnsEMBL::Exon->new
	  ( $start, $end, $hr->{'chr_strand'});
	$exon->contig( $slice );
	$exon->adaptor( $core_DBAdaptor->get_ExonAdaptor() );
	$exon_cache{"$start-$end"} = $exon;
	print STDERR "\tCreated Exon\n";
      } else {
	$exon = $exon_cache{"$start-$end"};
      }

      push( @exons, $exon );
    }
    
    # create the transcript 
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->adaptor( $core_DBAdaptor->get_TranscriptAdaptor() );
    $transcript->dbID( $hr->{'transcript_id'});
    $transcript->coding_start( $hr->{'coding_start'} );
    $transcript->coding_end( $hr->{'coding_end'} );


    $transcript->stable_id( $hr->{ 'transcript_name' });
    $transcript->type( $hr->{ 'transcript_type' } );
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
    $translation->adaptor( $core_DBAdaptor->get_TranslationAdaptor() );
    $translation->stable_id( $hr->{'translation_name'} );
    $translation->dbID( $hr->{'translation_id'} );
    $transcript->translation( $translation ); 

    $gene->add_Transcript($transcript);

    # we need start and end Exon
    # hope they are lazy loaded ... nope they are not!!!
  }

  my @out = values( %gene_cache );

  return @out;
}


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
   $self->warn( "Use the GeneAdaptor for this query" );
   return ();
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
   $self->warn( "Use the GeneAdaptor for this query" );
   return ();
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
   $self->warn( "Use the GeneAdaptor for this query" );
  return undef;
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
  $self->warn( "Use the GeneAdaptor for this query" );
  return undef;
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
  $self->warn( "Use the GeneAdaptor for this query" );
  return undef;
}


sub fetch_by_Transcript_id {
    my $self = shift;
    $self->warn( "fetch_by_Transcript_id not YET supported from this adaptor" );
    return undef;
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
    $self->warn( "fetch_by_Peptide_id not YET supported from this adaptor" );
    return undef;
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
    $self->warn( "fetch_by_maximum_DBLink not supported by this adaptor" );
    return undef;
}


=head2 get_description

 Title   : get_description
 Usage   : $geneAdptor->get_description($dbID)
 Function: gets gene description line 
 Returns : a string
 Args    : 


=cut

sub get_description {
    my $self = shift;
    $self->warn( "get_description not supported by this adaptor" );
    return undef;
}

=head2 get_stable_entry_info

 Title   : get_stable_entry_info
 Usage   : $geneAdptor->get_stable_entry_info($gene)
 Function: gets stable info for gene and places it into the hash
 Returns : 
 Args    : 


=cut

sub get_stable_entry_info {
    my $self = shift;
    $self->warn( "get_stable_entry_info not supported by this adaptor" );
    return undef;
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
    $self->warn( "fetch_by_DBEntry not supported by this adaptor" );
    return undef;
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
    my $self = shift;
    $self->warn( "store not supported by this adaptor" );
    return undef;
}


sub remove {
    my $self = shift;
    $self->warn( "remove not supported by this adaptor" );
    return undef;
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
    my $self = shift;
    $self->warn( "get_Interpro_by_geneid not supported by this adaptor" );
    return undef;
}


sub create_tables {
# read sql from geneAdaptor.sql
}



1;
__END__

