=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::SlicingEventAdaptor - Database adaptor for the retrieval and
storage of SplicingEvent objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
  );

  $se_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "SplicingEvent" );

  $se = $se_adaptor->fetch_by_dbID(12);

  @ses = $se_adaptor->fetch_by_stable_id('ENSG00000184129');

  @sess = @{ $se_adaptor->fetch_by__name('...-CNE-1') };

  $slice_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );

  $slice =
    $slice_adaptor->fetch_by_region( 'chromosome', '1', 1, 1000000 );

  @ase = @{ $se_adaptor->fetch_all_by_Slice($slice) };

=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage ofSlicingEvents
objects.

=head1 METHODS

=cut
package Bio::EnsEMBL::DBSQL::SplicingEventAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SplicingEvent;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 list_dbIDs

  Example    : @gene_ids = @{$gene_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all genes in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_dbIDs {
  my ($self,$ordered) = @_;

  return $self->_list_dbIDs("splicing_event", undef, $ordered);
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on.
  Arg [2]    : type of Transcript event
  Arg [3]    : (optional) boolean $load_features
               if true, transcript will be loaded immediately rather than
               lazy loaded later.
=cut

sub fetch_all_by_Slice{
  my $self  = shift;
  my $slice = shift;
  my $type  = shift;
  my $load_features = shift;
  
  my $constraint = "";

  if(defined($type)){
    $constraint .= " and se.type = '$type'";
  }

  my $tes = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);


  # Is there any use in having a splice event without the paris and features??


  if(!$load_features || @$tes < 2) {
    return $tes;
  }

  ## do someother stuff..

  foreach my $te (@$tes){
    $te->get_all_Features();
    $te->get_all_Pairs();
  }

  return $tes;
}
 
 



sub fetch_all_by_Gene{
  my $self = shift;
  my $gene = shift;
  
  my ($splicing_event_id,$seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $name, $gene_id, $type);
  
  my $sth = $self->dbc->prepare("select se.splicing_event_id , se.seq_region_id, se.seq_region_start, se.seq_region_end, se.seq_region_strand, se.name, se.type  from splicing_event se where se.gene_id = ".$gene->dbID);
  
  $sth->execute();
  
  $sth->bind_columns(\$splicing_event_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$name, \$type);   
  
  my @splicing_events;

  my $sa = $self->db->get_SliceAdaptor();

  while($sth->fetch){
    
    my $slice = $sa->fetch_by_seq_region_id($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand);
    
    push( @splicing_events,
	  $self->_create_feature_fast( 'Bio::EnsEMBL::SplicingEvent', {
								       'start'     => $seq_region_start,
								       'end'       => $seq_region_end,
								       'strand'    => $seq_region_strand,
								       'adaptor'   => $self,
								       'slice'     => $slice,
								       'dbID'      => $splicing_event_id,
								       'name'      => $name,
								       'gene_id'   => $gene->dbID,
								       'type'      => $type
								      } ) );
  }
  $sth->finish;

  foreach my $te (@splicing_events){
    $te->get_all_Features();
    $te->get_all_Pairs();
  }
  return \@splicing_events;
}

sub fetch_all_by_Exon{
  my $self = shift;
  my $exon = shift;
  
  my $list = "(";
  my $sth= $self->dbc->prepare("Select distinct(splicing_event_id)  from splicing_event_feature where exon_id = ".$exon->dbID);
  $sth->execute;
  my ($se_id);
  $sth->bind_columns(\$se_id);
  while($sth->fetch){
    $list .= $se_id." ,";
  }
  $sth->finish;
  chop $list; # get rid of last ","
  $list .= ")";

  my ($splicing_event_id,$seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $name, $gene_id, $type);
  
  $sth = $self->dbc->prepare("select se.splicing_event_id , se.seq_region_id, se.seq_region_start, se.seq_region_end, se.seq_region_strand, se.name, se.type, se.gene_id from splicing_event se where se.splicing_event_id in ".$list);
  
  $sth->execute();
  
  $sth->bind_columns(\$splicing_event_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$name, \$type, \$gene_id);   
  
  my @splicing_events;

  my $sa = $self->db->get_SliceAdaptor();

  while($sth->fetch){
    
    my $slice = $sa->fetch_by_seq_region_id($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand);
    
    push( @splicing_events,
	  $self->_create_feature_fast( 'Bio::EnsEMBL::SplicingEvent', {
								       'start'     => $seq_region_start,
								       'end'       => $seq_region_end,
								       'strand'    => $seq_region_strand,
								       'adaptor'   => $self,
								       'slice'     => $slice,
								       'dbID'      => $splicing_event_id,
								       'name'      => $name,
								       'gene_id'   => $gene_id,
								       'type'      => $type
								      } ) );
  }
  $sth->finish;

  foreach my $te (@splicing_events){
    $te->get_all_Features();
    $te->get_all_Pairs();
  }
  return \@splicing_events;
}

sub fetch_all_by_Transcript{
  my $self = shift;
  my $tran = shift;
  
  my $list = "(";
  my $sth= $self->dbc->prepare("Select distinct(splicing_event_id)  from splicing_event_feature where transcript_id = ".$tran->dbID);
  $sth->execute;
  my ($se_id);
  $sth->bind_columns(\$se_id);
  while($sth->fetch){
    $list .= $se_id." ,";
  }
  $sth->finish;
  chop $list; # get rid of last ","
  $list .= ")";

  my ($splicing_event_id,$seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $name, $gene_id, $type);
  
  $sth = $self->dbc->prepare("select se.splicing_event_id , se.seq_region_id, se.seq_region_start, se.seq_region_end, se.seq_region_strand, se.name, se.type, se.gene_id from splicing_event se where se.splicing_event_id in ".$list);
  
  $sth->execute();
  
  $sth->bind_columns(\$splicing_event_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$name, \$type, \$gene_id);   
  
  my @splicing_events;

  my $sa = $self->db->get_SliceAdaptor();

  while($sth->fetch){
    
    my $slice = $sa->fetch_by_seq_region_id($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand);
    
    push( @splicing_events,
	  $self->_create_feature_fast( 'Bio::EnsEMBL::SplicingEvent', {
								       'start'     => $seq_region_start,
								       'end'       => $seq_region_end,
								       'strand'    => $seq_region_strand,
								       'adaptor'   => $self,
								       'slice'     => $slice,
								       'dbID'      => $splicing_event_id,
								       'name'      => $name,
								       'gene_id'   => $gene_id,
								       'type'      => $type
								      } ) );
  }
  $sth->finish;

  foreach my $te (@splicing_events){
    $te->get_all_Features();
    $te->get_all_Pairs();
  }
  return \@splicing_events;
}

# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _tables {
  my $self = shift;

  return ([ 'splicing_event', 'se' ]);
}

# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _columns {
  my $self = shift;

#  my $created_date = $self->db->dbc->from_date_to_seconds("gsi.created_date");
#  my $modified_date = $self->db->dbc->from_date_to_seconds("gsi.modified_date");

  return ( 'se.splicing_event_id', 'se.seq_region_id', 'se.seq_region_start', 'se.seq_region_end', 'se.seq_region_strand', 'se.name', 'se.gene_id', 'se.type' );

}

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;
  
  my ($splicing_event_id,$seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $name, $gene_id, $type);
  
  
  $sth->bind_columns(\$splicing_event_id, \$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$name, \$gene_id, \$type);

  my $sa = $self->db()->get_SliceAdaptor();


  my @splicing_events;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

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
  my $dest_slice_cs;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  my $asma;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_cs     = $dest_slice->coord_system();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id = $dest_slice->get_seq_region_id();
    $asma = $self->db->get_AssemblyMapperAdaptor();
  }

  FEATURE: while($sth->fetch()) {
    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);

    my $slice = $slice_hash{"ID:".$seq_region_id};
    my $dest_mapper = $mapper;

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #obtain a mapper if none was defined, but a dest_seq_region was
    if(!$dest_mapper && $dest_slice && 
       !$dest_slice_cs->equals($slice->coord_system)) {
      $dest_mapper = $asma->fetch_by_CoordSystems($dest_slice_cs,
                                                 $slice->coord_system);
      $asm_cs = $dest_mapper->assembled_CoordSystem();
      $cmp_cs = $dest_mapper->component_CoordSystem();
      $asm_cs_name = $asm_cs->name();
      $asm_cs_vers = $asm_cs->version();
      $cmp_cs_name = $cmp_cs->name();
      $cmp_cs_vers = $cmp_cs->version();
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};
    #
    # remap the feature coordinates to another coord system 
    # if a mapper was provided
    #
    if($dest_mapper) {

      ($seq_region_id,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $dest_mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                              $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($seq_region_id));

      #get a slice in the coord system we just mapped to
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
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

      #throw away features off the end of the requested slice
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	 ( $dest_slice_sr_id != $seq_region_id )) {
	next FEATURE;
      }

      $slice = $dest_slice;
    }

    # Finally, create the new splicing_event.
    push( @splicing_events,
          $self->_create_feature_fast( 'Bio::EnsEMBL::SplicingEvent', {
                                    'start'     => $seq_region_start,
                                    'end'       => $seq_region_end,
                                    'strand'    => $seq_region_strand,
                                    'adaptor'   => $self,
                                    'slice'     => $slice,
                                    'dbID'      => $splicing_event_id,
                                    'name'      => $name,
                                    'gene_id'   => $gene_id,
                                    'type'      => $type
                                  } ) );

  }

  return \@splicing_events;
}


1;


