# EnsEMBL Gene Adaptor Interface
#
# Copyright EMBL-EBI 2002
#
# Author: Graham McVicker
# 
# Date : 06.08.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::GeneAdaptorI

=head1 SYNOPSIS

The interface definition which all Gene adaptors should implement

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::GeneAdaptorI;

=head2 fetch_by_Slice

  Arg  1    : Bio::EnsEMBL::Slice $slice
              The slice we want genes on
  Function  : retrieve all the genes on this slice. 
  Returntype: list of Bio::EnsEMBL::Gene
  Exceptions: none
  Caller    : Bio::EnsEMBL::Slice

=cut

sub fetch_by_Slice {}

=head2 list_geneIds

 Title   : list_geneIds
 Usage   : $geneAdaptor->list_geneIds
 Function: Gets an array of internal ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_geneIds {}

=head2 list_stable_geneIds

 Title   : list_stable_geneIds
 Usage   : $geneAdaptor->list_stable_geneIds
 Function: Gets an array of stable ids for all genes in the current db
 Example : 
 Returns : array of ids
 Args    : none

=cut

sub list_stable_geneIds {}


=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $geneobj->fetch_by_dbID( $geneid)
 Function: gets one gene out of the db
 Example : $obj->get($dbID)
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag

=cut


sub fetch_by_dbID {}


=head2 fetch_by_stable_id

 Title   : fetch_by_stable_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_stable_id{}

=head2 fetch_by_contig_list

 Title   : fetch_by_contig_list
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_contig_list {}


sub fetch_by_Transcript_id {}


=head2 fetch_by_Peptide_id 

 Title   : fetch_by_Peptide_id
 Usage   : $geneAdaptor->fetch_by_Peptide_id($peptideid)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : peptide id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_Peptide_id {}


=head2 fetch_by_maximum_DBLink

 Title   : fetch_by_maximum_DBLink
 Usage   : $geneAdaptor->fetch_by_maximum_DBLink($ext_id)
 Function: gets one gene out of the db with 
 Returns : gene object (with transcripts, exons)
 Args    : 


=cut

sub fetch_by_maximum_DBLink {}


=head2 get_description

 Title   : get_description
 Usage   : $geneAdptor->get_description($dbID)
 Function: gets gene description line 
 Returns : a string
 Args    : 


=cut

sub get_description {}

=head2 get_stable_entry_info

 Title   : get_stable_entry_info
 Usage   : $geneAdptor->get_stable_entry_info($gene)
 Function: gets stable info for gene and places it into the hash
 Returns : 
 Args    : 


=cut

sub get_stable_entry_info {}

=head2 fetch_by_DBEntry

 Title   : fetch_by_DBLink
 Usage   : $geneAdptor->fetch_by_DBLink($ext_id)
 Function: gets one gene out of the db with or without supporting evidence
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified,
assumes without
           Note that it is much faster to get genes without supp.evidence!


=cut

sub fetch_by_DBEntry {}




=head2 store

 Title   : store
 Usage   : $geneAdaptor->store($gene)
 Function: writes a particular gene into the database. Assumes that everything 
           has dbIDs ....
 Example :
 Returns : nothing
 Args    : $gene object


=cut

sub store {}


sub remove {}


=head2 get_Interpro_by_geneid

 Title   : get_Interpro_by_geneid
 Usage   : @interproid = $geneAdaptor->get_Interpro_by_geneid($gene->id);
 Function: gets interpro accession numbers by geneid. A hack really -
           we should have a much more structured system than this
 Example :
 Returns : 
 Args    :


=cut

sub get_Interpro_by_geneid {}


sub create_tables {}



1;
__END__

