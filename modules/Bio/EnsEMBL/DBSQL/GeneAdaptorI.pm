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


sub fetch_by_Slice {}

#sub list_geneIds {}

#sub list_stable_geneIds {}

sub fetch_by_dbID {}

sub fetch_by_stable_id{}

sub fetch_by_contig_list {}

#sub fetch_by_Transcript_id {}

#sub fetch_by_Peptide_id {}

#sub fetch_by_maximum_DBLink {}

sub get_stable_entry_info {}

sub fetch_by_DBEntry {}

sub store {}

sub remove {}

#sub get_Interpro_by_geneid {}




1;
__END__

