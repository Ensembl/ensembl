# EnsEMBL Gene Adaptor Interface
#
# Copyright EMBL-EBI 2002
#
# Author: Graham McVicker
# 
# Date : 06.08.2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI

=head1 SYNOPSIS

The interface definition which all RepeatFeature adaptors should implement

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk
  Graham McVicker : mcvicker@ebi.ac.uk

=head1 APPENDIX

=cut

use strict;

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI;

sub fetch_by_Slice {}

sub fetch_by_RawContig {}

sub fetch_by_contig_id {}

sub fetch_by_dbID {}

sub fetch_by_assembly_location {}

sub fetch_by_assembly_location_constraint {}

sub store {}

1;

__END__

