#
# BioPerl module for DBSQL::DBAdaptorI
#
# Cared for by Frank Visser <fvisser@hgmp.mrc.ac.uk>
#
# Copyright Frank Visser
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::DBAdaptorI - Object representing an instance of
an EnsEMBL DB, independent of the underlying database.

=head1 SYNOPSIS

=head1 DESCRIPTION 

This object represents a database independent
module. All code shared between the different DBAdaptors is stored in
this module. For the moment this module is empty and is only used to
make sure that the isa tests in the Adaptors test for an database
independent module

=cut

package Bio::EnsEMBL::DBSQL::DBAdaptorI;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DB::ObjI;

@ISA = qw(Bio::EnsEMBL::DB::ObjI Bio::EnsEMBL::Root);

sub new {};
1;
