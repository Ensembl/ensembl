
#
# BioPerl module for VersionedSeqAdaptor
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

VersionedSeqAdaptor - DB Adaptor for ArchiveSeq objects

=head1 SYNOPSIS

    my $asad = Bio::EnsEMBL::Archive::DBSQL->new($db);
    my @aseqs = $asad->fetch_by_ensembl_id('ENSE00000023423');

=head1 DESCRIPTION

The VersionedSeqAdaptor contains all the SQL needed to fetch/store ArchiveSeq objects in the Archive Database

=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::Overlap::DBSQL::VersionedSeqAdaptor;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Archive::DBSQL::ArchiveSeq;

use Bio::Root::RootI;

#The method ->new is inherited from the BaseAdaptor
@ISA = qw(Bio::EnsEMBL::Archive::DBSQL::BaseAdaptor);


1;








