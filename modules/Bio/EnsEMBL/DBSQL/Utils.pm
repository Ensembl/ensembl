# Ensembl module for Bio::EnsEMBL::DBSQL::Utils
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code



=head1 NAME

Bio::EnsEMBL::DBSQL::Utils - Module having the fset2transcript subroutines

=head1 SYNOPSIS

    use Bio::EnsEMBL::DBSQL::Utils;

    &Bio::EnsEMBL::DBSQL::Utils::fset2transcript($fset_id);

=head1 DESCRIPTION

Module containing the fset2transcript* subroutines, which create
transcripts from features. Note that this is now deprecated - use
the routines in Bio::EnsEMBL::TranscriptFactory instead.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut



package Bio::EnsEMBL::DBSQL::Utils;

use strict;
use Bio::EnsEMBL::TranscriptFactory;



sub fset2transcript {
    my ($genscan,$contig)=@_;

    warn("fset2transcript subroutines have been moved: please use Bio::EnsEMBL::TranscriptFactory::fset2transcript() instead");
    
    return &Bio::EnsEMBL::TranscriptFactory::fset2transcript($genscan, $contig);
}

sub fset2transcript_guess_phases {
    my ($fset,$contig) = @_;

    warn("fset2transcript subroutines have been moved: please use Bio::EnsEMBL::TranscriptFactory::fset2transcript_guess_phases() instead");
    
    return &Bio::EnsEMBL::TranscriptFactory::fset2transcript_guess_phases($fset, $contig);
}


sub fset2transcript_3frame {
  my ($fset,$contig) = @_;

    warn("fset2transcript subroutines have been moved: please use Bio::EnsEMBL::TranscriptFactory::fset2transcript_3frame() instead");
    
    return &Bio::EnsEMBL::TranscriptFactory::fset2transcript_3frame($fset, $contig);
}


1;
