
#
# BioPerl module for Bio::EnsEMBL::EMBL_Dump.pm
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::EMBL_Dump.pm - A collection of variables, subroutines and modifers for EMBL dumping

=head1 SYNOPSIS

    use Bio::EnsEMBL::EMBL_Dump;

    # adds standard comments, feature table source lines etc for embl
    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($aseq); 
    
    # adds standard dumping information to the aseqstream to drive off the
    # the annseq objects
    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($aseqstream);

    #
    # Complete dumping from a clone object
    #

    my $as = $clone->get_AnnSeq();
    &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($as);
    
    my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);
    &Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($emblout);    
    
    $emblout->write_annseq($as);


=head1 DESCRIPTION

We need this information in a module because more than one script uses it. However,
the information here is really just comments etc to add onto or manipulate embl
dumping. There is no "object" here - it is just a placeholder for all this information.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::EMBL_Dump;
use vars qw( @ISA @EXPORT_OK );
use Exporter;
use strict;
use Carp;
use Bio::Species;

@ISA = ('Exporter');

# Some of these subroutines are used
# by Bio::EnsEMBL::EMBL_Dump_Sanger
@EXPORT_OK = qw(
                add_ensembl_comments
                ensembl_annseq_output
                id_EnsEMBL
                kw_EnsEMBL
                sv_EnsEMBL
                ac_EnsEMBL
                sort_FTHelper_EnsEMBL
                );



=head2 add_ensembl_comments

 Title   : add_ensembl_comments
 Usage   : &Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($aseq);
 Function: adds ensembl comments for
             - de line - reannotated clone
             - comment about ensembl
             - comment interpreting ensembl flat files
 Example :
 Returns : 
 Args    :


=cut

sub add_ensembl_comments {
   my ($aseq) = @_;

   if( !$aseq  ) {
       confess "Ok. No aseq passed into EMBL_Dump and I can't even throw a nice exception!";
   }

   if( !$aseq->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $aseq->throw("not got a EnsEMBL annseq but a $aseq. Not going to add comments");
   }

   $aseq->seq->desc("Reannotated sequence via Ensembl");
   my $comment = Bio::Annotation::Comment->new();
   
   $comment->text("This sequence was reannotated via the Ensembl system. Please visit the Ensembl web site, http://www.ensembl.org/ for more information.");
   $aseq->annotation->add_Comment($comment);

   $comment = Bio::Annotation::Comment->new();
   $comment->text("The reference, comment, description and feature table of the original entry can be found in the DDBJ/EMBL/GenBank database with the identical accession number.");
   $aseq->annotation->add_Comment($comment);
   
   $comment = Bio::Annotation::Comment->new();
   $comment->text("The /gene_id indicates a unique id for a gene, /cds_id a unique id for a translation and a /exon_id a unique id for an exon. These ids are maintained wherever possible between versions. For more information on how to interpret the feature table, please visit http://www.ensembl.org/docs/embl.html.");
   $aseq->annotation->add_Comment($comment);

   $comment = Bio::Annotation::Comment->new();
   $comment->text("All the exons and transcripts in Ensembl are confirmed by similarity to either protein or cDNA sequences.");
   $aseq->annotation->add_Comment($comment); 

   $comment = Bio::Annotation::Comment->new();
   $comment->text("In unfinished, rough draft DNA sequence gene structures can cross fragments and, in these cases, the order and orientation of the fragments is likely to be different from the order in the the International Nucleotide Sequence Databases DDBJ/EMBL/GenBank.");

   $aseq->annotation->add_Comment($comment); 
  
   my $sf = Bio::SeqFeature::Generic->new();
   $sf->start(1);
   $sf->end($aseq->seq->seq_len());
   $sf->strand(1);
   $sf->primary_tag('source');
   $sf->add_tag_value('organism','Homo sapiens');
   $aseq->add_SeqFeature($sf);
   
   my $species = new Bio::Species;
   $species->common_name("Human");
   $species->classification(qw( sapiens Homo Hominidae
				Catarrhini Primates Eutheria
				Mammalia Vertebrata Chordata
				Metazoa Eukaryota ));

   $aseq->species($species);
   # done!
}

=head2 ensembl_annseq_output

 Title   : ensembl_annseq_output
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub ensembl_annseq_output {
   my ($aseqstream) = @_;

   if( !$aseqstream  ) {
       confess "Ok. No aseq passed into EMBL_Dump and I can't even throw a nice exception!";
   }

   if( !$aseqstream->isa('Bio::SeqIO') ) {
       $aseqstream->throw("not got EMBL IO but a $aseqstream. Not going to add output functions");
   }

   if( $aseqstream->can('_post_sort') ) {
       $aseqstream->_post_sort(\&sort_FTHelper_EnsEMBL);
   }

   # attach ensembl specific dumping functions
   $aseqstream->_id_generation_func(\&id_EnsEMBL);
   $aseqstream->_kw_generation_func(\&kw_EnsEMBL);
   $aseqstream->_sv_generation_func(\&sv_EnsEMBL);
   $aseqstream->_ac_generation_func(\&ac_EnsEMBL);
   
}

#########################
# sub routines
#########################

sub id_EnsEMBL {
    my $annseq = shift;

    return $annseq->id;

    # JGRG - is this correct?  I thought phase 3 was HUM.
    my $division = $annseq->htg_phase == 4 ? 'HUM' : 'HTG';
    my $length = $annseq->seq->seq_len();
    my $id = $annseq->embl_id();

    return sprintf("%-9s  ENSEMBL; DNA; %s; %d BP.", $id, $division, $length );
}


sub kw_EnsEMBL {
   my ($annseq) = @_;

   if( !$annseq->can('htg_phase') ) {
       return "ENSEMBL";
   }

   if( $annseq->htg_phase == 4 ) {
       return "HTG.";
   }

   return "HTG; HTGS_PHASE" . $annseq->htg_phase() . ".";
}

sub sv_EnsEMBL {
   my ($annseq) = @_;

   if( ! $annseq->can('sv') || ! $annseq->sv ) {
       return "NO_SV_NUMBER";
   }
   if( $annseq->sv == -1 ) {
       return undef;
   }

   return $annseq->seq->id() . "." . $annseq->sv
}

sub ac_EnsEMBL {
   my ($annseq) = @_;

   return $annseq->id() . ";";
}


BEGIN {
    # A value of 0 is illegal in %sort_order
    my %sort_order = (
        source          => 1,
        CDS             => 2,
        exon            => 3,
        repeat_region   => 4,
    );

    # $last is one more than the largest value in %sort_order
    my $last = (sort {$b <=> $a} values %sort_order)[0] + 1;

    sub sort_FTHelper_EnsEMBL {
        my $a = shift;
        my $b = shift;

        my $a_ord = $sort_order{$a->key} || $last;
        my $b_ord = $sort_order{$b->key} || $last;

        my ($a_5prime) = $a->loc =~ /(\d+)/;
        my ($b_5prime) = $b->loc =~ /(\d+)/;

        # Features are sorted by location if they don't
        # sort by their keys.
        return $a_ord <=> $b_ord  || $a_5prime <=> $b_5prime;
    }
}


1;
