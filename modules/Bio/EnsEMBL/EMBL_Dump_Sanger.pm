
#
# BioPerl module for Bio::EnsEMBL::EMBL_Dump.pm
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code


=pod

=head1 NAME

Bio::EnsEMBL::EMBL_Dump_Sanger - A collection of variables, subroutines and modifers for EMBL dumping

=head1 SYNOPSIS

    use Bio::EnsEMBL::EMBL_Dump_Sanger;

    # adds standard comments, feature table source lines etc for embl
    &Bio::EnsEMBL::EMBL_Dump_Sanger::add_ensembl_comments($aseq); 
    
    # adds standard dumping information to the aseqstream to drive off the
    # the annseq objects
    &Bio::EnsEMBL::EMBL_Dump_Sanger::ensembl_annseq_output($aseqstream);

    #
    # Complete dumping from a clone object
    #

    my $as = $clone->get_AnnSeq();
    &Bio::EnsEMBL::EMBL_Dump_Sanger::add_ensembl_comments($as);
    
    my $emblout = Bio::AnnSeqIO->new( -format => 'EMBL', -fh => \*STDOUT);
    &Bio::EnsEMBL::EMBL_Dump_Sanger::ensembl_annseq_output($emblout);    
    
    $emblout->write_annseq($as);


=head1 DESCRIPTION

Based on B<Bio::EnsEMBL::EMBL_Dump> with comments
and special lines for submitting Sagner Centre
sequence to EMBL.

=head1 CONTACT

James Gilbert B<email> jgrg@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::EMBL_Dump_Sanger;

use strict;
use vars qw( @ISA @EXPORT_OK );
use Exporter;
use Carp;
use Bio::Annotation::Reference;
use Bio::Annotation::DBLink;

# We inherit a number of functions from Bio::EnsEMBL::EMBL_Dump
use Bio::EnsEMBL::EMBL_Dump qw(
                               id_EnsEMBL
                               sort_FTHelper_EnsEMBL
                               );

@ISA = ('Exporter');
@EXPORT_OK = qw(
                add_ensembl_de_cc
                add_source_feature
                add_contig_comments
                add_accession_info
                add_reference
                ensembl_annseq_output
                add_phase_keywords
                add_phase_info
                add_db_link
                );

=head2 add_ensembl_comments

 Title   : add_ensembl_comments
 Usage   : &Bio::EnsEMBL::EMBL_Dump_Sanger::add_ensembl_comments($aseq);
 Function: adds ensembl comments for
             - de line - reannotated clone
             - comment about ensembl
             - comment interpreting ensembl flat files
 Example :
 Returns : 
 Args    :


=cut

BEGIN { my @sang_comments = (

"IMPORTANT: This sequence is unfinished and does not necessarily
represent the correct sequence.  Work on the sequence is in progress and
the release of this data is based on the understanding that the sequence
may change as work continues.  The sequence may be contaminated with
foreign sequence from E.coli, yeast, vector, phage etc.",

"This sequence was annotated via the Ensembl system.  Please visit the
Ensembl web site, http://ensembl.ebi.ac.uk for more information.",

"All the exons and transcripts in Ensembl are confirmed by similarity to
either protein or cDNA sequences.",

    );
    
    # Replace newlines with spaces
    foreach (@sang_comments) { s/\n/ /g };

    sub add_ensembl_de_cc {
        my ($aseq, $ext_clone) = @_;

        # Sanity checks
        unless ($aseq) {
            confess("Error: No \$aseq passed into EMBL_Dump_Sanger so I can't even throw a nice exception!");
        }
        unless ($aseq->isa('Bio::EnsEMBL::AnnSeq')) {
            $aseq->throw("not got a EnsEMBL annseq but a $aseq. Not going to add comments");
        }

        # DE line
        confess("external clone not specified") unless $ext_clone;
        $aseq->seq->desc("Human DNA sequence *** SEQUENCING IN PROGRESS *** from clone $ext_clone");

        # CC lines
        foreach my $text (@sang_comments) {
            my $comment = Bio::Annotation::Comment->new();
            $comment->text($text);
            $aseq->annotation->add_Comment($comment);
        }
    }
}

sub add_source_feature {
    my( $aseq, $chr, $map, $lib, $cln ) = @_;
   
    my $sf = Bio::SeqFeature::Generic->new();
    $sf->start(1);
    $sf->end($aseq->seq->seq_len());
    $sf->strand(1);
    $sf->primary_tag('source');
    
    $sf->add_tag_value('organism',   'Homo sapiens');
    $sf->add_tag_value('chromosome', $chr ) if $chr;
    $sf->add_tag_value('map',        $map ) if $map;
    $sf->add_tag_value('clone',      $cln ) if $cln;
    $sf->add_tag_value('clone_lib',  $lib ) if $lib;
    
    $aseq->add_SeqFeature($sf);
}

=head2 add_contig_comments

 Title   : add_contig_comments 
 Usage   : add_contig_comments( ANNSEQ, CLONE )
 Function: Adds comment lines listing the contigs in the
           sequence in the order in which they appear.
 Example :
 Returns : 
 Args    :


=cut

sub add_contig_comments {
    my( $aseq, $clone ) = @_;

    # Order and ID of the contigs
    my( @cc ) = (

"The order of contigs shown below is Ensembl's best guess.  Contigs which
have been flipped relative to the original data are marked 'Reversed'. 
800 N's separate the contigs:\n");

    foreach my $contig ($clone->get_all_Contigs) {
        my $id          = $contig->id;
        my $length      = $contig->length;
        my $orientation = $contig->orientation;

        $id =~ s/^([^\.]+)\.//
            or confess "Can't remove clone id from '$id'";

        my $text = "Contig_ID: $id  Length: ${length}bp";
        $text .= "  Reversed" if $orientation == -1;
        push( @cc, "$text\n" );
    }
    my $comment = Bio::Annotation::Comment->new();
    $comment->text(join '', @cc);
    $aseq->annotation->add_Comment($comment);

}


=head2 add_accession_info

 Usage   : add_accession_info( ANNSEQ, ACC, ID );
 Function: Adds EMBL accession and entryname information
           to an AnnSeq object in the right places for
           EMBL dumping.
 Example : add_accession_info( $aseq, 'AC000345', 'HS45D13' );
 Returns : 
 Args    : ANNSEQ - a Bio::EnsEMBL::AnnSeq object
           ACC    - Primary EMBL accession number
           ID     - EMBL entryname

=cut

sub add_accession_info {
    my( $aseq, $acc, $id, $project ) = @_;
    
    $aseq->seq->id($acc);
    $aseq->embl_id($id);
    $aseq->project_name($project);
}

sub add_phase_info {
    my( $aseq, $is_finished ) = @_;
    
    if ($is_finished) {
        $aseq->htg_phase(4);
    } else {
        $aseq->htg_phase(1);
    }
}

sub add_db_link {
    my( $aseq, $db, $clone_id ) = @_;
    
    my $dbl = Bio::Annotation::DBLink->new;
    $dbl->database($db);
    # FIXME - HUMAN is hard coded
    $dbl->primary_id("HUMAN-Clone-$clone_id");
    
    $aseq->annotation->add_DBLink($dbl);
}

sub add_reference {
    my( $aseq, $author, $date ) = @_;
    
    confess("No author supplied") unless $author;
    my $reference = Bio::Annotation::Reference->new();
    $reference->authors("$author;");
    $reference->location(
"Submitted ($date) to the EMBL/Genbank/DDBJ databases.
Sanger Centre, Hinxton, Cambridgeshire, CB10 1SA, UK.
E-mail enquiries: humquery\@sanger.ac.uk
Clone requests: clonerequest\@sanger.ac.uk\n");
    
    $aseq->annotation->add_Reference($reference);
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

   unless ($aseqstream) {
       confess("Error: No \$aseqstream passed into EMBL_Dump_Sanger so I can't even throw a nice exception!");
   }

   unless ($aseqstream->isa('Bio::AnnSeqIO::EMBL')) {
       $aseqstream->throw("not got EMBL IO but a $aseqstream. Not going to add output functions");
   }

   $aseqstream->_post_sort(\&sort_FTHelper_EnsEMBL);
   
   # attach ensembl specific dumping functions
   $aseqstream->_id_generation_func(\&id_EnsEMBL);
   $aseqstream->_ac_generation_func(\&ac_EnsEMBL);
   $aseqstream->_sv_generation_func(\&sv_EnsEMBL);
}

sub ac_EnsEMBL {
    my ($annseq) = @_;

    my $acc = $annseq->seq->id()
        or $annseq->throw("No accession in \$annseq");
    my $proj = $annseq->project_name()
        or $annseq->throw("No project_name in \$annseq");
    $proj = '_'. uc $proj;
    
    return "$acc;\nXX\nAC * $proj";
}

# We don't add SV lines to submissions
sub sv_EnsEMBL {  }


sub add_phase_keywords {
    my ($annseq, $draft) = @_;

    my @kw = ('HTG');

    unless ($annseq->htg_phase == 4) {
        push(@kw, 'HTGS_PHASE' . $annseq->htg_phase);
        push(@kw, $draft) if $draft;
    }
    
    $annseq->keywords( join('; ', @kw).'.' );
}



1;

__END__


BEGIN {
    my $eight_hundred_Ns = 'n' x 800;

    sub make_unfinished_embl {
        my( $info, $contigs ) = @_;

        my $project = $info->{'Sequence'};
        my $acc     = $info->{'Accno'};
        my $embl_id = $info->{'embl_id'} || 'ENTRYNAME';
        my $author  = $info->{'Author'};
        my $species = $info->{'speciesname'};
        my $status  = $info->{'status_code'};
        my( $ext_clone );
        {
            my $e = external_clone_name($project);
            $ext_clone = $e->{$project}
                or confiess "Can't make external clone name";
        }
        my $date = EMBLdate();
        my $binomial = species_binomial($species)
            or confiess "Can't get latin name for '$species'";

        # Make the sequence
        my( $dna, %contig_lengths );
        foreach my $contig (sort keys %{$contigs->{DNA}}) {
            my $con = \$contigs->{DNA}{$contig};
            $contig_lengths{$contig} = length($$con);
            $dna .= $eight_hundred_Ns if $dna;
            $dna .= $$con;
        }
        my $seqlength = length($dna);

        # New embl file object
        my $embl = Hum::EMBL->new();
        
        # ID line
        my $id = $embl->newID;
        $id->entryname($embl_id);
        $id->dataclass('standard');
        $id->molecule('DNA');
        $id->division('HTG'); ### Assume this is the same for other organisms
        $id->seqlength($seqlength);
        $embl->newXX;
        
        # AC line
        if ($acc) {
            my $ac = $embl->newAC;
            $ac->primary($acc);
            $embl->newXX;
        }
        
        # AC * line
        my $ac_star = $embl->newAC_star;
        my $identifier = '_'. uc $project;
        $ac_star->identifier($identifier);
        $embl->newXX;
    
        # DE line
        my $de = $embl->newDE;
        $de->list("$species DNA sequence *** SEQUENCING IN PROGRESS *** from clone $ext_clone");
        $embl->newXX;
        
        # KW line
        my $kw = $embl->newKW;
        my @kw_list = ('HTG', 'HTGS_PHASE1');
        push( @kw_list, 'HTGS_DRAFT' ) if $status == 30;
        $kw->list(@kw_list);
        $embl->newXX;
    
        # Organism
        add_Organism($embl, $species);
        $embl->newXX;
        
        # Reference
        my $ref = $embl->newReference;
        $ref->number(1);
        $ref->authors($author);
        $ref->locations("Submitted ($date) to the EMBL/Genbank/DDBJ databases.",
                        'Sanger Centre, Hinxton, Cambridgeshire, CB10 1SA, UK.',
                        'E-mail enquiries: humquery@sanger.ac.uk',
                        'Clone requests: clonerequest@sanger.ac.uk');
        $embl->newXX;
        
        # Comments
        my $unfin_cc = $embl->newCC;
        $unfin_cc->list(
"IMPORTANT: This sequence is unfinished and does not necessarily
represent the correct sequence.  Work on the sequence is in progress and
the release of this data is based on the understanding that the sequence
may change as work continues.  The sequence may be contaminated with
foreign sequence from E.coli, yeast, vector, phage etc.");
        $embl->newXX;
        
        my $contig_cc = $embl->newCC;
        $contig_cc->list(
            "Order of segments is not known; 800 n's separate segments.",
            map "Contig_ID: $_  Length: $contig_lengths{$_}bp",
            sort keys %contig_lengths );
        $embl->newXX;
    
        # Feature table source feature
        my( $chr, $map   ) = localisation_data(  $project );
        my( $libraryname ) = library_and_vector( $project );
        add_source_FT( $embl, $seqlength, $binomial, $ext_clone,
                       $chr, $map, $libraryname );            
        $embl->newXX;
    
        # Sequence
        my $sq = $embl->newSequence;
        $sq->seq($dna);
        
        $embl->newEnd;
        
        return $embl;
    }   
}
