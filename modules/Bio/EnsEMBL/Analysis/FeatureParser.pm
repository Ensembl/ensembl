#
# BioPerl module for Bio::EnsEMBL::Analysis::FeatureParser
#
# Cared for by Tim Hubbard <th@sanger.ac.uk>
#
# Copyright Tim Hubbard, Michele Clamp
#
# (based on ensembl/scripts/test_write_features
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::FeatureParser - Perl wrapper over feature flat files

=head1 SYNOPSIS

    $sfobj = Bio::EnsEMBL::Analysis::FeatureParser->new($dir,$contig,$gs,$seq);

    foreach my $sf ($sfobj->each_feature){
	my($start,$end,$strand,$score,
	   $name2,$start2,$end2,$pid,$method)=@$sf;
    }

=head1 DESCRIPTION

For each transcript in genscan_object passed, reads feature files, remaps from
gs peptide to dna coordinates.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::Analysis::FeatureParser;
use vars qw($AUTOLOAD @ISA);
use strict;

use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Homol;

use Bio::Tools::HMMER::Results;

use Bio::EnsEMBL::Analysis::GenscanPeptide;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::EnsEMBL::Analysis::Repeat;
use Bio::EnsEMBL::Analysis::GFF;

use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
    my($self,@args) = @_;
  
    my $make = $self->SUPER::_initialize;
    my ($id,$clone_dir,$disk_id,$gs,$seq,$debug)=@args;

    $id                   || $self->throw("Cannot make contig_feature object without id");
    $clone_dir            || $self->throw("Cannot make contig_feature object without clone_dir");
    $disk_id              || $self->throw("Cannot make contig_feature object without disk_id");
    $gs                   || $self->throw("Cannot make contig_feature object without gs object");

    $gs->isa('Bio::EnsEMBL::Analysis::Genscan') ||   $self->throw("$gs is not a gs object in new contig_feature");
    $seq                  || $self->throw("Cannot make contig_feature object without seq object");
    $seq->isa('Bio::Seq') || $self->throw("$seq is not a seq object in new contig_feature");

    $self->_debug(1) if $debug;


    $self->{_features} = [];   # This stores the features.

    $self->id($id);

    # read 2 types of features
    # 1. features aligned against contigs (no remapping required)
    # 2. features aligned against transcripts (remapping required)
    

    # DEBUG
    print_genes($gs,$seq) if $self->_debug;
  
    # mapping of data to filenames
    my $msptype = [['swir_p',  'blastp',  'swir',     'pep', '.blastp_swir.msptmp',   'msp' ,'PEP-PEP' ],
		   ['ce_p',    'tblastn', 'ce',       'dna', '.tblastn_ce.msptmp',    'msp'  ,'PEP-DNA'],
		   ['vert_p',  'tblastn', 'vert',     'dna', '.tblastn_vert.msptmp',  'msp'  ,'PEP-DNA' ],
		   ['sh_p',    'tblastn', 'sh',       'dna', '.tblastn_sh.msptmp',    'msp'  ,'PEP-DNA' ],
		   ['dbest_p', 'tblastn', 'dbest',    'dna', '.tblastn_dbest.msptmp', 'msp'  ,'PEP-DNA' ],
		   ['pfam_p',  'hmmpfam', 'PfamFrag', 'pep', '.hmmpfam_frag',         'pfam' ,'PEP-PEP' ],
		   ['repeat',  'RepeatMasker', '',    'dna', '.RepMask.out.gff',      'gff'  ,'DNA-DNA'],
		   ];

    # MC. I've bypassed the rest of this code as 
    #        1. The feature reading doesn't work
    #        2. The feature reading is being rewritten in as slightly different way
    #        3. We're only reading repeats at the moment.

    $self->read_Genscan($gs);

    # loop over transcripts
    my $count = 1;

    foreach my $g ($gs->each_Transcript) {
	my $genpep     = new Bio::EnsEMBL::Analysis::GenscanPeptide($g);
	print_gene_details($g,$count) if $self->_debug;
    
	foreach my $msp (@$msptype) {
      
	    my $mspfile        = "$clone_dir/$disk_id.$count".$msp->[4];
	    my $pid            = "$id.$count";

	    if ($msp->[5]     eq 'msp'){
		$self->read_MSP($mspfile,$genpep,$msp);
	    } elsif ($msp->[5] eq 'pfam'){
		$self->read_Pfam($genpep,$clone_dir,$disk_id,$count,$msp);
	    } elsif ($msp->[5] eq 'gff'){
		$self->read_Repeats($clone_dir,$disk_id,$msptype->[6]);
	    } else {
		$self->throw("no parser for $$msp[5] defined");
	    }
	    
	}

	$count++;    
    }

    return $make;
}

sub id {
    my ($self,$id) = @_;

    if (defined($id)) {
	$self->{_id} = $id;
    }
    return $self->{_id};
}

sub read_Repeats {
    my ($self,$clone_dir,$disk_id,$msp) = @_;

    my $gfffile    = "$clone_dir/$disk_id".$msp->[4];    

    if (! -e $gfffile) {
	print(STDERR "   - No repeat file $gfffile  exists - Skipping repeats\n");
	return;
    } else {
	print(STDERR "   - Reading RepeatMasker file $gfffile\n");
    }

    my $GFF        = new Bio::EnsEMBL::Analysis::GFF(-file => $gfffile,
						     -type => 'Repeat');
	    
    foreach my $f ($GFF->each_Feature) {
	$self->add_Feature($f);
    }
}




sub read_Pfam {
    my ($self,$genscan_peptide,$clone_dir,$disk_id,$count,$pfam) = @_;

    my $pfamfile  = "$clone_dir/$disk_id.$count". $pfam->[4];    

    if (! -e $pfamfile) {
	print(STDERR "   - No pfam file $pfamfile exists - Skipping pfam\n");
	return;
    } else {
	print(STDERR "   - Reading pfam file $pfamfile\n");
    }
    
    my $pfamobj = new Bio::Tools::HMMER::Results(-file => $pfamfile,
						 -type => 'hmmpfam');

    my @homols;

    foreach my $dom ($pfamobj->each_Domain) {
	my $align     = new Bio::EnsEMBL::Analysis::PairAlign;
	
	$dom->source_tag ('hmmpfam');
	$dom->primary_tag('similarity');
	$dom->strand    (1);

	my $dom2 = $dom->homol_SeqFeature;

	$dom2->source_tag('hmmpfam');
	$dom2->primary_tag('similarity');
	$dom2->strand(1);

	$genscan_peptide->add_pepHit($dom);

    }
    my @newh = $genscan_peptide->each_Homol;

    foreach my $h (@newh) {
	$self->add_Feature($h);
    }

}
sub read_Genscan {
    my ($self,$genscan) = @_;

    foreach my $trans($genscan->each_Transcript) {
	foreach my $ex ($trans->each_Exon) {
	    my $f = new Bio::SeqFeature::Homol(-start  => $ex->start,
					       -end    => $ex->end,
					       -strand => $ex->strand);

	    $f->source_tag($ex->source_tag);
	    $f->primary_tag($ex->primary_tag);
	    $f->seqname($ex->seqname);

	    if (defined($ex->score)) {
		$f->score($ex->score);
	    }

	    $self->add_Feature($f);
	}
    }
}
	
sub add_Feature {
    my ($self,$f) = @_;

    $self->throw("Feature must be Bio::SeqFeature::Generic in add_Feature") unless $f->isa("Bio::SeqFeature::Generic");

    if (!defined($self->{_features})) {
	$self->{_features} = [];
	$self->warn("The feature array does not exist!! Creating an empty one");
    }
    my @f = @{$self->{_features}};

    push(@{$self->{_features}},$f);
}

sub each_Feature {
    my ($self) = @_;

    if (defined($self->{_features})) {
	return @{$self->{_features}};
    } 
}


		
sub print_gene_details {
    my ($g,$count) = @_;
  
    print STDERR "Genscan predicted gene number $count\n";
    print STDERR "DNA exon coordinates (phase, peptide coords, exon peptide seq) are :\n";

    my ($starts,$ends) = $g->pep_coords;
    my $c = 0;
    foreach my $ex ($g->each_Exon) {
	my $seq = $ex->translate;
	print STDERR $ex->start . "\t" . $ex->end . "\t" . $ex->phase . "\t" . 
	    $starts->[$c] . "\t" . $ends->[$c] . "\t" . $seq->seq() . "\n";
	$c++;
    }
    print STDERR "\n";
}


sub read_MSP {

    my($self,$mspfile, $genpep, $msp) = @_;

    my $type = $msp->[6];

    if (! -e $mspfile) {
	print(STDERR "   - MSPcrunch file $mspfile doesn't exist. Skipping\n");
	return;
    } else {
	print(STDERR "   - Reading MSPcrunch file $mspfile\n");
    }
    
    print(STDERR "Setting source tag to " . $msp->[1] . " for $mspfile\n");
    my $mspobj  = new Bio::EnsEMBL::Analysis::MSPcrunch(-file => $mspfile,
							-type => $type,
							-source_tag => $msp->[1]);

    my ($type1,$type2) = $mspobj->get_types;

    foreach my $homol ($mspobj->each_Homol) { 
	if ($type1 eq "PEP") {
	    if ($type2 eq "DNA") {
		$genpep->add_dnaHit($homol);
	    } elsif ($type2 eq "PEP") {
		$genpep->add_pepHit($homol);
	    } else {
		print(STDERR "      - unrecognised homol type $type2\n");
	    }
	} elsif ($type1 eq "DNA") {
	    $self->add_Feature($homol);
	    
	}else {
	    print(STDERR "      - unrecognised query type $type1\n");
	}
    }

    my @homols = $genpep->each_Homol;      # Converts the hits from peptide into genomic coordinates

    foreach my $h (@homols) {
	if ($mspobj->source_tag eq "") {
	    print(STDERR "ERROR: Empty string for source tag for $mspfile\n");
	}
	#if ($h->homol_SeqFeature->seqname eq "TR:Q14839") {
	#    print(STDERR "Setting source tag for " . $h->homol_SeqFeature->seqname . " to " . $mspobj->source_tag . "\n");
	#}
	$h->source_tag($mspobj->source_tag);
	$self->add_Feature($h);
    }
}


sub print_genes {
    my ($gs,$seq) = @_;
  
    print("\nSequence id : " . $seq->id() . "\n");
  
    my $count = 1;
    foreach my $gene ($gs->each_Transcript) {
	$gene->contig_dna($seq);
	print("\nGene number  $count\n");
	my $trans = $gene->translate();
	print("GENE $count " . $seq->id() . " " . $trans->seq() . "\n");
	
	$count++;
    }
}


=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
    my $obj = shift;
    if( @_ ) {
	my $value = shift;
	$obj->{'_debug'} = $value;
    }
    return $obj->{'_debug'};
}
1;

