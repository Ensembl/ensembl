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

use Bio::Tools::HMMER::Results;

use Bio::EnsEMBL::Analysis::GenscanPeptide;
use Bio::EnsEMBL::Analysis::Analysis;
use Bio::EnsEMBL::Analysis::MSPcrunch;
use Bio::EnsEMBL::Analysis::GFF;

use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Repeat;

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
    $seq                  || $self->throw("Cannot make contig_feature object without seq object");
    
    $self->id       ($id);
    $self->clone_dir($clone_dir);
    $self->disk_id  ($disk_id);
    $self->gs       ($gs);
    $self->seq      ($seq);

    $self->_debug(1) if $debug;
    $self->{_features} = [];   # This stores the features.
    $self->{_repeats}  = [];   # This stores the features.
    $self->{_genscan}  = [];   # This stores the genscan.
    
    # DEBUG
    print_genes($gs,$seq) if $self->_debug;

    return $make;
}

sub read_Similarities {
    my ($self) = @_;

    my $gs        = $self->gs;
    my $msptype   = $self->msptype;
    my $clone_dir = $self->clone_dir;
    my $disk_id   = $self->disk_id;
    my $id        = $self->id;

    # loop over transcripts
    my $count = 1;

    foreach my $g ($gs->each_Transcript) {
	my $genpep     = new Bio::EnsEMBL::Analysis::GenscanPeptide($g);
	
	MSP: foreach my $msp (@$msptype) {

	    next MSP if ($msp->[0] eq "ce_p" || $msp->[0] eq "sh_p");

	    my $mspfile        = "$clone_dir/$disk_id.$count".$msp->[4];
	    my $pfamfile       = "$clone_dir/$disk_id.$count".$msp->[4];

	    my $pid            = "$id.$count";
	    
	    if ($msp->[5]     eq 'msp'){
		eval {
		    $self->read_MSP($mspfile,$genpep,$msp);
		};
		if ($@) {
		    $self->warn("Error reading MSPFile $mspfile\n");
		}
		
	    } elsif ($msp->[5] eq 'pfam'){
		eval {
		    $self->read_Pfam($pfamfile,$genpep,$msp);
		};
		if ($@) {
		    $self->warn("Error reading Pfam file\n");
		}

	    } elsif ($msp->[5] eq 'gff'){
		
		
	    } else {
		$self->throw("no parser for $$msp[5] defined");
	    }
	    
	}
	
	my @homols = $genpep->each_Homol;      # Converts the hits from peptide into genomic coordinates
	
	foreach my $homol (@homols) {
	    $self->add_Feature($homol);
	}
	
	$count++;    
    }
}



sub read_Repeats {
    my ($self) = @_;
    
    my $clone_dir = $self->clone_dir;
    my $disk_id   = $self->disk_id;
    my $msp       = $self->msptype->[6];

    my $gfffile    = "$clone_dir/$disk_id".$msp->[4];    

    if (! -e $gfffile) {
	print(STDERR "   - No repeat file $gfffile  exists - Skipping repeats\n");
	return;
    } elsif (! -r $gfffile) {
	print(STDERR "   - Repeat file $gfffile unreadable - Skipping repeats\n");
	return;
    } else {
	print(STDERR "   - Reading RepeatMasker file $gfffile\n");
    }

    my $analysis   = new Bio::EnsEMBL::Analysis::Analysis( -program         => $msp->[1],
							   -program_version => 1,
							   -gff_source      => 'RepeatMasker',
							   -gff_feature     => 'repeat');
    

    my $GFF        = new Bio::EnsEMBL::Analysis::GFF(-file => $gfffile,
						     -type => 'Repeat');
	  

    foreach my $f ($GFF->each_Feature) {
	$f->analysis($analysis);
	$self->add_Feature($f);
    }
}




sub read_Pfam {
    my ($self,$pfamfile,$genscan_peptide,$pfam) = @_;

    if (! -e $pfamfile) {
	print(STDERR  "   - No pfam file $pfamfile exists - Skipping pfam\n");
	return;
    } elsif (! -r $pfamfile) {
	print(STDERR  "   - Pfam file $pfamfile unreadable - Skipping pfam\n");
	return;
    } else {
	print(STDERR "   - Reading pfam file $pfamfile\n");
    }
    
    my $pfamobj = new Bio::Tools::HMMER::Results(-file => $pfamfile,
						 -type => 'hmmpfam');

    my $analysis = new Bio::EnsEMBL::Analysis::Analysis(-db              => $pfam->[2],
							-db_version      => 1,
							-program         => $pfam->[1],
							-program_version => 1,
							-gff_source      => 'hmmpfam',
							-gff_feature     => 'similarity');


    my @homols;

    foreach my $dom ($pfamobj->each_Domain) {
	my $dom2 = $dom->feature2;

	my $f1 = new Bio::EnsEMBL::SeqFeature(-seqname => $self->id,
					      -start   => $dom->start,
					      -end     => $dom->end,
					      -score   => $dom->score,
					      -source_tag  =>'hmmpfam',
					      -primary_tag =>'similarity',
					      -strand      => 1,
					      -analysis    => $analysis,
					      );

	my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $dom2->seqname,
					      -start   => $dom2->start,
					      -end     => $dom2->end,
					      -score   => $dom2->score,
					      -source_tag  =>'hmmpfam',
					      -primary_tag =>'similarity',
					      -strand      => 1,
					      -analysis    => $analysis,
					      );

	my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
					       -feature2 => $f2);

	$genscan_peptide->add_pepHit($fp);

    }

}

sub read_Genscan {
    my ($self,) = @_;

    my $analysis = new Bio::EnsEMBL::Analysis::Analysis(-program     => "Genscan",
							-program_version => 1,
							-gff_source  => "genscan",
							-gff_feature => "exon",
							);
		    
    foreach my $trans ($self->gs->each_Transcript) {
	my $gene = new Bio::EnsEMBL::SeqFeature(-primary_tag => 'prediction');
	$gene->source_tag ('genscan');
	
	$gene->score(-1000);
	$gene->analysis($analysis);

	foreach my $ex ($trans->each_Exon) {

	    my $f = new Bio::EnsEMBL::SeqFeature(-start  => $ex->start,	
						 -end    => $ex->end,
						 -strand => $ex->strand);

	    $f->source_tag ('genscan');
	    $f->primary_tag('prediction');
	    if( ! defined $ex->seqname ) {
		$f->seqname($ex->entire_seq->id());
	    } else {
		$f->seqname($ex->seqname);
	    }

	    if (defined($ex->score)) {
		$f->score($ex->score);
	    } else {
		$f->score(-1000);
	    }
	    $f->analysis($analysis);
	    if( !defined $gene->seqname ) {
		$gene->seqname($f->seqname);
	    }

	    $gene->add_sub_SeqFeature($f,'EXPAND');
	}

	push(@{$self->{_genscan}},$gene);
    }
}
	
sub each_Genscan {
    my ($self) = @_;

    if (defined($self->{_genscan})) {
	return @{$self->{_genscan}};
    }
}

sub add_Feature {
    my ($self,$f) = @_;


    $self->throw("Feature must be Bio::EnsEMBL::SeqFeatureI in add_Feature") 
	unless $f->isa("Bio::EnsEMBL::SeqFeatureI");
    
    if ($f->isa("Bio::EnsEMBL::Repeat")) {
	push(@{$self->{_repeats}},$f);
    } else {
	push(@{$self->{_features}},$f);
    }

}

sub each_Feature {
    my ($self) = @_;

    if (defined($self->{_features})) {
	return @{$self->{_features}};
    } 
}

sub each_Repeat {
    my ($self) = @_;
    
    if (defined($self->{_repeats})) {
	return @{$self->{_repeats}};
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
    } elsif (! -r $mspfile) {
	print(STDERR "   - MSPcrunch file $mspfile unreadable. Skipping\n");
	return;
    } else {
	print(STDERR "   - Reading MSPcrunch file $mspfile\n");
    }
    

    my $mspobj  = new Bio::EnsEMBL::Analysis::MSPcrunch(-contig_id  => $self->id,
							-file       => $mspfile,
							-type       => $type,
							-source_tag => $msp->[1]);


    my ($type1,$type2) = $mspobj->get_types;

    foreach my $homol ($mspobj->each_Homol) { 

	$homol->source_tag($mspobj->source_tag);
	$homol->feature2->source_tag($mspobj->source_tag);

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
	    
	} else {
	    print(STDERR "      - unrecognised query type $type1\n");
	}
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

=head2 id

 Title   : id
 Usage   : $obj->id($clonedir)
 Function: 
 Returns : 
 Args    : string 


=cut

sub id {
    my ($self,$id) = @_;

    if (defined($id)) {
	$self->{_id} = $id;
    }
    return $self->{_id};
}

=head2 clone_dir

 Title   : clone_dir
 Usage   : $obj->clone_dir($clonedir)
 Function: 
 Returns : 
 Args    : string 


=cut

sub clone_dir {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_clone_dir} = $arg;
    }

    return $self->{_clone_dir};
}

=head2 disk_id

 Title   : disk_id
 Usage   : $obj->disk_id($diskid)
 Function: 
 Returns : 
 Args    : string 


=cut

sub disk_id {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_disk_id} = $arg;
    }

    return $self->{_disk_id};
}

=head2 gs

 Title   : gs
 Usage   : $obj->gs($genscan)
 Function: 
 Returns : 
 Args    : Bio::EnsEMBL::Analysis::Genscan


=cut

sub gs {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$arg->isa('Bio::EnsEMBL::Analysis::Genscan') ||   $self->throw("$arg is not a gs object in new contig_feature");	
	$self->{_gs} = $arg;
    }

    return $self->{_gs};
}

=head2 seq

 Title   : seq
 Usage   : $obj->seq($seq)
 Function: 
 Returns : 
 Args    : 


=cut

sub seq {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$arg->isa('Bio::Seq') || $self->throw("$arg is not a seq object in new contig_feature");
	$self->{_seq} = $arg;
    }

    return $self->{_seq};
}

=head2 msptype

 Title   : msptype
 Usage   : $obj->msptype
 Function: 
 Returns : [][]
 Args    : 


=cut

sub msptype {
    my ($self) = @_;
    
    # mapping of data to filenames
    my $msptype = [['swir_p',  'blastp',  'swir',     'pep', '.blastp_swir.msptmp',   'msp' ,'PEP-PEP' ],
		   ['ce_p',    'tblastn', 'ce',       'dna', '.tblastn_ce.msptmp',    'msp'  ,'PEP-DNA'],
		   ['vert_p',  'tblastn', 'vert',     'dna', '.tblastn_vert.msptmp',  'msp'  ,'PEP-DNA' ],
		   ['sh_p',    'tblastn', 'sh',       'dna', '.tblastn_sh.msptmp',    'msp'  ,'PEP-DNA' ],
		   ['dbest_p', 'tblastn', 'dbest',    'dna', '.tblastn_dbest.msptmp', 'msp'  ,'PEP-DNA' ],
		   ['pfam_p',  'hmmpfam', 'PfamFrag', 'pep', '.hmmpfam_frag',         'pfam' ,'PEP-PEP' ],
		   ['repeat',  'RepeatMasker', '',    'dna', '.RepMask.out.gff',      'gff'  ,'DNA-DNA'],
		   ];
    return $msptype;

}

1;

