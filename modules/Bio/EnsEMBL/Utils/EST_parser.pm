#
# BioPerl module for Bio::EnsEMBL::Utils::EST_parser
#
# Cared for by Elia Stupka <elia@sanger.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::EST_parser 

=head1 SYNOPSIS

    EST_parser->new();

    #To parse a EST file, build feature pairs:
    my @features=EST_parser->parse_file($file_handle);

=head1 DESCRIPTION

This module parses EST alignments produced by Jim Kent and generates 
feature pairs.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Utils::EST_parser;

use strict;
use vars qw(@ISA);
use Bio::Root::RootI;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
@ISA = qw(Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : EST_parser->new()
 Function: Constructor
 Example : 
 Returns : Reference to an object
 Args    : 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;

    $self->{'_feature_array'} = [];
    return $self;
}

=head2 parse_file

 Title   : parse_file
 Usage   : EST_parser->parse_file(filename)
 Function: Parses a EST file, reading in features
 Example : EST_parser->parse_file(gens.EST)
 Returns : array of Bio::EnsEMBL::Gene objects
 Args    : name of the file

=cut

sub parse_file {
    my ($self,$fh) = @_;

    $fh || $self->throw("You have to pass a filehandle in order to parse a file!");
    
    while( <$fh> ) {
	(!/^\d/) && next;
	#First, let's split all the columns in the file
	my @array= split;

	my $match=$array[0];
	my $mismatch=$array[1];

	#Not sure about this figure...
	if ($match == 0 && $mismatch == 0) {
	    print STDERR "Could not determine score for line $_\n";
	    print STDERR "Skipping!\n";
	    next;
	}
	my $id=($match/($match+$mismatch))*100;
	my $id=sprintf("%.4f",$id);
        #Then we get rid of the first 8 columns
	@array = splice(@array,8);

	#Then we get the strand
	my $strand;
	if ($array[0] eq '+') {
	    $strand = 1;
	}
	elsif ($array[0] eq '-') {
	    $strand = -1;
	}

	#The hseqname
	my $hseqname = $array[1];

	#Length of the gapped matches (needed for reverse strand features)
	my $e_match = $array[2];
	my $g_match = $array[6];
	
	#International Contig id
	my $seqname = $array[5];
	$seqname =~ s/\_/\~/g;
	#We get rid of the next columns which refer to gapped data
	@array = splice(@array,9);

	#Find out how many feature pairs we are going to build
	my $fn = @array[0];
	print STDERR "Number of feature pairs: $fn\n";
	shift @array;
	my $c;
	my @lengths=split(/\,/,$array[0]);
	my @e_starts=split(/\,/,$array[1]);
	my @g_starts=split(/\,/,$array[2]);
	for ($c=0;$c<$fn;$c++) {
	    my $length=$lengths[$c];
	    my $g_start=$g_starts[$c]+1;
	    my $g_end=$g_start+$length;
	    my $e_start;
	    my $e_end;
	    
	    if ($strand == 1) {
		$e_start=$e_starts[$c]+1;
		$e_end=$e_start+$length;
	    }
	    else {
		$e_start=$e_match-$e_starts[$c]+1;
		$e_end=$e_start-$length;
	    }
	    my $analysis = new Bio::EnsEMBL::Analysis(-db              => 'dbest',
						      -db_version      => 1,
						      -program         => '0-0greedy',
						      -program_version => 1,
						      -gff_source      => '0-0greedy',
						      -gff_feature     => 'similarity');
	    my $f1 = new Bio::EnsEMBL::SeqFeature(-seqname => $seqname,
						  -start   => $g_start,
						  -end     => $g_end,
						  -score   => $id,
						  -source_tag  =>'0-0greedy',
						  -primary_tag =>'similarity',
						  -strand      => 1,
						  -analysis    => $analysis,
						  -percent_id     => $id
						  );
	    
	    my $f2 = new Bio::EnsEMBL::SeqFeature(-seqname => $hseqname,
						  -start   => $e_start,
						  -end     => $e_end,
						  -score   => $id,
						  -source_tag  =>'est',
						  -primary_tag =>'similarity',
						  -strand      => $strand,
						  -analysis    => $analysis,
						  -percent_id     => $id
						  );
	    
	    print STDERR "Feature Pair ".($c+1).":\n";
	    print STDERR " seqname: $seqname\n";
	    print STDERR "hseqname: $hseqname\n";
	    print STDERR "  Strand: $strand\n";
	    print STDERR "Identity: $id\n";
	    print STDERR " Genomic: start $g_start, end $g_end\n";
	    print STDERR "     EST: start $e_start, end $e_end\n\n";
	    eval {
		$f1->validate();
		$f2->validate();
	    };
	    if ($@) {
		print STDERR "Could not validate $@\n";
		next;
	    }
	    	    
	    my $fp = new Bio::EnsEMBL::FeaturePair(-feature1 => $f1,
						   -feature2 => $f2);
	    push @{$self->{'_feature_array'}},$fp;
	}
    }
    return  @{$self->{'_feature_array'}};
}

=head2 print_ests

 Title   : print_genes
 Usage   : EST_parser->print_genes(@genes)
 Function: Prints gene structures to STDOUT (for tests)
 Example : 
 Returns : nothing
 Args    : 

=cut

sub print_ests {
    my ($self) = @_;
    
    my @features= @{$self->{'_feature_array'}};

    foreach my $feature (@features) {
	print STDOUT "Feature ".$feature->seqname." ".$feature->start." ".$feature->end." ".$feature->score." ".$feature->strand." ".$feature->analysis->db."\n";
    }
}


