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
    $id || $self->throw("Cannot make contig_feature object without id");
    $clone_dir || $self->throw("Cannot make contig_feature object without clone_dir");
    $disk_id || $self->throw("Cannot make contig_feature object without disk_id");
    $gs || $self->throw("Cannot make contig_feature object without gs object");
    $gs->isa('Bio::EnsEMBL::Analysis::Genscan') || 
	$self->throw("$gs is not a gs object in new contig_feature");
    $seq || $self->throw("Cannot make contig_feature object without seq object");
    $seq->isa('Bio::Seq') || 
	$self->throw("$seq is not a seq object in new contig_feature");
    $self->_debug(1) if $debug;

    # read 2 types of features
    # 1. features aligned against contigs (no remapping required)
    # 2. features aligned against transcripts (remapping required)
    
    # DEBUG
    print_genes($gs,$seq) if $self->_debug;
  
    # mapping of data to filenames
    my $msptype = [['swir_p',  'blastp',  'swir',     'pep', '.blastp_swir.msptmp',   'msp'  ],
		   ['ce_p',    'tblastn', 'ce',       'dna', '.tblastn_ce.msptmp',    'msp'  ],
		   ['vert_p',  'tblastn', 'vert',     'dna', '.tblastn_vert.msptmp',  'msp'  ],
		   ['sh_p',    'tblastn', 'sh',       'dna', '.tblastn_sh.msptmp',    'msp'  ],
		   ['dbest_p', 'tblastn', 'dbest',    'dna', '.tblastn_dbest.msptmp', 'msp'  ],
		   ['pfam_p',  'hmmpfam', 'PfamFrag', 'pep', '.hmmpfam_frag',         'pfam' ],
		   ];
  
    # loop over transcripts
    my $count=1;
    foreach my $g ($gs->each_Transcript) {
	
	print_gene_details($g,$count) if $self->_debug;
    
	foreach my $msp (@$msptype) {
      
	    my $mspfile        = "$clone_dir/$disk_id.$count".$$msp[4];
	    my $pid            = "$id.$count";
	    my @homols;
	    if($$msp[5] eq 'msp'){
		@homols        = $self->_read_MSP($mspfile,$seq,$pid,$msp);
	    }elsif($$msp[5] eq 'pfam'){
		@homols        = $self->_read_pfam($mspfile,$seq,$pid,$msp);
	    }else{
		$self->throw("no parser for $$msp[5] defined");
	    }
	    my $offset         = $self->_get_offset($gs,$g,$count);
      
	    print STDERR "read MSP file $mspfile:\n  $#homols homols. Offset $offset\n" 
		if $self->_debug;
      
	    # these are all matches to transcripts, so need to remap to contig coordinates
	    foreach my $h (@homols) {

		my $hsf=$h->homol_SeqFeature;
		my $seq;
		print STDERR $h->seqname . "\t" . $h->start . "\t" . $h->end . "\t" . 
		    $hsf->start . "\t" . $hsf->end . "\t" . 
			$hsf->strand . "\t" . join(',',$hsf->each_tag_value('title')) . "\t" . 
			$seq ."\n";

		# Converts peptide to dna coords
		my @newhomols=$self->map_homols($h,$g,$offset,$seq,\*STDOUT); 
		
		print STDERR "Homols mapped to dna coords are :\n" if $self->_debug;

		foreach my $hh (@newhomols){
		    my $hsf=$hh->homol_SeqFeature;
		    my $seq;
		    #my $tmp=$g->translate_region($hh->start,$hh->end);
		    #$seq=$tmp->[0]->seq;
		    print STDERR $hh->seqname . "\t" . $hh->start . "\t" . $hh->end . "\t" . 
			$hsf->start . "\t" . $hsf->end . "\t" . $hsf->strand . "\t" . 
			    $seq ."\n";
		}	
	    }
	}
	$count++;    
    }
    return $make;
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

sub map_coords {
    my ($self,$gene,$homol,$start,$end,$realseq) = @_;
  
    my @homols;
  
    my $foundstart = 0;
    my $foundend   = 0;

    my @exons = $gene->each_Exon ; 
    my $strand = $exons[0]->strand;
    my $count = 0;
    
    my ($starts,$ends) = $gene->pep_coords;

    my $hoffset;

    # We have the start and end points for the peptide on the DNA
    # We now have to return exon sized chunks that make up
    # the peptide homology

    if ($strand == 1) {
	foreach my $ex ($gene->each_Exon) {
	    my $tmpstart;
	    my $tmpend;
	    
	    # look for start if not found
	    if ($foundstart == 0 ) {
		if ($start >= $ex->start && $start <= $ex->end) {
		    $foundstart = 1;
		    $tmpstart = $start;
		}
	    }
	    
	    # look for end if found start
	    if (!$foundend && $foundstart) {
	
		# Adjust to be in the right frame
		$tmpstart +=  (3 - ($start - $ex->phase - $ex->start)%3 ) % 3;
	
		$tmpend   = $ex->end;
	
		# Before making a new exon check we haven't also
		# Got the end point in the same exon
	
		if ($end <= $ex->end) {
		    $tmpend = $end;
		    $foundend = 1;
		}
	
		# Make sure the end coordinate is in the right frame.
		$tmpend -= ($tmpend - $ex->phase - 2 - $ex->start)%3;
	    }
      
	    if (defined($tmpstart) && defined($tmpend)) {
	
		# Find the peptide coord
		my $pstart = int(($tmpstart - $ex->start - $ex->phase)/3 + $starts->[$count]);
		my $pend   = int(($tmpend   - $ex->start - $ex->phase)/3 + $starts->[$count]);
	
		$hoffset = $pstart unless $hoffset;
	
		print STDERR "hoffset = $hoffset\n" if $self->_debug;
		print STDERR "pstart $pstart : $start " . $ex->start . " " . $ex->phase . " " . 
		    $starts->[$count] . " " . $ends->[$count] . "\n" if $self->_debug;
	
		my $h=Bio::SeqFeature::Homol->new(-start=>$tmpstart,
						  -end=>$tmpend,
						  -strand=>1,
						  );
		# hsf object contains matched sequence, title
		my $homolsf=$homol->homol_SeqFeature;
		my $hsf;
		$hsf=Bio::SeqFeature::Generic->new(-start=>
						   ($homolsf->start - $hoffset + $pstart),
						   -end=>
						   ($homolsf->end,  + $pend    - $pstart),
						   -strand=>$homolsf->strand,
						   );

		# standard values (this time 'dna' not 'pep')
		$h->add_tag_value('seqtype','dna');
		$h->source_tag('ensembl');
		$h->primary_tag('similarity');

		# copy standard tags
		$h->seqname($homol->seqname);
		$h->score($homol->score);
		$h->add_tag_value('percentid',
				  join(',',
				       $homol->each_tag_value('percentid')
				       )
				  );
		$h->add_tag_value('sim_label',
				  join(',',
				       $homol->each_tag_value('sim_label')
				       )
				  );
		$h->add_tag_value('method',
				  join(',',
				       $homol->each_tag_value('method')
				       )
				  );
		$h->attach_seq($homol->entire_seq);

		# copy standard tags
		$hsf->add_tag_value('title',
				    join(',',
					 $homolsf->each_tag_value('title')
					 )
				    );
		$hsf->add_tag_value('db',
				    join(',',
					 $homolsf->each_tag_value('db')
					 )
				    );
		$hsf->add_tag_value('seqtype',
				    join(',',
					 $homolsf->each_tag_value('seqtype')
					 )
				    );

		# link features
		$h->homol_SeqFeature($hsf);

		# add to array
		push(@homols,$h);
	
	    }
	    $count++;
	}
	return @homols;
    } else {

    foreach my $ex ($gene->each_Exon) {
      my $tmpstart;
      my $tmpend;
      
      if ($foundstart == 0 ) {
	if ($end <=  $ex->end && $end >= $ex->start) {
	  $foundstart = 1;
	  $tmpstart   = $end;
	  
	}
      }
      
      if (!$foundend && $foundstart) {
	$tmpstart = $ex->end unless $tmpstart;
	$tmpend   = $ex->start;
	
	# Adjust tmpstart to be in the right frame
	$tmpstart -=  (3 - ($ex->end - $end - $ex->phase)%3 ) % 3;

	# Before making a new exon check we haven't also
	# Got the end point in the same exon
      
	if ($start >= $ex->start) {
	  $tmpend = $start;
	  $foundend = 1;
	} elsif ($count < $#exons && $start >= $exons[$count+1]->end) {
	  $tmpend = $ex->start;
	  $foundend = 1;
	}
	# Make sure the end coordinate is in the right frame.
	$tmpend += ($tmpend - $ex->phase - 2 - $ex->end)%3;
      }
      

      if (defined($tmpstart) && defined($tmpend)) {
#	print("Exon $count : " . $ex->start . "\t" . $ex->end . "\t" . $tmpstart . "\t" . $tmpend . "\n");
	# Find the peptide coord

	my $pstart = int(($ex->end - $tmpstart  - $ex->phase)/3 + $starts->[$count]);
	my $pend   = int(($ex->end - $tmpend    - $ex->phase)/3 + $starts->[$count]);
	
	$hoffset = $pstart unless $hoffset;
	
#	print "hoffset = $hoffset\n";
#	print("pstart $pstart : $tmpstart " . $ex->end . " " . $ex->phase . " " . $starts->[$count] . " " . $ends->[$count] . "\n");
	
	# Now make new homol and add to the array
	my $h = Spangle::Homol->new();
	
	$h->id    ($homol->id);
	$h->score ($homol->score);
	
	$h->start ($tmpend);
	$h->end   ($tmpstart);
	
	$h->hend  ($homol->hstart - $hoffset + $pstart);
	$h->hstart($h->hend     + $pend    - $pstart);
	
	$h->strand($homol->strand);
	
	push(@homols,$h);
	
      }
      $count++;
    }
    return @homols;
  }
}

sub map_homols {
    my($self,$h,$g,$offset,$seq,$fh) = @_;

    # The start and end coords we have in the homol objects
    # are genscan peptide start and end points
  
    # First convert them using $offset into transcript peptide
    # coordinates

    my @exon   = $g->each_Exon;
    my $strand = $exon[0]->strand;

    my $wholeseq = $g->translate->seq;

    $h->start ($h->start  - $offset);
    $h->end   ($h->end    - $offset);
    my $hsf=$h->homol_SeqFeature;
  
    # Adjust the strand
    # (two strand: one for h object, one for hsf object)
    if ($strand == -1) {
	print STDERR "Strand is -1 " . $hsf->strand . "\n" if $self->_debug;
	$h->strand($strand);
	if($hsf->strand==1){
	    $hsf->strand(-1);
	}elsif($hsf->strand==-1){
	    $hsf->strand(1);
	}
    }
    print STDERR "New strand is " . $hsf->strand . "\n" if $self->_debug;

    # Now we need to convert the 'true' peptide coords into 
    # genomic dna coords.

    my $startdna;
    my $enddna;


    print STDERR "Finding dna coords for ". $h->start . "\t" . $h->end . "\n"
	if $self->_debug;
    if ($strand == 1) {
	$startdna = $g->find_coord($h->start,"start");
	$enddna   = $g->find_coord($h->end,"end");
    } else {
	$enddna   = $g->find_coord($h->start,"start");
	$startdna = $g->find_coord($h->end  ,"end");
    }
    print STDERR "Translating region $startdna $enddna\n"
	if $self->_debug;

    #my $wholetr = $g->translate_region($startdna,$enddna);
    #for (my $i = 0; $i < 1; $i++) {
	#print("Whole peptide is " . $wholetr->[$i]->seq . "\n");
    #}
    
#  print("Real seq is      " . $realseq   . "\n");
#  print("Real start end   " . $realstart . " " . $realend .  "\n");  

    # Translate this match into a peptide for each exon
    my @newh = $self->map_coords($g,$h,$startdna,$enddna,$seq);
    return @newh;
}


sub _get_offset {
    my ($self,$gs,$g,$count) = @_;

    # translation of transcript (
    my $pep    = $g->translate()->seq;
    # peptide that was actually searched (genscan) [may be longer due to initial X]
    my $peps   = $gs->{_peptides}[$count-1]->seq;
    my $offset = index($peps,$pep);

    # error if does not match
    $self->throw("$pep not found in $peps") if $offset==-1;
    
#   print("PEP " . $peps . " " . (length($peps)) . "\n\n");
  
    return $offset;
}  

sub _read_MSP {

    my($self,$mspfile, $seq, $id, $msp) = @_;
    my($sim_label,$method,$db,$seqtype,$ext)=@$msp;
    my @homols;

    open(MSP,$mspfile) || $self->throw("Couldn't open $mspfile");

    while (my $line = <MSP>) {
	unless ($line =~ /^\#/) {
	    my ($score,$pid,$start,$end,$id,$hstart,$hend,$hid,$title) = split(' ',$line,9);

	    # using Bioperl objects to represent homology
	    # homol object contains feature on peptide, score, pid
	    # and is attached to underlying sequence
	    # strand is always 1 since sequence was peptide
	    my $homol=Bio::SeqFeature::Homol->new(-start=>$start,
						  -end=>$end,
						  -strand=>1,
						  );
	    $homol->seqname($id);
	    $homol->add_tag_value('seqtype','pep');

	    $homol->source_tag('ensembl');
	    $homol->primary_tag('similarity');
	    $homol->score($score);
	    $homol->add_tag_value('percentid',$pid);
	    $homol->add_tag_value('sim_label',$sim_label);
	    $homol->add_tag_value('method',$method);
	    $homol->attach_seq($seq);

	    # hsf object contains matched sequence, title
	    my $hsf;
	    if($hstart>$hend){
		$hsf=Bio::SeqFeature::Generic->new(-start=>$hend,
						   -end=>$hstart,
						   -strand=>-1,
						   );
	    }else{
		$hsf=Bio::SeqFeature::Generic->new(-start=>$hstart,
						   -end=>$hend,
						   -strand=>1,
						   );
	    }
	    $hsf->add_tag_value('title',$title);
	    $hsf->add_tag_value('db',$db);
	    $hsf->add_tag_value('seqtype',$seqtype);
	    $homol->homol_SeqFeature($hsf);
	    
	    push(@homols,$homol);
	}
    }
    return @homols;
}

sub _read_pfam {

    my($self,$mspfile, $seq, $id, $msp) = @_;
    my($sim_label,$method,$db,$seqtype,$ext)=@$msp;
    my @homols;
    return @homols;

    open(MSP,$mspfile) || $self->throw("Couldn't open $mspfile");

    while (my $line = <MSP>) {
	unless ($line =~ /^\#/) {
	    my ($score,$pid,$start,$end,$id,$hstart,$hend,$hid,$title) = split(' ',$line,9);

	    # using Bioperl objects to represent homology
	    # homol object contains feature on peptide, score, pid
	    # and is attached to underlying sequence
	    # strand is always 1 since sequence was peptide
	    my $homol=Bio::SeqFeature::Homol->new(-start=>$start,
						  -end=>$end,
						  -strand=>1,
						  );
	    $homol->seqname($id);
	    $homol->add_tag_value('seqtype','pep');

	    $homol->source_tag('ensembl');
	    $homol->primary_tag('similarity');
	    $homol->score($score);
	    $homol->add_tag_value('percentid',$pid);
	    $homol->add_tag_value('sim_label',$sim_label);
	    $homol->add_tag_value('method',$method);
	    $homol->attach_seq($seq);

	    # hsf object contains matched sequence, title
	    my $hsf;
	    if($hstart>$hend){
		$hsf=Bio::SeqFeature::Generic->new(-start=>$hend,
						   -end=>$hstart,
						   -strand=>-1,
						   );
	    }else{
		$hsf=Bio::SeqFeature::Generic->new(-start=>$hstart,
						   -end=>$hend,
						   -strand=>1,
						   );
	    }
	    $hsf->add_tag_value('title',$title);
	    $hsf->add_tag_value('db',$db);
	    $hsf->add_tag_value('seqtype',$seqtype);
	    $homol->homol_SeqFeature($hsf);
	    
	    push(@homols,$homol);
	}
    }
    return @homols;
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

