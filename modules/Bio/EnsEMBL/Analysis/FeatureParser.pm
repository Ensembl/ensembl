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

    $sfobj = Bio::EnsEMBL::Analysis::FeatureParser->new($dir,$contig,$genscan_object);

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
use Spangle::Homol;
use FileHandle;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;
  my ($clone_dir,$disk_id,$gs,$seq)=@args;
  $clone_dir || $self->throw("Cannot make contig object without clone_dir");
  $disk_id || $self->throw("Cannot make contig object without disk_id");
  $gs || $self->throw("Cannot make contig object without gs object");
  $gs->isa('Bio::EnsEMBL::Analysis::Genscan') || 
      $self->throw("Cannot make contig object with a $gs object");

  # DEBUG
  #print_genes($gs,$seq);
  
  # mapping of data to filenames
  my $msphash = { blastp        => ".blastp_swir.msptmp",
		  tblastn       => ".tblastn.msptmp",
		  tblastn_ce    => ".tblastn_ce.msptmp",
		  tblastn_vert  => ".tblastn_vert.msptmp",
		  tblastn_sh    => ".tblastn_sh.msptmp",
		  tblastn_dbest => ".tblastn_dbest.msptmp",
	      };
  
  # loop over transcripts
  my $count=1;
  foreach my $g ($gs->each_Transcript) {
    
      #print_gene_details($g,$count);
    
      foreach my $msp (keys %$msphash) {
      
	  my $mspfile        = "$clone_dir/$disk_id.$count" . $msphash->{$msp} ;
	  my @homols         = read_MSP  ($mspfile);    
	  my $offset         = get_offset($gs,$g,$count);
      
	  print("read MSP file $mspfile: $#homols homols. Offset $offset\n");
      
      
	  foreach my $h (@homols) {
	      # Converts peptide to dna coords
	      my @newhomols =  map_homols($h,$g,$offset,$seq,\*STDOUT); 
	
	      print("Homols mapped to dna coords are :\n");
	      foreach my $hh (@newhomols) {
		  my $tmp = $g->translate_region($hh->start,$hh->end);
		  print($hh->id . "\t" . $hh->start . "\t" . $hh->end . "\t" . 
			$hh->hstart . "\t" . $hh->hend . "\t" . $hh->frame . "\t" . 
			$tmp->[0]->seq ."\n");
	      }	
	  }
      }
      $count++;    
  }
  return $make;
}

sub print_gene_details {
  my ($g,$count) = shift;
  
  print("Genscan predicted gene number $count\n");
  print("DNA exon coordinates (phase,peptide coords, exon peptide seq) are :\n");

  my ($starts,$ends) = $g->pep_coords;
  my $c = 0;
  
  foreach my $ex ($g->each_Exon) {
    my @seq = $g->translate_exon($ex);
    
    print($ex->start . "\t" . $ex->end . "\t" . $ex->phase . "\t" . $starts->[$c] . "\t" . $ends->[$c] . "\t" . $seq[$ex->phase]->seq . "\n");
    $c++;
  }

  print("\n");
}

sub map_coords {
  my ($gene,$homol,$start,$end,$realseq) = @_;
  
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
      
      if ($foundstart == 0 ) {
	if ($start >= $ex->start && $start <= $ex->end) {
	  $foundstart = 1;
	  $tmpstart = $start;
	  
	  
	}
      }
      
      if (!$foundend && $foundstart) {
	$tmpstart = $ex->start unless $tmpstart;
	
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
	
#	print "hoffset = $hoffset\n";
#	print("pstart $pstart : $start " . $ex->start . " " . $ex->phase . " " . $starts->[$count] . " " . $ends->[$count] . "\n");
	
	# Now make new homol and add to the array
	my $h = Spangle::Homol->new();
	
	$h->id    ($homol->id);
	$h->score ($homol->score);
	
	$h->start ($tmpstart);
	$h->end   ($tmpend);
	
	$h->hstart($homol->hstart - $hoffset + $pstart);
	$h->hend  ($h->hstart     + $pend    - $pstart);
	
	$h->strand($homol->strand);
	
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
  my ($h,$g,$offset,$seq,$fh) = @_;

  # The start and end coords we have in the homol objects
  # are genscan peptide start and end points
  
  # First convert them using $offset into transcript peptide
  # coordinates

  my @exon   = $g->each_Exon;
  my $strand = $exon[0]->strand;

  my $wholeseq = $g->translate->seq;

  $h->start ($h->start  - $offset);
  $h->end   ($h->end    - $offset);
  
  # Adjust the strand
  if ($strand == -1) {
    print("Strand is -1 " . $h->strand . "\n");
    if ($h->strand ==  1) {
      $h->strand(-1);
    } elsif ($h->strand == -1) {
      $h->strand( 1);
    }
  }

  print("New strand is " . $h->strand . "\n");
  # Now we need to convert the 'true' peptide coords into 
  # genomic dna coords.

  my $startdna;
  my $enddna;


  print("Finding dna coords for ". $h->start . "\t" . $h->end . "\n");
  if ($strand == 1) {
    $startdna = $g->find_coord($h->start,"start");
    $enddna   = $g->find_coord($h->end,"end");
  } else {
    $enddna   = $g->find_coord($h->start,"start");
    $startdna = $g->find_coord($h->end  ,"end");
  }

  print("Translating region $startdna $enddna\n");
  my $wholetr = $g->translate_region($startdna,$enddna);

  for (my $i = 0; $i < 1; $i++) {
    print("Whole peptide is " . $wholetr->[$i]->seq . "\n");
  }

#  print("Real seq is      " . $realseq   . "\n");
#  print("Real start end   " . $realstart . " " . $realend .  "\n");  

  # Translate this match into a peptide for each exon
  my @newh = map_coords($g,$h,$startdna,$enddna,$seq);
  return @newh;
}


sub get_offset {
  my ($gs,$g,$count) = @_;

  my $pep    = $g->translate()->seq;
  my $peps   = $gs->{_peptides}[$count-1]->seq;
  my $offset = index($peps,$pep);

#  print("PEP " . $peps . " " . (length($peps)) . "\n\n");
  
  return $offset;
}  

sub read_MSP {

  my ($mspfile) = @_;
  my @homols;

  open(MSP,"<$mspfile") || die ("Couldn't open $mspfile");

  while (my $line = <MSP>) {
    unless ($line =~ /^\#/) {
      my ($score,$pid,$start,$end,$id,$hstart,$hend,$hid,$title) = split(' ',$line,9);

      # using Bioperl objects
      my $homol=new Bio::SeqFeature::Homol(-start=>$start,
					   -end=>$end,
					   -strand=>1,
					   );
      my $hsf;
      if($hstart>$hend){
	  $hsf=new Bio::SeqFeature(-start=>$hend,
				   -end=>$hstart,
				   -strand=>-1,
				   );
      }else{
	  $hsf=new Bio::SeqFeature(-start=>$hstart,
				   -end=>$hend,
				   -strand=>1,
				   );
      }
      $homol->homol_SeqFeature($hsf);

      #my $strand;
      #if ($hstart > $hend) {
	#$strand = -1;
      #} else {
	#$strand = 1;
      #}
      #my $homol = Spangle::Homol->new(-start  => $start,
	#			      -end    => $end,
	#			      -strand => $strand,
	#			      -id     => $hid,
	#			      -hstart => $hstart,
	#			      -hend   => $hend,
	#			      -score  => $score);
      #$homol->title($title);

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

1;

