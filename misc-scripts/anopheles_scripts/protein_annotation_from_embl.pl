#!/usr/local/ensembl/bin/perl -w

# script to pull selected protein annotation from an EMBL scaffold entry
# point to file with all 69724 anopheles scaffolds/contigs
# STDOUT links ENSANGP# and EAA# protein id, inter alia

use strict;
use Bio::SeqIO;
use Bio::SeqFeatureI;

my $infile = "/ecs2/work6/mh4/data/wgs_aaab01_embl";
my $stream = new Bio::SeqIO (-file => $infile,
				-format => 'EMBL');
my ($inputcount, $id);
my $scaffoldcount=0;
while (my $seq = $stream->next_seq) {
  $inputcount +=1;
  last if $scaffoldcount > 9000;         #temp limit for testing
  $id=$seq->id;
  if ($id=~/AAAB01(\d{6})/) {
    next if $1 > 8987;
  }
  else {print STDERR "Scaffold $id has an invalid id\n";}
  $scaffoldcount +=1;
  if ($seq -> feature_count > 1) {
    print STDERR "$id\n";
    my @features = $seq -> all_SeqFeatures;
    foreach my $feature (@features) {
#      print STDERR "\t".$feature -> primary_tag."\n";
      if ($feature -> primary_tag eq 'CDS') {
# ugly use of x's to help printing - better to make '.' or deal with this in the printing?
	my $note = 'x';
	my $protein_id = 'x';
	my $old_id = 'x';
	my (@product_values);
	if ($feature->has_tag('product')) {
	  @product_values = $feature -> each_tag_value('product');
# may be CDS with >1 product in EMBL, if aa sequence identical
	  print STDERR scalar(@product_values)." product values! in $id\n" if scalar(@product_values)>1;
	}
	if ($feature->has_tag('protein_id')) {
	  my @protein_values = $feature -> each_tag_value('protein_id');
# should be no CDS features with >1 protein_id
	  print STDERR scalar(@protein_values)." protein values!\n" if scalar(@protein_values)>1;
	  my $protein_id_v = join ' ',@protein_values;
	  $protein_id = $1 if ($protein_id_v =~ /^(\w+)\.\w+$/);
	}
	if ($feature->has_tag('note')) {
	  my @note_values = $feature -> each_tag_value('note');
# should be no CDS features with >1 note
# but this may change next time
	  print STDERR scalar(@note_values)." note values!\n" if scalar(@note_values)>1;
	  $note = join ' ',@note_values;
# parsing of this will change if I print just the old id as a note next time
	  if ($note=~/^gnl\|WGS:AAAB\|(\w+)\|gb\|(\w+)/) {
	    $old_id = $1;
	    my $temp_prot_id = $2;
	    print STDERR "ids don't match\n" if ($protein_id ne $temp_prot_id);
	  }
	  else {print STDERR "Can't parse old id from #$note#\n";}
	}
# print a line for each product (ENSANGPxx id)
# no old id if there is a note but it fails to parse (misses 42)
	foreach my $product_value (@product_values) {
	  print "$id\t$product_value\t$protein_id\t$old_id\t$note\n";
	}
      }
    }
  }
}
print STDERR "\nRead $inputcount seqs and $scaffoldcount selected as scaffolds\n";
