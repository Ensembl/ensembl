#!/usr/local/ensembl/bin/perl -w

# Written by Jan-Hinnerk Vogel, modifed by mh4
#
# Copyright GRL/EBI 2004
# You may distribute this script under the same terms as perl itself
#
## Original design:
## Script reads a large fasta-file containg the sequence of a whole chromosome
## and splits this file into smaller chunks ($chunk_size). The chunks are
## written in a new file FILE.fasta.splitted. Additionally an AGP-file 
## is written (containing the chr-name). The sequence can be loaded in db 
## using load_seq_region.pl and load_agp.pl.
## - script was used in droso-build DROM3B

# Modified by mh4 to split Anopheles scaffolds into chunks
# Used for assembly AgamP3
# Takes a single fasta file of all the scaffolds
# produces a single fasta file of all the chunks and a single agp file describing the chunking
# default is split into equally-sized ~50kb chunks if scaffold is >= 100kb
# Chunking part adapted from misc_scripts/anopheles_scripts load_dna.pl
# Command line needs input file name in form -file name.fasta
# Writes a single multifasta chunks file and a single agp file to location of input file given in the command line


use strict;
use Bio::SeqIO;
use Getopt::Long;


my $file;

my $format = 'FASTA';
my $outformat = 'FASTA';


my $chunk_size = 50_000;

my $type = "chunk";


&GetOptions(
            'file:s'   => \$file,
            'size:i'  => \$chunk_size,
            'type:s'   => \$type,
           ) ;



(my $basename = $file) =~ s/\.fasta//gi;


my $outfile = $basename."_chunks.fasta";
my $agpfile = $basename.".agp";

if( !$file){
  print "You have to supply name of fasta-file like this:\n\tperl_script.pl  -file scaffseqs.fasta\n\n";
  exit(0);
}

my $inseq = Bio::SeqIO->new(
                             '-file'   => "<$file",
                             '-format' => $format,
                            );


my $outseq = Bio::SeqIO->new(
                             -file     => ">>$outfile", # changed to appending
                             -format => $outformat,
                            );

my $chunks_total = 0;
my $scaff_count=0;
my @agp_file_data;

# process each sequence-object
while ( my $seq_obj = $inseq->next_seq() ) {
  $scaff_count ++;
  my $scaff_l = $seq_obj->length();
  my $ac = $seq_obj->id;
  print STDERR "Scaffold $ac has length: $scaff_l\n";

#set length of chunks
  my $div;
  if ($scaff_l < $chunk_size) {
    $div = 1;
  }
  else {
    $div = int ($scaff_l/$chunk_size);
  }
  my $chunk_l = int ($scaff_l/$div);
  print STDERR "AC: $ac\tDIV: $div\tCHUNK_L: $chunk_l\n";

#make chunks
  my $chunk_count = 1;
  my $next_start;
  my $total_length=0;

  while ($chunk_count <= $div) {

  #first chunk
    if ($chunk_count == 1) {
      my $chunk_id = "$ac".'_'."$chunk_count";
      my $subseq = $seq_obj->subseq(1,$chunk_l);
      my $chunk_obj = Bio::PrimarySeq->new(
                                           -seq => $subseq,
                                           -id  => $chunk_id,
                                          );
      my $subseq_l = length($subseq);

      # print chunk seq in fasta format
      print STDERR "CHUNK_ID: $chunk_id\tL: $subseq_l\n";
      $outseq->write_seq($chunk_obj);

      # store data for agp file in array
      push @agp_file_data, [$ac,1,$chunk_l,$chunk_count,"chunk",$chunk_id,1,$subseq_l,'1'];

      $total_length += $subseq_l;
      last if $div==1;

      $next_start = $total_length +1;
      $chunk_count++;
    }

  #chunks between first and last
    elsif (($chunk_count > 1) && ($chunk_count < $div)) {
      my $chunk_id = "$ac".'_'."$chunk_count";
      my $end = $next_start + $chunk_l -1;
      my $subseq = $seq_obj->subseq($next_start,$end);
      my $chunk_obj = Bio::PrimarySeq->new(
                                           -seq => $subseq,
                                           -id  => $chunk_id,
                                          );
      my $subseq_l = length($subseq);

      # print chunk seq in fasta format
      print STDERR "CHUNK_ID: $chunk_id\tL: $subseq_l\n";
      $outseq->write_seq($chunk_obj);

      # store data for agp file in array
      push @agp_file_data, [$ac,$next_start,$end,$chunk_count,"chunk",$chunk_id,1,$subseq_l,'1'];

      $total_length += $subseq_l;
      $next_start = $total_length +1;
      $chunk_count++;
    }


  #last chunk
    elsif ($chunk_count == $div) {
      my $chunk_id = "$ac".'_'."$chunk_count";
      my $subseq = $seq_obj->subseq($next_start,$scaff_l);
      my $subseq_l = length($subseq);
      my $chunk_obj = Bio::PrimarySeq->new(
                                           -seq => $subseq,
                                           -id  => $chunk_id,
                                          );

      # print chunk seq in fasta format
      print  "CHUNK_ID: $chunk_id\tL: $subseq_l\n";
      $outseq->write_seq($chunk_obj);

      # store data for agp file in array
      push @agp_file_data, [$ac,$next_start,$scaff_l,$chunk_count,$type,$chunk_id,1,$subseq_l,'1'];

      $total_length += $subseq_l;
      last;
    }

    else {
      print STDERR "DUD RUN - internal counting screwed up\n";
    }
  }
  if ($total_length != $scaff_l) {
    print STDERR  "ERROR:  CALCULATED TOTAL: $total_length does not = REAL TOTAL LENGTH: $scaff_l\n";
    die;
  }

  print "SCAFF $ac\tCHUNKS $chunk_count\tTotal length $total_length\n\n";
  $chunks_total += $chunk_count;
}


# write agp-file

open AGP, ">>$agpfile" || die "write-error $agpfile\n";
my $agp_lines = 0;
foreach my $line (@agp_file_data) {
  my @line = @$line;
  print AGP join("\t",@line)."\n";
  $agp_lines ++;
}
close AGP;


print "$scaff_count scaffolds in file $file split into $chunks_total chunks and agp-file written with $agp_lines lines\n\tNow, import seqs with load_seqregion.pl script and agp with load_agp.pl script\n";


# Example of agp output with a chunk size of 12
#asm_name asm_start asm_end nr type cmp_seq_region_name cmp_start cmp_end cmp_strand

#scaffA  1       12      1       chunk   scaffA_1        1       12      1
#scaffA  13      24      2       chunk   scaffA_2        1       12      1
#scaffA  25      36      3       chunk   scaffA_3        1       12      1
#scaffA  37      49      4       chunk   scaffA_4        1       13      1
#scaffB  1       12      1       chunk   scaffB_1        1       12      1
#scaffB  13      24      2       chunk   scaffB_2        1       12      1
#scaffB  25      36      3       chunk   scaffB_3        1       12      1
#scaffB  37      48      4       chunk   scaffB_4        1       12      1
#scaffB  49      60      5       chunk   scaffB_5        1       12      1
#scaffC  1       23      1       chunk   scaffC_1        1       23      1
