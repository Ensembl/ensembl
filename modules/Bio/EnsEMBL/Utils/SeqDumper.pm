use strict;

package Bio::EnsEMBL::Utils::SeqDumper;

use IO::File;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Root);

#keys must be uppercase
my $DUMP_HANDLERS = 
  { 'FASTA'     => \&dump_fasta,
    'EMBL'      => \&dump_embl,
    'GENBANK'   => \&dump_genbank };

my @COMMENTS = 
  ('This sequence was reannotated via the Ensembl system. Please visit ' .
   'the Ensembl web site, http://www.ensembl.org/ for more information.',

   'All feature locations are relative to the first (5\') base ' .
   'of the sequence in this file.  The sequence presented is '.
   'always the forward strand of the assembly. Features ' .
   'that lie outside of the sequence contained in this file ' .
   'have clonal location coordinates in the format: ' .
   '<clone accession>.<version>:<start>..<end>',

   'The /gene indicates a unique id for a gene, /cds a unique id for a ' .
   'translation and a /exon a unique id for an exon. These ids are ' .
   'maintained wherever possible between versions.',

   'All the exons and transcripts in Ensembl are confirmed by ' .
   'similarity to either protein or cDNA sequences.');


=head2 new

  Arg [1]    : none
  Example    : $seq_dumper = Bio::EnsEMBL::Utils::SeqDumper->new;
  Description: Creates a new SeqDumper 
  Returntype : Bio::EnsEMBL::Utils::SeqDumper
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ($caller, $slice) = @_;

  my $class = ref($caller) || $caller;

  my $feature_types = {'gene'        => 1,
		       'genscan'     => 1,
		       'repeat'      => 1,
		       'similarity'  => 1,
		       'variation'   => 1,
		       'contig'      => 1};

  my $self = bless {'feature_types' => $feature_types}, $class;

  return $self;
}



=head2 enable_feature_type

  Arg [1]    : string $type
  Example    : $seq_dumper->enable_feature_type('similarity');
  Description: Enables the dumping of a specific type of feature
  Returntype : none
  Exceptions : warn if invalid feature type is passed,
               thrown if no feature type is passed
  Caller     : general

=cut

sub enable_feature_type {
  my ($self, $type) = @_;

  $type || $self->throw("type arg is required");

  if(exists($self->{'feature_types'}->{$type})) {
    $self->{'feature_types'}->{$type} = 1;
  } else {
    $self->warn("unknown feature type '$type'\n" .
	  "valid types are: " . join(',', keys %{$self->{'feature_types'}})); 
  }
}



=head2 disable_feature_type

  Arg [1]    : string $type
  Example    : $seq_dumper->disable_feature_type('genes');
  Description: Disables the dumping of a specific type of feature
  Returntype : none
  Exceptions : warn if an invalid feature type is passed,
               thrown if no feature type is passed
  Caller     : general

=cut

sub disable_feature_type {
  my ($self, $type) = @_;
  
  $type || $self->throw("type arg is required");

  if(exists($self->{'feature_types'}->{$type})) {
    $self->{'feature_types'}->{$type} = 0;
  } else {
    $self->warn("unknown feature type '$type'\n" .
	    "valid types are: " . join(',', keys %{$self->{'feature_types'}}));
  }
}



=head2 is_enabled

  Arg [1]    : string $type 
  Example    : do_something() if($seq_dumper->is_enabled('gene'));
  Description: checks if a specific feature type is enabled
  Returntype : none
  Exceptions : warning if invalid type is passed, 
               thrown if no type is passed 
  Caller     : general

=cut

sub is_enabled {
  my ($self, $type) = @_;

  $type || $self->throw("type arg is required");

  if(exists($self->{'feature_types'}->{$type})) {
    return $self->{'feature_types'}->{$type};
  } else {
    $self->warn("unknown feature type '$type'\n" .
	   "valid types are: " . join(',', keys %{$self->{'feature_types'}}));
  }
}


=head2 dump

  Arg [1]    : Bio::EnsEMBL::Slice slice
               The slice to dump
  Arg [1]    : string $format
               The name of the format to dump
  Arg [2]    : (optional) $outfile
               The name of the file to dump to. If no file is specified STDOUT
               is used
  Example    : $seq_dumper->dump($slice, 'EMBL');
  Description: Dumps a region of a genome specified by the slice argument into
               an outfile of the format $format
  Returntype : none
  Exceptions : thrown if slice or format args are not supplied
  Caller     : general

=cut


sub dump {
  my ($self, $slice, $format, $outfile) = @_;

  $format || $self->throw("format arg is required");
  $slice  || $self->throw("slice arg is required");

  my $dump_handler = $DUMP_HANDLERS->{uc($format)};

  unless($dump_handler) {
    $self->throw("No dump handler is defined for format $format\n");
  }


  my $FH = IO::File->new;;
  if($outfile) {
    $FH->open(">$outfile") or $self->throw("Could not open file $outfile");
  } else {
    $FH = \*STDOUT;
    #mod_perl did not like the following
    #$FH->fdopen(fileno(STDOUT), "w") 
    #  or $self->throw("Could not open currently selected output filehandle " .
    #		      "for writing");
  }

  
  &$dump_handler($self, $slice, $FH);

  $FH->close if ($outfile); #close if we were writing to a file
}



=head2 dump_embl

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub dump_embl {
  my $self = shift;
  my $slice = shift;
  my $FH   = shift;

  my $id = $slice->name;
  my $len = $slice->length;


  #line breaks are allowed near the end of the line on ' ', "\t", "\n", ',' 
  $: = (" \t\n-,");

  #############
  # dump header
  #############

  my $EMBL_HEADER = 
'@<   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~
';
  
  #ID and moltype
  my $VALUE = "$id    ENSEMBL; DNA; PLN; $len BP.";
  $self->write($FH, $EMBL_HEADER, 'ID', $VALUE);  
  print $FH "XX\n";
  
  #Accession
  $self->write($FH, $EMBL_HEADER, 'AC', $id);
  print $FH "XX\n";
  
  #Version
  $self->write($FH, $EMBL_HEADER, 'SV', "$id.ENSEMBL_DB:".
	       $slice->adaptor->db->dbname);
  print $FH "XX\n";

  #Date
  $self->write($FH, $EMBL_HEADER, 'DT', $self->_date_string);
  print $FH "XX\n";

  #Description
  $self->write($FH, $EMBL_HEADER, 'DE', "Reannotated sequence via EnsEMBL");
  print $FH "XX\n";

  #key words
  $self->write($FH, $EMBL_HEADER, 'KW', '.');
  print $FH "XX\n";

  #Species
  my $species   = $slice->adaptor->db->get_MetaContainer->get_Species();
  my $species_name = $species->binomial();
  if(my $cn = $species->common_name()) {
    $species_name .= " ($cn)";
  }

  $self->write($FH, $EMBL_HEADER, 'OS', $species_name);

  #Classification
  my @cls = $species->classification;
  shift @cls; #shift off species name
  $self->write($FH, $EMBL_HEADER, 'OC', join('; ', reverse(@cls)) . '.');
  print $FH "XX\n";
  
  #References (we are not dumping refereneces)

  #Database References (we are not dumping these)

  #comments
  foreach my $comment (@COMMENTS) {
    $self->write($FH, $EMBL_HEADER, 'CC', $comment);
    print $FH "XX\n";
  }

  ####################
  #DUMP FEATURE TABLE
  ####################
  print $FH "FH   Key             Location/Qualifiers\n";

  my $FEATURE_TABLE = 
'FT   ^<<<<<<<<<<<<<<<^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~
';
  $self->_dump_feature_table($slice, $FH, $FEATURE_TABLE);  

  #write an XX after the feature tables
  print $FH "XX\n";

  ###################
  #DUMP SEQUENCE
  ###################

  my $SEQ     = $slice->seq();
  my $length  = length($SEQ);
  my $a_count = $SEQ =~ tr/aA/aA/;
  my $c_count = $SEQ =~ tr/cC/cC/;
  my $t_count = $SEQ =~ tr/tT/tT/;
  my $g_count = $SEQ =~ tr/gG/gG/;
  my $other_count = $length - $a_count - $c_count - $t_count - $g_count;

  my $value = "Sequence $length BP; $a_count A; $c_count C; " .
    "$g_count G; $t_count T; $other_count other;";
  $self->write($FH, $EMBL_HEADER, 'SQ', $value);

  $self->write_embl_seq($FH, \$SEQ);

  print $FH "//\n";

  # Set formatting back to normal
  $: = " \n-";
}




=head2 dump_genbank

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub dump_genbank {
  my ($self, $slice, $FH) = @_;

  #line breaks are allowed near the end of the line on ' ', "\t", "\n", ',' 
  $: = " \t\n-,";

  my $GENBANK_HEADER = 
'^<<<<<<<<<  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
';

  my $GENBANK_SUBHEADER =
'  ^<<<<<<<  ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
';

  my $GENBANK_FT =
'     ^<<<<<<<<<<<<<< ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<~~
';

  my $id = $slice->name;
  my $length = $slice->length;

  my ($name_str, $start, $end);
  if($slice->isa('Bio::EnsEMBL::Slice')) {
    $name_str  = "chromosome " . $slice->chr_name;
    $start = $slice->chr_start;
    $end  = $slice->chr_end;
  } else {
    $name_str = $slice->name;
    $start = 1;
    $end   = $slice->length();
  }

  my $date = $self->_date_string;

  my $species = $slice->adaptor->db->get_MetaContainer->get_Species;

  #LOCUS
  my $tag   = 'LOCUS';
  my $value = "ENS:$id $length bp DNA HTG $date";
  $self->write($FH, $GENBANK_HEADER, $tag, $value);


  #DEFINITION
  $tag   = "DEFINITION";
  $value = $species->binomial . ' ' . $slice->adaptor->db->assembly_type . 
    " assembly reannotated via EnsEMBL DNA, $name_str $start..$end";
  $self->write($FH, $GENBANK_HEADER, $tag, $value);

  #ACCESSION
  $self->write($FH, $GENBANK_HEADER, 'ACCESSION', $id);

  #VERSION
  $self->write($FH, $GENBANK_HEADER, 'VERSION', "$id.ENSEMBL_DB:".
	       $slice->adaptor->db->dbname);

  # KEYWORDS
  $self->write($FH, $GENBANK_HEADER, 'KEYWORDS', '.');

  # SOURCE
  $self->write($FH, $GENBANK_HEADER, 'SOURCE', $species->common_name());

  #organism
  my @cls = $species->classification();
  shift @cls;
  $self->write($FH, $GENBANK_SUBHEADER, 'ORGANISM', $species->binomial);
  $self->write($FH, $GENBANK_SUBHEADER, '', join('; ', reverse @cls) . ".");

  #refereneces

  #comments
  foreach my $comment (@COMMENTS) {
    $self->write($FH, $GENBANK_HEADER, 'COMMENT', $comment);
  }
  
  
  ####################
  # DUMP FEATURE TABLE
  ####################
  print $FH "FEATURES             Location/Qualifiers\n";
  $self->_dump_feature_table($slice, $FH, $GENBANK_FT);


  ####################
  # DUMP SEQUENCE
  ####################

  my $SEQ       = $slice->seq();
  my $a_count = $SEQ =~ tr/aA/aA/;
  my $c_count = $SEQ =~ tr/cC/cC/;
  my $t_count = $SEQ =~ tr/tT/tT/;
  my $g_count = $SEQ =~ tr/gG/gG/;
  my $length = length($SEQ);
  my $other_count = $length - $a_count - $c_count - $t_count - $g_count;

  $tag   = 'BASE COUNT';
  $value = "$a_count a $c_count c $g_count g $t_count t";
  $value .= " $other_count n" if($other_count);
  $self->write($FH, $GENBANK_HEADER, $tag, $value);
  print $FH "ORIGIN\n";

  $self->write_genbank_seq($FH, \$SEQ);

  print $FH "//\n";

  # Set formatting back to normal
  $: = " \n-";
}



=head2 _dump_feature_table

  Arg [1]    : Bio::EnsEMBL::Slice slice
  Example    : none
  Description: Helper method used to dump feature tables used in EMBL, FASTA,
               GENBANK.  Assumes formating of file handle has been setup
               already to use $FEAT and $VALUE values.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub _dump_feature_table {
  my $self   = shift;
  my $slice  = shift;
  my $FH     = shift;
  my $FORMAT = shift;

  #use only the core database to dump features (except for bloody snps)
  my $lite = $slice->adaptor->db->remove_db_adaptor('lite');

  my $meta = $slice->adaptor->db->get_MetaContainer;
  my $species = $meta->get_Species;

  #lump file handle and format string together for simpler method calls
  my @ff = ($FH, $FORMAT);
  my $value;

  #source
  my $classification = join(', ', $species->classification);
  $self->write(@ff,'source', "1.." . $slice->length());
  $self->write(@ff,''      , '/classification="'.$classification.'"');
  $self->write(@ff,''      , '/organism="'.$species->binomial . '"');
  $self->write(@ff,''      , '/db_xref="taxon:'.$meta->get_taxonomy_id().'"');

  #
  # Transcripts & Genes
  #
  if($self->is_enabled('gene') && $slice->can('get_all_Genes')) {

    foreach my $gene (@{$slice->get_all_Genes}) {
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
	my $translation = $transcript->translation;
	$value = $self->features2location($transcript->get_all_Exons);
	$self->write(@ff,'CDS', $value);
	$self->write(@ff,''   , '/gene="'.$gene->stable_id().'"');
	$self->write(@ff,''   , '/cds="'.$translation->stable_id().'"');
	$self->write(@ff,''   , '/transcript="'.$transcript->stable_id().'"');

	foreach my $dbl (@{$transcript->get_all_DBLinks}) {
	  $value = '/db_xref="'.$dbl->dbname().':'.$dbl->primary_id().'"';
	  $self->write(@ff, '', $value);
	}
	$value = '/translation="'.$transcript->translate()->seq().'"';
	$self->write(@ff, '', $value);
      }
    }
 
    # exons
    foreach my $gene (@{$slice->get_all_Genes}) {
      foreach my $exon (@{$gene->get_all_Exons}) {
	$self->write(@ff,'exon', $self->features2location([$exon]));
	$self->write(@ff,''    , '/exon_id="'.$exon->stable_id().'"');
      }
    }
  }

  #
  # genscans
  #
  if($self->is_enabled('genscan')) {
    my @genscan_exons;
    foreach my $transcript (@{$slice->get_all_PredictionTranscripts}) {
      my $exons = $transcript->get_all_Exons();
      push @genscan_exons, @$exons;
      $self->write(@ff, 'CDS_gscan', $self->features2location($exons));
      $self->write(@ff, '', '/transcript="'.$transcript->stable_id().'"');
      $self->write(@ff, '', 
		   '/translation="'.$transcript->translate()->seq().'"');
    }
    foreach my $exon (@genscan_exons) {
      $self->write(@ff, 'exon_gscan', $self->features2location([$exon]));
      $self->write(@ff, ''          , '/start_phase="'.$exon->phase.'"'); 
      $self->write(@ff, ''          , '/end_phase="'.$exon->end_phase.'"');
    }
  }

  #
  # snps
  #
  if($self->is_enabled('variation') && $slice->can('get_all_SNPs')) {
    $slice->adaptor->db->add_db_adaptor('lite', $lite) if $lite;

    foreach my $snp (@{$slice->get_all_SNPs}) {
      my $ss = $snp->start;
      my $se = $snp->end;
      #skip snps that hang off edge of slice
      next if($ss < 1 || $se > $slice->length); 

      $self->write(@ff, 'variation', "$ss..$se");
      $self->write(@ff, ''         , '/replace="'.$snp->alleles.'"'); 
      #$self->write(@ff, ''         , '/evidence="'.$snp->status.'"'); 
      foreach my $link ($snp->each_DBLink) {
	my $id = $link->primary_id;
	my $db = $link->database;
	$self->write(@ff, '', "/db_xref=\"$db:$id\""); 
      }
    }

    $slice->adaptor->db->remove_db_adaptor('lite') if $lite;
  }

  #
  # similarity features
  #     
  if($self->is_enabled('similarity')) {
    foreach my $sim (@{$slice->get_all_SimilarityFeatures}) {
      $self->write(@ff, 'similarity', $self->features2location([$sim]));
      $self->write(@ff, ''       , '/hit_name="'.$sim->hseqname.'"');
      $self->write(@ff, ''       , '/hit_score="'.$sim->score.'"');
      $self->write(@ff, ''       , '/percent_ident="'.$sim->percent_id.'"');
      $self->write(@ff, ''       , '/hit_start="'.$sim->hstart.'"');
      $self->write(@ff, ''       , '/hit_end="'.$sim->hend.'"');
      $self->write(@ff, ''       , '/hit_strand="'.$sim->hstrand.'"');
    }
  }

  #
  # repeats
  #
  if($self->is_enabled('repeat')) {
    foreach my $repeat (@{$slice->get_all_RepeatFeatures}) {
      $self->write(@ff, 'repeat', $self->features2location([$repeat]));
      $self->write(@ff, ''    , '/name="'.$repeat->repeat_consensus->name.'"');
    }
  }

  #
  # contigs
  #
  if($self->is_enabled('contig')) {
    foreach my $tile (@{$slice->get_tiling_path}) {
      $self->write(@ff, 'contig', 
		   $tile->assembled_start .'..'. $tile->assembled_end);
      $self->write(@ff, '', '/name="'.$tile->component_Seq->name.'"');
      $self->write(@ff, '', '/contig_start="'.$tile->component_start.'"');
      $self->write(@ff, '', '/contig_end="'.$tile->component_end.'"');
      $self->write(@ff, '', '/contig_orientation="'.$tile->component_ori.'"');
    }
  }

  $slice->adaptor->db->add_db_adaptor('lite', $lite) if $lite;

}



=head2 dump_fasta

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub dump_fasta {
  my $self = shift;
  my $slice = shift;
  my $FH   = shift;


  my $species = 
    $slice->adaptor->db->get_MetaContainer->get_Species->binomial();
  
  my $name = $slice->name;
  my $start = 1;
  my $end = $slice->length();

  my $header = ">$species|$name\n";
  print $FH $header;

  #set the formatting to FASTA
  my $FORMAT = '^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
';

  #chunk the sequence in 60kb chunks to use less memory
  my $cur = $start;
  while($cur < $end) {
    my $to = $cur + 59_999;
    $to = $end if($to > $end); 
    my $seq = $slice->subseq($cur, $to);
    $cur = $to + 1;
    $self->write($FH, $FORMAT, $seq);
  }
}



=head2 features2location

  Arg [1]    : listref of Bio::EnsEMBL::SeqFeatures
  Example    : $location = $self->features2location(\@features);
  Description: Constructs an EMBL location string from a list of features
  Returntype : string
  Exceptions : none
  Caller     : internal

=cut

sub features2location {
  my $self = shift;
  my $features = shift;
  
  my @join = ();

  foreach my $f (@$features) {
    my $ctg = $f->contig;

    if($ctg->isa('Bio::EnsEMBL::Slice')) {
      if($f->start >= 1 && $f->end <= $ctg->length) {
	#this feature in on a slice and doesn't lie outside the boundary
	
	if($f->strand() == 1) {
	  push @join, $f->start()."..".$f->end();
	} else {
	  push @join, "complement(".$f->start()."..".$f->end().")";
	}
      } else {
	my @fs = ();
	#this feature is outside the boundary of the dump, 
	#convert to clone coords
	if($ctg->isa('Bio::EnsEMBL::Slice')) {
	  #copy the feature, and transform to contig coords
	my $new_f;
	%$new_f = %$f;
	bless $new_f, ref $f;
	push @fs, $new_f->transform;
      }
	
	if($fs[0]->isa('Bio::EnsEMBL::StickyExon')) {
	  @fs = @{$fs[0]->get_all_component_Exons};
	}
	
	#use the accession:x-y format
	foreach my $feat (@fs) {
	  my $contig = $feat->contig;
	  my $clone = $contig->clone;
	  my $acc;
	  if($clone->embl_id) {
	    $acc = $clone->embl_id .'.'. $clone->embl_version;
	  } else {
	    $acc = $clone->id;
	  }
	  my $start = $feat->start + $contig->embl_offset - 1;
	  my $end   = $feat->end + $contig->embl_offset - 1;
	  my $strand = $feat->strand;
	  if($strand == 1) {
	    push @join, "$acc:$start..$end";
	  } else {
	    push @join, "complement($acc:$start..$end)";
	  }
	}
      }
    } else {
      #handle features in contig coordinates

      my $start = $f->start;
      my $end = $f->end;

      if($start < 1 || $end > $ctg->length) {
	$self->throw("Feature off end of contig boundary");
      }

      if($f->strand == -1) {
	push @join, "complement($start..$end)";
      } else {
	push @join, "($start..$end)";
      }
    }
  }
    
  my $out = join ',', @join;

  if(scalar @join > 1) {
    $out = "join($out)";
  }

  return $out;
}


sub _date_string {
  my $self = shift;

  my ($sec, $min, $hour, $mday,$mon, $year) = localtime(time());

  my $month = ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
	       'AUG', 'SEP', 'OCT', 'NOV', 'DEC')[$mon];
  $year += 1900;
  
  return "$mday-$month-$year";
}


sub write {
  my ($self, $FH, $FORMAT, @values) = @_;
  
  #while the last value still contains something
  while($values[-1] ne '') {
    formline($FORMAT, @values);
    print $FH "$^A";
    $^A = '';
  }
}

sub write_genbank_seq {
  my $self = shift;
  my $FH  = shift;
  my $seq = shift;
  my $base_total = shift;

  $base_total ||= 0;

  my $GENBANK_SEQ = 
'@>>>>>>>> ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<<~
';

  my $total = -59 + $base_total;
  #keep track of total and print lines of 60 bases with spaces every 10bp
  while($$seq) {
    $total += 60; 
    formline($GENBANK_SEQ,$total, $$seq, $$seq, $$seq, $$seq, $$seq, $$seq);
    print $FH $^A;
    $^A = '';
  }
}

sub write_embl_seq {
  my $self = shift;
  my $FH   = shift;
  my $seq  = shift;
  my $base_total = shift;

  $base_total ||= 0;

  my $EMBL_SEQ = 
'     ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<< ^<<<<<<<<<@>>>>>>>>>~
';
  #keep track of total and print lines of 60 bases with spaces every 10bp
  my $length = length($$seq);
  my $total = $length - $base_total;
  while($$seq) {
    $total -= 60;
    $total = 0 if($total < 0);
    formline($EMBL_SEQ, 
	     $$seq, $$seq, $$seq, $$seq, $$seq, $$seq, 
	     $length - $total);
    print $FH $^A;
    $^A = '';
  }
}

1;   
