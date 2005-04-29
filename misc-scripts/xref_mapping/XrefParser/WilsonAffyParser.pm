package XrefParser::WilsonAffyParser;

use strict;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);
my $xref_sth ;
my $dep_sth;
my $syn_sth;



sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs = $self->create_xrefs($source_id, $species_id, $file);

  # upload
  XrefParser::BaseParser->upload_xref_object_graphs(@xrefs);

}

sub create_xrefs {

  my ($self, $source_id, $species_id, $file) = @_;

  my ($count, $noseq, $direct) = (0,0,0);

  $| = 1; # don't buffer

  my @xrefs;

  open(FILE,"<".$file) || die "Could not open $file\n";

  <FILE>; # skip first line

  while (<FILE>) {

    #last if ($count > 200);
    my $xref;

    my @fields = split /\t/;

    # first field (probe_set) is accession
    my $acc = $fields[0];
    chomp($acc);
    $acc =~ s/\"//g;


    # get linked accession (may be RefSeq or EMBL or ensembl)
    my $target = $fields[2];
    chomp($target);
    $target =~ s/\"//g;

    # Create direct xrefs for mappings to Ensembl transcripts
    if ($target =~ /ENSGALT/) {

      # remove version if present
      ($target) = $target =~ /([^.]*)\.([^.]*)/;

      # add xref - not we're assuming it doesn't already exist;
      # may need to check like in CCDS parser
      my $xref_id = $self->add_xref($acc, 0, $acc, "", $source_id, $species_id);
      $self->add_direct_xref($xref_id, $target, "transcript", "");
      $direct++;

    } else {

      # fetch sequence for others (EMBL ESTs and RefSeqs - pfetch will handle these)
      system ("pfetch -q $target > seq.txt");
      open(SEQ, "<seq.txt");
      my $seq = <SEQ>;
      chomp($seq);
      close(SEQ);

      if ($seq && $seq !~ /no match/) {

	$xref->{ACCESSION} = $acc;
	$xref->{SEQUENCE} = $seq;
	$xref->{LABEL} = $acc;
	$xref->{SOURCE_ID} = $source_id;
	$xref->{SPECIES_ID} = $species_id;
	$xref->{SEQUENCE_TYPE} = 'dna';
	$xref->{STATUS} = 'experimental';

	# Add description noting where the mapping came from
	$xref->{DESCRIPTION} = $target . " used as mapping target";

	#print $xref->{ACCESSION} . " " . $target . " " . $? . "\n";

	$count++;

	print "$count " if ($count % 100 == 0);

	push @xrefs, $xref;

      } else {

	print "Couldn't get sequence for $target\n";
	$noseq++;

      }

    }

  }

  close(FILE);

  print "\n\nParsed $count primary xrefs.\n";
  print "Couldn't get sequence for $noseq primary_xrefs\n" if ($noseq);
  print "Added $direct direct xrefs.\n";

  return \@xrefs;

}

sub new {

  my $self = {};
  bless $self, "XrefParser::WilsonAffyParser";
  return $self;

}

1;
