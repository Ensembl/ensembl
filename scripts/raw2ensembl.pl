# raw2ensembl
#
# Cared for by Simon Potter
# (C) GRL/EBI 2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod

=head1 NAME

raw2ensembl.pl

=head1 SYNOPSIS

raw2ensembl.pl -clone AB000381.1 -fasta <file> -contigs <file>
<db options> -write -replace

=head1 DESCRIPTION

Load clone into EnsEMBL database from raw sequence. Takes a fasta file,
splits into contigs - either from coords in a file or by looking for
runs of 'n' in the sequence - and loads in to the DB

=head1 OPTIONS

    -dbhost  DB host
    -dbuser  DB user
    -dbname  DB name
    -contigs file of contig start/end pairs
    -phase   clone phase
    -clone   clone accession and version (separated by '.')
    -fasta   fasta file of clone sequence (fasta header must contain accession)
    -single  only one seq in fasta file (no need for header to match)
    -write   write clone
    -replace replace existing clone
    -v       print info about clones and contigs

=head1 CONTACT

Simon Potter: scp@sanger.ac.uk

=head1 BUGS

Insert list of bugs here!

=cut


use strict;
use Getopt::Long;
use Bio::Root::RootI;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::PerlDB::Clone;
use Bio::EnsEMBL::PerlDB::Contig;


my($id, $acc, $ver, $phase, $contigs);
my($fasta, $single, $seqio, $seq);
my($dbname, $dbhost, $dbuser);
my($help, $info, $write, $replace, $verbose);


$Getopt::Long::autoabbrev = 0;   # personal preference :)
$dbuser = 'ensadmin';            # default

my $ok = &GetOptions(
    "clone=s"   => \$id,
    "phase=s"   => \$phase,
    "fasta=s"   => \$fasta,
    "single"    => \$single,
    "contigs=s" => \$contigs,
    "dbname=s"  => \$dbname,
    "dbhost=s"  => \$dbhost,
    "dbuser=s"  => \$dbuser,
    "help"      => \$help,
    "info"      => \$info,
    "write"     => \$write,
    "replace"   => \$replace,
    "v"         => \$verbose
);

if ($help || not $ok) {
    &usage;
    exit 0;
} elsif ($info) {
    exec("perldoc $0");
}

if (defined $id && $id =~ /(\S+)\.(\d+)/) {
    ($acc, $ver) = ($1, $2);
}
else {
    print STDERR "Must specify -clone <clone.version>\n";
    exit 1;
}

unless ($dbname && $dbuser && $dbhost) {
    print STDERR "Must specify all DB parameters\n";
    exit 1;
}

unless ($fasta) {
    print STDERR "Must specify -fasta\n";
    exit 1;
}

if ($phase < 0 && $phase > 4) {
    print STDERR "Phase should be 1, 2, 3 or 4\n";
    exit 1;
}
$phase = -1 unless defined $phase;

if ($phase == 4 and defined $contigs) {
    print STDERR "Don't need contig info for phase 4 clone - ignoring\n";
    undef $contigs;
}


my $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    '-host'   => $dbhost,
    '-user'   => $dbuser,
    '-dbname' => $dbname
) or die "Can't connect to DB $dbname on $dbhost as $dbuser";


open DNA, "< $fasta" or die "Can't open file $fasta for DNA";
$seqio = new Bio::SeqIO(
    -fh     => \*DNA,
    -format => 'fasta'
);
while ($seq = $seqio->next_seq) {
    last if $single or $seq->display_id =~ /$id/;
    undef $seq;
}
close DNA;
die "Couldn't get $id from $fasta" unless defined $seq;

print STDERR "display id  ", $seq->display_id, "\n";
print STDERR "description ", $seq->desc, "\n";

# my $seq = $seqio->next_seq or die "Canna get seq from file $fasta";

print "Loaded dna: length ", $seq->length, "\n";

my $clone = new Bio::EnsEMBL::PerlDB::Clone;
$clone->id($acc);
$clone->htg_phase($phase);
$clone->embl_id($acc);
$clone->version(1);
$clone->embl_version($ver);

print "Clone ", $clone->id, "\n";
if ($verbose) {
    print "\tembl_id     ", $clone->embl_id, "\n";
    print "\tversion     ", $clone->version, "\n";
    print "\temblversion ", $clone->embl_version, "\n";
    print "\thtg_phase   ", $clone->htg_phase, "\n";
}

if ($phase == 4) {
    my $contig = new Bio::EnsEMBL::PerlDB::Contig;
    my $length = $seq->length;
    $contig->id("$acc.$ver.1.$length");
    $contig->embl_offset(1);
    $contig->length($length);
    $contig->seq($seq);
    $contig->version(1);
    $contig->embl_order(1);

    print "Contig ", $contig->id, "\n";
    if ($verbose) {
	print "\toffset: ", $contig->embl_offset, "\n";
	print "\tlength: ", $contig->length, "\n";
	print "\tend:    ", ($contig->embl_offset + $contig->length - 1), "\n";
	print "\tversion:", $contig->version, "\n";
	print "\torder:  ", $contig->embl_order, "\n";
	print "\tlength: ", $contig->length, "\n";
    }

    $clone->add_Contig($contig);
}
else {
    my @split;
    if (defined $contigs) {
	@split = &getContigs($contigs);
    }
    else {
        @split = &scanClone($seq->seq);
    }
    my $order = 1;
    foreach my $startend (@split) {
	my $offset = $startend->[0];
	my $length = $startend->[1] - $startend->[0] + 1;
	my $id = join '.', ($acc, $ver, $offset, $startend->[1]);
	my $subseq = $seq->subseq($offset, $startend->[1]);

	my $contig = new Bio::EnsEMBL::PerlDB::Contig;
	$contig->id($id);
	$contig->embl_offset($offset);
	$contig->length($length);
	$contig->seq(new Bio::Seq(-id => $id, -seq =>$subseq));
	$contig->version(1);
	$contig->embl_order($order);
	$order++;

	print "Contig ", $contig->id, "\n";
	if ($verbose) {
	    print "\toffset  ", $contig->embl_offset, "\n";
	    print "\tlength  ", $contig->length, "\n";
	    print "\tend     ", ($contig->embl_offset + $contig->length - 1), "\n";
	    print "\tversion ", $contig->version, "\n";
	    print "\torder   ", $contig->embl_order, "\n";
	}

        $clone->add_Contig($contig);
    }
}

if ($write) {
    my $dbclone;
    eval {
	$dbclone = $dbobj->get_Clone("$acc.$ver");
    };
    if ($replace && $dbclone) {
	$dbclone->delete;
	$dbobj->write_Clone($clone);
    }
    elsif ($dbclone) {
	print "$acc.$ver already exists - ignoring\n";
    }
    else {
	$dbobj->write_Clone($clone);
    }
}



sub usage {
    print <<EOF
$0 [options]
Options:
  -clone    accession.version of clone to load
  -dbname
  -dbhost
  -dbuser
  -fasta    fasta file of clone sequence
  -single   only read first seq in file
  -phase    clone phase
  -contigs  file of contig start/end pairs
  -clone    clone accession and version (separated by '.')
  -write    write clone
  -replace  replace existing clone
  -v        print info about clones and contigs
EOF
}



=head2 scanClone

    Title   :   scanClone
    Usage   :   @contigs = $obj->scanClone($seq)
    Function:   Scans the clone sequence to find positions of contigs
		by assuming at least x n's  between contigs
    Returns :   list of lists (start, end)
    Args    :   string

=cut

sub scanClone {
  my($seq) = @_;
  my(@gaps, @contig, $start, $gap);

  # get a list of gaps - at least 50 bp
  my $pos = 0;
  while ($pos < length $seq) {
    my $unused = substr $seq, $pos;
    ($gap) = $unused =~ /(n{50,})/;
    last unless $gap;
    $start = 1 + index $seq, $gap, $pos;
    push @gaps, [ $start, $start + length($gap) - 1 ];
    $pos = $start + length $gap;
  }

  # calc coords of contigs

  if (@gaps){
    # 1st contig before 1st gap unless the sequence starts off with a gap
    push @contig, [1, $gaps[0]->[0] - 1] unless $gaps[0]->[0] == 1;

    # contigs other than 1st and last are between gaps
    foreach my $i (0 .. $#gaps - 1) {
      push @contig, [$gaps[$i]->[1] + 1, $gaps[$i + 1]->[0] - 1];
    }

    # last contig after last gap unless the sequence ends with a gap
    push @contig, [$gaps[$#gaps]->[1] + 1, length($seq)]
     unless $gaps[$#gaps]->[1] == length($seq);
  }
  else {
    # no gaps
    push @contig, [1, length($seq)];
  }

  return sort {$a->[0] <=> $b->[0]} @contig;
}



=head2 getContigs

    Title   :   getContigs
    Usage   :   @contigs = $obj->getContigs($seq)
    Function:   Reads contig info from file
    Returns :   list of lists (start, end)
    Args    :   string

=cut

sub getContigs {
    my($file) = @_;
    my($start, $end, @contigs);

    open CONTIG, "< $file" or die "Can't open contigs file $file";
    while (<CONTIG>) {
	chomp;
	($start, $end) = split;
	unless ($start && $end && $start < $end) {
	    die "Illegal line in contig file: $_\n";
	}
	push @contigs, [ $start, $end ];
    }
    close CONTIG;

    return @contigs;
}
