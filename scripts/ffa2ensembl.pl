# ffa2ensembl
#
# Cared for by Simon Potter
# (C) GRL/EBI 2001
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod

=head1 NAME

ffa2ensembl.pl

=head1 SYNOPSIS

ffa2ensembl.pl -clones clones.txt -ffa draft.ffa.gz

=head1 DESCRIPTION

Load clones into EnsEMBL database from an ffa file. Takes a list of
clones and an ffa file and creates clones and contigs. The clones in
the ffa file are already split into contigs. It is asssumed that all
contigs for a particular clone will be together.

=head1 OPTIONS

    -dbhost   DB host
    -dbuser   DB user
    -dbname   DB name
    -clones   file containing list of "clone.version" or literal 'all'
    -clobber  overwrite existing clone
    -ffa      location of ffa.gz file

=head1 CONTACT

Simon Potter: scp@sanger.ac.uk

=head1 BUGS

Insert list of bugs here!

=cut


use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::PerlDB::Contig;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Root::RootI;
use Bio::EnsEMBL::PerlDB::Clone;
use Bio::Seq;
use Bio::SeqIO;


my ($clone, $contig);
my (%required_clones, %written_clone);
my ($clonelist, $ffa, $seqfile, $write, $clobber);
my ($dbname, $dbhost, $dbuser);
my ($help, $info);
my (%inDB);


$Getopt::Long::autoabbrev = 0;   # personal preference :)
$dbuser = 'ensadmin';  # default

&GetOptions(
            "clones=s"  => \$clonelist,
            "ffa=s"     => \$ffa,
            "dbname=s"  => \$dbname,
            "dbhost=s"  => \$dbhost,
            "dbuser=s"  => \$dbuser,
            "help"      => \$help,
            "clobber"   => \$clobber,
            "info"      => \$info,
            "write"     => \$write
);

if ($help) {
    &usage;
    exit 0;
} elsif ($info) {
    exec("perldoc $0");
}

if (! $clonelist) {
    print STDERR "Must specify -clones\n";
    exit 1;
}

unless ($dbname && $dbuser && $dbhost) {
    print STDERR "Must specify all DB parameters\n";
    exit 1;
}

if ($clonelist eq 'all') {
    undef $clonelist;
}
else {
    open CLONES, "< $clonelist" or die "Can't open clone list $clonelist";
    while (<CLONES>) {
	chomp;
	$required_clones{$_} = 1;
    }
    close CLONES;
}


my $dbobj = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    '-host'   => $dbhost,
    '-user'   => $dbuser,
    '-dbname' => $dbname
) or die "Can't connect to DB";

my $seq = new FileHandle;;
if ($ffa =~ /\.gz$/) {
    open $seq, "gzcat $ffa |" or die "Can't open zipped file $ffa";
}
else {
    open $seq, "$ffa |" or die "Can't open plain file $ffa";
}

my $seqio = Bio::SeqIO->new(
    '-format' => 'fasta',
    '-fh'   => $seq
) or die "Can't read seq from file $ffa";


SEQ: while (my $seqobj = $seqio->next_seq) {

    my ($acc, $ver, $count) = _parse_id($seqobj->id);
    my ($offset, $end) = _parse_desc($seqobj->desc);
    my $length = $end - $offset + 1;

    # write the clone if this contig is from a different clone
    if (defined $clone && $acc ne $clone->id) {
	if ($write) {
	    $dbobj->write_Clone($clone);
	    $written_clone{$clone->id} = 1;
	    print STDERR "Written clone $acc.$ver\n";
	}
	undef $clone;
    }

    # skip unwanted clones
    unless (defined $clonelist && defined $required_clones{"$acc.$ver"}) {
	# print STDERR "Clone $acc.$ver not required - skip\n";
	next SEQ;
    }
    if (defined $inDB{"$acc.$ver"}) {
	print STDERR "Clone $acc.$ver already in DB - skip\n";
	next SEQ;
    }

    my $contigid = "$acc.$ver.$offset.$end";
    print STDERR "Found contig ", $seqobj->id, " ", $seqobj->desc, "\n";
    print STDERR "\tclone  $acc.$ver\n";
    print STDERR "\tid     $contigid\n";
    print STDERR "\toffset $offset\n";
    print STDERR "\tend    $end\n";
    print STDERR "\tlen    ", length($seqobj->seq), "\n";

    $contig = new Bio::EnsEMBL::PerlDB::Contig;
    $contig->embl_offset($offset);
    $contig->id($contigid);
    $contig->length($length);
    $contig->seq(new Bio::Seq('-id' => $acc, '-seq' => $seqobj->seq));
    $contig->version(1);
    $contig->embl_order($count);

    # this is for the case where a clone is not contiguous in the file
    # get clone -> add contig -> delete old clone -> write new one
    # not tested
    if ($written_clone{$acc}) {
	print STDERR "Fetching clone $acc for rewrite\n";
	my $tmpclone = $dbobj->get_Clone($acc);
	$clone = $tmpclone;
	$clone->add_Contig($contig);
	$tmpclone->delete;
	$dbobj->write_Clone($clone);
	next SEQ;
    }

    unless (defined $clone) {
	my $dbclone;
	# first check to see if we have already written that clone
	eval {
	    $dbclone = $dbobj->get_Clone_by_version($acc, $ver);
	};
	if (defined $dbclone) {
	    if ($clobber) {
		print STDERR "Overwriting $acc.$ver\n";
		$dbclone->delete;
	    }
	    else {
		print STDERR "Already have $acc.$ver - skip\n";
		$inDB{"$acc.$ver"} = 1;
		undef $contig;
		next SEQ;
	    }
	}
	$clone = new Bio::EnsEMBL::PerlDB::Clone;
	$clone->id($acc);
	$clone->embl_id($acc);
	$clone->version(1);
	$clone->embl_version($ver);
	print STDERR "Created clone $acc $ver\n";
    }

    $clone->add_Contig($contig);
}



sub _parse_desc {
# parse things like "(Z82169.1:1..110441)"
    my($desc) = @_;

    my($acc, $ver, $offset, $end) = $desc =~ /\((.+)\.(\d+):(\d+)\.\.(\d+)\)/;
    return $offset, $end;
}


sub _parse_id {
# parse things like "Z82169.1~1"
    my($id) = @_;

    my($acc, $ver, $count) = $id =~ /(.+)\.(\d+)~(\d+)/;
    return $acc, $ver, $count;
}


sub usage {
    print <<EOF
$0 [options]
Options:
  -clones   file containing list of clones
  -all      load all clones
  -dbname
  -dbhost
  -dbuser
  -ffa      location of ffa file
  -clobber  overwrite existing clone
EOF
}
