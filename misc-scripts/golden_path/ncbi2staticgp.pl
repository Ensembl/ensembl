use strict;

# ncbi assembl -> ensembl static golden path
# author: simon potter, 09/01
# based loosely on the ucsc assembly parser
# (c) ebi/grl 2001

# takes 3 file names from cmd line:
# 1. loc of nt_contigs on chromosome (seq_contig.md)
# 2. agp - clones in nt_contigs (allcontig.agp.build25)
# 3. ensembl raw_contigs (contig.txt).
# [examples at end of file]

# "stategy":
# read nt_contig/chromosome file and store start/end/orientation
# of nt_contigs on the chromosome;
# parse list of ensembl raw contig coords and store in uberhash;
# go through agp linewise and look for matching raw contigs


my($chrom, $agp, $raw) = @ARGV or die "horribly";
my @chr_len;
my %nt_contig;
my %clone;

my $default_ori = '+';
my $chr_offset = 1;    # chr coordinates appear to be off by 1


open CHR, "< $chrom" or die "Can't open chromosome file $chrom";
open AGP, "< $agp"   or die "Can't open agp file $agp";
open CTG, "< $raw"   or die "Can't open raw contig file $raw";

open SGP, "> static.txt" or die "Can't open static.txt for output";
open INF, "> idlist.txt" or die "Can't open idlist.txt for output";


while (<CHR>) {
    my ($chr, $start, $end, $ori, $nt_ctg, $obj) = (split)[1, 2, 3, 4, 5, 7];
    if ($obj eq 'start') {
	$chr_len[$chr]->[0] = $start;
	next;
    }
    if ($obj eq 'end') {
	$chr_len[$chr]->[1] = $end;
	next;
    }
    $ori = $default_ori unless $ori =~ /^[+-]$/;
    $nt_contig{$nt_ctg}->{'chr'}   = $chr;
    $nt_contig{$nt_ctg}->{'start'} = $start + $chr_offset;  # nt_ctg loc on chromosome
    $nt_contig{$nt_ctg}->{'end'}   = $end + $chr_offset;
    $nt_contig{$nt_ctg}->{'ori'}   = $ori;
}
close CHR;


while (<CTG>) {
    my ($contig_id, $internal_id, $start, $length) = split;
    $contig_id =~ /([^.]+\.[^.]+)\./ || die "Bad id $contig_id";
    my $clone_id = $1;
    if (!defined $clone{$clone_id}) {
        $clone{$clone_id} = [];
    }

    push @{$clone{$clone_id}}, {
	"id"    => $contig_id,
	"iid"   => $internal_id,
	"start" => $start,
	"end"   => $start + $length - 1
    };
}
close CTG;


while (<AGP>) {
    my ($nt_ctg, $nt_start, $nt_end, $sv, $raw_start, $raw_end, $raw_ori) =
     (split)[0, 1, 2, 5, 6, 7, 8];
    next if $raw_start eq 'fragment';
    next unless $sv =~ m{^\S+\.\d+$};
    unless (($raw_end - $raw_start) == ($nt_end - $nt_start)) {
	die "Raw contig and nt contig coords don't match";
    }

    my $nt_ori = $nt_contig{$nt_ctg}->{'ori'};
    my $chr    = $nt_contig{$nt_ctg}->{'chr'};
    my ($chr_start, $chr_end);

    if ($nt_ori eq '+') {
	# forward oriented nt contig: raw contigs forward from nt contig start
	$chr_start = $nt_contig{$nt_ctg}->{'start'} + $nt_start - 1;
	$chr_end   = $nt_contig{$nt_ctg}->{'start'} + $nt_end   - 1;
    }
    else {
	# reverse oriented nt contig: raw contigs back from nt contig end
	$raw_ori *= -1;  # flip raw contig if nt contig is in reverse ori
	$chr_start = $nt_contig{$nt_ctg}->{'end'} - $nt_end   + 1;
	$chr_end   = $nt_contig{$nt_ctg}->{'end'} - $nt_start + 1;
    }

    my $seen = 0;
    unless (exists $clone{$sv}) {
	print STDERR "No clone found: $sv\n";
	next;
    }
    foreach my $raw_ctg (@{$clone{$sv}}) {
	my $start = $raw_ctg->{'start'};
	my $end   = $raw_ctg->{'end'};
	if ($raw_start >= $start && $raw_start <= $end) {
	    if ($raw_end > $end) {
		print STDERR $raw_ctg->{'id'}, " doesn't fit into ",
		 "$raw_start:$raw_end\n";
	    }
	    else {
		$seen = 1;
		$raw_start = $raw_start - $start + 1;
		$raw_end   = $raw_end   - $start + 1;

		print SGP "$nt_ctg\t$chr\t", $raw_ctg->{'iid'};
		print SGP "\t$chr_start\t$chr_end\t$nt_start\t$nt_end";
		print SGP "\traw_start\t$raw_end\t$raw_ori\tNCBI\n";

		print INF $raw_ctg->{'id'}, "\t", $raw_ctg->{'iid'}, "\n";

		# print "$chr $nt_ctg:$nt_start:$nt_end ";
		# print $raw_ctg->{'id'}, ":$raw_start:$raw_end:$raw_ori ";
		# print "$chr_start:$chr_end\n";

		last;
	    }
	}
    }
    if ($seen == 0) {
	print STDERR "Can't fit clone $sv: $raw_start - $raw_end\n";
    }
}
close AGP;


__END__

seq_contig.md:

9606    1       0       0       +       start   -1      contig  10
9606    1       0       644775  ?       NT_021903.5     21903   contig  1
9606    1       694776  1755501 -       NT_004350.5     4350    contig  1
9606    1       1805502 2002021 -       NT_028057.1     28057   contig  1

allcontig.agp.build25:

NT_004615.5          1       4637   1  D  AL390868.5     144054   148690  -
NT_004615.5       4638       4737   2  N       100 fragment no
NT_004615.5       4738      15035   3  D  AC021696.3      87558    97855  +
NT_004615.5      15036      15135   4  N       100 fragment no
NT_004615.5      15136      20250   5  D  AC021696.3      19497    24611  -

contig.txt (id, internal_id, offset, length):

AB000381.1.1.35863      1       1       35863
AB012723.1.1.40850      11      1       40850
AB015355.1.1.43999      24      1       43999
AB015752.1.1.116160     25      1       116160
AB016897.1.1.331211     26      1       331211
