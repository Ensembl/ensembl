#!/bin/ksh

# $Id$
#
# Author:  Andreas Kahari, <andreas.kahari@ebi.ac.uk>
#
# This is a wrapper around exonerate that produces a comma
# separated file of alignment information, including cigar lines
# and separate percentage of identity for query and target
# sequences.  Run with -h or -? for usage information.
#
#----------------------------------------------------------------------
#
# Format of output:
#
# Col	Value				Type
# 1	Query sequence ID		double quoted string
# 2	Query percentage of identity	%g-type float <= 100
# 3	Alignment start in query	integer
# 4	Alignment end in query		integer
# 5	Target sequence ID		double quoted string
# 6	Target percentage of identity	%g-type float <= 100
# 7	Alignment start in target	integer
# 8	Alignment end in target		integer
# 9	Score				integer
# 10	Cigar line			double quoted string
#
#----------------------------------------------------------------------
#
# Variables:
#
# q_fa:   Query FastA file (-q flag).
# t_fa:   Target FastA file (-t flag).
# e_cmd:  Which exonerate executable to pick up (-e flag).
# e_ryo:  Exonerate output format specification.
# [qts]_min: Minimum values on query, target, and score percentages
#            (min score is as a percentage of best score per query)
#            (-Q, -T, and -S flags).
#
#----------------------------------------------------------------------

q_fa='/acari/work4/mongin/anopheles_mai/qc/peptides.fa'
t_fa='/acari/work4/mongin/anopheles_mai/qc/submitted_genes.fa'
e_cmd='/acari/work2/gs2/gs2/local/OSF1/bin/exonerate'
q_min=25
t_min=25
s_min=95

function usage
{
	typeset -R${#0} padding=' '

	cat >&2 <<-EOT
	Usage:	$0 [-h?]
	        $0 [-e path] [-q path] [-t path] [-v]
	        $padding [-Q val] [-T val] [-S val] [-- opts]	

	-h, -?   Show usage information (this text).

	-e path  Explicit path to exonerate executable.

	-q path  Path to FastA file containing query sequences.

	-t path  Path to FastA file containing target sequences.

	-f path  Don't run exonerate, read from this properly formatted
	         file instead.

	-v       Be verbose.

	-Q val   Don't consider alignments which have less than this
	         query percentage of identity.

	-T val   Don't consider alignments which have less than this
	         target percentage of identity.

	-S val   Don't consider alignments which have a score less than
	         this many percent of the best score for each query
	         sequence.

	opts     Extra options to pass to exonerate, e.g. -M 512 -s 200
	         (must not affect format of output!).

	Default values:
	  -e $e_cmd
	  -q $q_fa
	  -t $t_fa
	  -Q $q_min
	  -T $t_min
	  -S $s_min
	EOT
}

while getopts 'h?e:q:t:f:vQ:T:S:' opt; do
    case $opt in
	e) e_cmd=$OPTARG ;;
	q) q_fa=$OPTARG  ;;
	t) t_fa=$OPTARG  ;;
	f) e_out=$OPTARG ;;
	Q) q_min=$OPTARG ;;
	T) t_min=$OPTARG ;;
	S) s_min=$OPTARG ;;
	v) set -o xtrace ;;
	*) usage; exit 1 ;;
    esac
done
shift $(( OPTIND - 1 ))

if [[ ! -r $q_fa ]]; then
    print -u2 "Can't find or read file '$q_fa'"
    exit 1
elif [[ ! -r $t_fa ]]; then
    print -u2 "Can't find or read file '$t_fa'"
    exit 1
elif [[ ! -x $e_cmd && ! -x $(whence $e_cmd) ]]; then
    print -u2 "Can't find or execute command '$e_cmd'"
    exit 1
fi

# Set default options for exonerate using environment variables (may be
# overridden by command line options).
export EXONERATE_EXONERATE_FSMMEMORY=512
export EXONERATE_EXONERATE_HSPDROPOFF=5
export EXONERATE_EXONERATE_HSPTHRESHOLD=30
export EXONERATE_EXONERATE_MODEL='affine:local'
export EXONERATE_EXONERATE_PROTEINWORDLEN=6
export EXONERATE_EXONERATE_WORDTHRESHOLD=3

e_ryo='%qi\t%qal\t%ql\t%qab\t%qae\t%ti\t%tal\t%tl\t%tab\t%tae\t%p\t%s\t%C\n'
e_opt="--showalignment no --showvulgar no -q $q_fa -t $t_fa --ryo $e_ryo $@"

if [[ -z $e_out ]]; then
    cmd="$e_cmd $e_opt"
else
    cmd="cat $e_out"
fi

$cmd |
perl -ne '
    # Perl script to calculate the percentage of identity
    # and reformat cigar lines.  Takes tab-delimited list in
    # specific format from exonerate as input on stdin and
    # writes comma-separated output to stdout.

    BEGIN {
	# Pick out values from the shell script.
	$q_min = '$q_min';
	$t_min = '$t_min';
	$s_min = '$s_min' / 100;
    }

    next if (/^Message:/ || /^--/);
    chomp;

    # Pick out the individual fields (variable names correspond to
    # exonerate format parameters).
    ($qi, $qal, $ql, $qab, $qae,
     $ti, $tal, $tl, $tab, $tae,
     $p,  $s,   $C) = split /\t/;

    # Calculate the percent of identity (and skip to the next alignment
    # if it is not good enough).
    $qp = $qal * $p / $ql; next if ($qp < $q_min);
    $tp = $tal * $p / $tl; next if ($tp < $t_min);

    # Only store the incoming result if its score is larger than 95% (or
    # whatever is in $s_min) of the best score found so far for this
    # query ID.  If the score is larger than the best score, or if there
    # is no best score yet, then update the best score.

    if (!defined $r{$qi}{BS}) {
	# This is the first score.
	$r{$qi}{BS} = $s;
    }
    if ($s >= $s_min * $r{$qi}{BS}) {
	# This is a good enough score, as far as we know now.

	# Reformat the cigar line.
	$C =~ s/([MDI]+) ([0-9]+) ?/$2$1/g;  # flip
	$C =~ s/([MDI])1([MDI])/$1$2/g;      # no lone 1
	$C =~ s/^1([MDI])/$1/;               # no lone 1 at start of line

	push @{ $r{$qi}{ALGN} }, {
	    qi => $qi, qp => $qp, qab => $qab, qae => $qae,
	    ti => $ti, tp => $tp, tab => $tab, tae => $tae,
	    s  => $s,  C  => $C };

	if ($s > $r{$qi}{BS}) {
	    # This score was better than the best score so far.
	    $r{$qi}{BS} = $s;
	}
    }

    END {
	# Stuff to do at the end (output).

	foreach $qi (sort keys %r) {
	    foreach $t (@{ $r{$qi}{ALGN} }) {
		next if ($t->{s} < $s_min * $r{$qi}{BS});

		printf "\"%s\",%g,%d,%d,\"%s\",%g,%d,%d,%d,\"%s\"\n",
		    $t->{qi}, $t->{qp}, $t->{qab}, $t->{qae},
		    $t->{ti}, $t->{tp}, $t->{tab}, $t->{tae},
		    $t->{s}, $t->{C};
	    }
	}
    } # END
'
