#!/bin/ksh

# $Id$
#
# Author:  Andreas Kahari, <andreas.kahari@ebi.ac.uk>
#
# This is a wrapper around exonerate that produces a comma separated file
# of alignment information, including cigar lines and separate percentage
# of identity for query and target sequences.  Run with -h or -? for
# usage information.
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
# e_mod:  What alignment mode exonerate should use.
# e_fmt:  The format for the output from exonerate.
# e_opt:  The aggregated options to pass to exonerate.
#
#----------------------------------------------------------------------

q_fa='/acari/work4/mongin/anopheles_mai/qc/peptides.fa'
t_fa='/acari/work4/mongin/anopheles_mai/qc/submitted_genes.fa'
e_cmd='/acari/work2/gs2/gs2/local/OSF1/bin/exonerate'

function usage
{
	cat >&2 <<-EOT
	Usage:	$0 [-h?]
	        $0 [-q path] [-t path] [-v]

	-h, -?   Show usage information.
	-e path  Explicit path to exonerate executable.
	-q path  Path to FastA file containing query sequences.
	-t path  Path to FastA file containing target sequences.
	-v       Be verbose.

	Default values:
	  -e $e_cmd
	  -q $q_fa
	  -t $t_fa
	EOT
}

e_mod='affine:local'

e_fmt='%qi\t%qal\t%ql\t%qab\t%qae\t%ti\t%tal\t%tl\t%tab\t%tae\t%p\t%s\t%C\n'

while getopts 'h?e:q:t:v' opt; do
    case $opt in
	e) e_cmd=$OPTARG ;;
	q) q_fa=$OPTARG  ;;
	t) t_fa=$OPTARG  ;;
	v) set -o xtrace ;;
	*) usage; exit 1 ;;
    esac
done

if [[ ! -r $q_fa ]]; then
    print -u2 "Can't find file '$q_fa'"
    exit 1
elif [[ ! -r $t_fa ]]; then
    print -u2 "Can't find file '$t_fa'"
    exit 1
elif [[ ! -x $e_cmd && ! -x $(whence $e_cmd) ]]; then
    print -u2 "Can't execute '$e_cmd'"
    exit 1
fi

e_opt="--showalignment no --showsugar no --showcigar no --showvulgar no 
       --model $e_mod --ryo $e_fmt -q $q_fa -t $t_fa"

nice -n 19 $e_cmd $e_opt |
perl -ne '
    # Perl script to calculate the percentage of identity.
    # Takes tab-delimited list in specific format from exonerate as
    # input on stdin and writes comma-separated output to stdout.

    next if (/^Message:/ || /^--/);
    chomp;

    # Pick out the individual fields (variable names correspond to
    # exonerate format parameters).
    ($qi, $qal, $ql, $qab, $qae,
     $ti, $tal, $tl, $tab, $tae,
     $p,  $s,   $C) = split /\t/;

    # Calculate the percent of identity.
    $qp = 100 * $qal*($p/100)/$ql;
    $tp = 100 * $tal*($p/100)/$tl;

    # Reformat the cigar line.
    $C =~ s/([MDI]+) ([0-9]+) ?/$2$1/g;	# flip
    $C =~ s/([MDI])1([MDI])/$1$2/g;	# no lone 1
    $C =~ s/^1([MDI])/$1/;		# no lone 1 at start of line

    printf "\"%s\",%g,%d,%d,\"%s\",%g,%d,%d,%d,\"%s\"\n",
	$qi, $qp, $qab, $qae,
	$ti, $tp, $tab, $tae, $s, $C;
'
