#!/usr/local/bin/perl
# $Id$ 

# Script to translate fgenesh++ gtf files from chromosome+chr_coords to
# fpcctg_name + fpc_coords.  Is a bit tricky due to our golden path
# missing (parts of) fpc_contigs.
# This is used for the IGI work.

# cared for by Philip lijnzaad@ebi.ac.uk

use DBI;
use strict;

use vars qw($opt_h $opt_C);

use Getopt::Long;
# use Bio::EnsEMBL::Utils::GTF_Merge('gtf_merge');
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;

my $Usage = "Usage: $0 [ -h (help) ]  < chromo-coords.gtf > fpc-coords.gtf [ >& log ]\n";
my $opts = 'hC:';

my $connection = ($opt_C || 'host=ecs1c;user=ensro;database=ensembl080;');

#Database options
my $host   = 'ecs1c';
my $port   = undef;
my $dbname = 'ensembl080';
my $dbuser = 'ensro';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $help;

&GetOptions( 
	     'host:s'     => \$host,
	     'dbname:s'   => \$dbname, 
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'h|help'     => \$help,
	     );
die $Usage if $help;

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";

my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my %bighash = undef;

&read_mapping($db);

# $db->static_golden_path_type('UCSC');
# my $stadaptor = $db->get_StaticGoldenPathAdaptor();

while(<>)  {
    chomp($_);
    my ($chr, $source, $tag, $start, $end, $score, $strand, $frame, @rest) 
      = split("\t");
### depracated, much too slow (2 s per query)
###    ($chr, $start, $end) = $stadaptor->convert_chromosome_to_fpc($chr,$start, $end);

    #reformat a bit as well:
    grep( s/ ([^ ]+)/ "$1";/ , @rest);
    my $rest = join(' ', @rest);
    chomp $rest;

    my ($contig, $ctg_start, $ctg_end, $status) = map_coords($chr,$start, $end);

    if (! $status ) {
        print join("\t", ($contig, $source, $tag, 
                          $ctg_start, $ctg_end, $score, 
                          $strand, $frame, $rest)), "\n";
    } elsif ($status eq 'init' ) {
        warn "$chr: ignoring it, but it might lie on missing initial part of $contig: $_\n";
    } elsif ($status eq 'missing') { 
        warn "$chr: cannot find contig for $_; ignoring it\n";
    } elsif ($status eq 'bug') {
        warn "$chr: bug with $_; ignoring it\n";
    } else {
        die "Bugger";               # even bigger bug
    }
}                                       # while

sub read_mapping {
    my ($db) = shift;
    my ($q) = "SELECT chr_name, MIN(chr_start), MAX(chr_end), fpcctg_name, MIN(fpcctg_start), MAX(fpcctg_end)
               FROM static_golden_path
               GROUP BY fpcctg_name
               ORDER BY chr_name, chr_start";
    $q=$db->prepare($q) || die $db->errstr;
    $q->execute() || die $db->errstr;

    while ( my ($chr, $chr_start, $chr_end, $fpcctg_name, $fpcctg_start, $fpcctg_end )
                 = $q->fetchrow_array ) {
        die $db->errstr if $q->err;
        push @{$bighash{$chr}}, [ $chr_start, $chr_end, $fpcctg_name, $fpcctg_start ];
        # structure : hash{chr} contains sorted list of these arrays 

    }
    return;
}                                       # read_mapping

# return fpc_ctg_name, start and end, status, given the chromosome+start/end
sub map_coords {
    my ($chr, $start, $end) = @_;

    my $contigs = $bighash{$chr};
    my $n = 0;

    if (! defined $contigs) {
        warn "don't know $chr";
        return ('ctg_on_unknown_chromosome', -999, -999, 'bug');
    }

    foreach my $l ( @{$contigs} ) {
        if ( $start < $l->[1] && $start >= $l->[0] ) { 
            # found it
            my ( $chr_start, $chr_end, $fpcctg_name, $fpcctg_start)
              = @$l;
            return ($fpcctg_name,
                    $fpcctg_start-1+ $start - $chr_start,
                    $fpcctg_start-1+ $end - $chr_start, undef);
        }
    }
    # this happens (only?) if we have an fpccontig with start != 1; 
    # look for the missing bit. 
    my $prev = undef;
    foreach my $curr ( @{$contigs} ) {
        if (defined $prev 
            && $start > $prev->[1] && $start < $curr->[0] ) { 
            ## see if we were correct
            my ( $chr_start, $chr_end, $fpcctg_name, $fpcctg_start)
              = @$curr; 
            if ( $fpcctg_start == 1 )  {
                # this means we can't pretend that it's on the missing
                # piece of this contig; so this contig is simply missing
                # from our static_golden_path
                return ('missing_ctg', -999, -999, 'missing');
            }
            #else:
            # fpccontig_start is not 1; now just pretend this feature
            # is on the 'unknown' first bit of the current fpc contig,
            # and pretend our fpccontig actually did start at 1:
            my $status = 'init';
            $chr_start -= ($fpcctg_start -1);
            $fpcctg_start =1;
            my $newstart  = $fpcctg_start-1+ $start - $chr_start;
            $status = 'missing' if ($newstart < 1); 
            # (i.e., it is too far away from our fpcctg_start to still be
            # able to pretend it's still on this contig, so report as missing)

            return ($fpcctg_name,
                    $newstart, 
                    $fpcctg_start-1+ $end - $chr_start, $status);
        }
        $prev = $curr;
    }
    warn "really couldn find coords of $chr $start $end\n";
    return ('bug_ctg', -999, -999, 'bug');
}                                       # map_coords

