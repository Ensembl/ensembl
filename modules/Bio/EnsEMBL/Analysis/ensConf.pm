#
# BioPerl module for Bio::EnsEMBL::Analysis::ensConf;
#
# Cared for by Tim Hubbard <th@sanger.ac.uk>
#
# Copyright Tim Hubbard, James Gilbert
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::ensConf - imports global variables used by EnsEMBL sequence analysis

=head1 SYNOPSIS

    use ensConf;
    use ensConf qw( HUMACESERVER_HOST HUMACESERVER_PORT );

=head1 DESCRIPTION

ensConf is a pure ripoff of humConf written by James Gilbert.

humConf is based upon ideas from the standard perl Env environment
module.

It imports and sets a number of standard global variables into the
calling package, which are used in many scripts in the human sequence
analysis system.  The variables are first decalared using "use vars",
so that it can be used when "use strict" is in use in the calling
script.  Without arguments all the standard variables are set, and
with a list, only those variables whose names are provided are set.
The module will die if a variable which doesn't appear in its
C<%ensConf> hash is asked to be set.

The variables can also be references to arrays or hashes.

Edit C<%ensConf> to add or alter variables.

All the variables are in capitals, so that they resemble environment
variables.

=head1 CONTACT

B<Tim Hubbard> email th@sanger.ac.uk
B<James Gilbert> email jgrg@sanger.ac.uk

=cut

#'

package Bio::EnsEMBL::Analysis::ensConf;

use strict;
use vars qw( %ensConf );

# Could change user in future
#my $humpub = ( getpwnam('humpub') )[7];
#my $ftp    = ( getpwnam('ftp')    )[7];

my $sanger_path="/nfs/disk100/humpub1a/unfinished_ana";
my $ext_path="/nfs/disk100/humpub2a/unfinished_ana";
my %cgp_path=map {$_,"$sanger_path/$_"} qw ( SU SF );
$cgp_path{'EU'}="$ext_path/EU";
$cgp_path{'EF'}="$ext_path/EF";

# Hash containing config info
%ensConf = (
	    UNFIN_ROOT => $sanger_path,
	    UNFIN_DATA_ROOT => $sanger_path,
	    UNFIN_DATA_ROOT2 => $ext_path,
	    UNFIN_DATA_ROOT_CGP => \%cgp_path,
            CONFIRMED_EXON_FASTA => "$sanger_path/confirmed_exon",
	    EXON_ID_SUBSCRIPT => 'ENSE',
	    EXON_ID_DIGITS => 11,
	    TRANSCRIPT_ID_SUBSCRIPT => 'ENST',
	    TRANSCRIPT_ID_DIGITS => 11,
	    GENE_ID_SUBSCRIPT => 'ENSG',
	    GENE_ID_DIGITS => 11,
	    PROTEIN_ID_SUBSCRIPT => 'ENSP',
	    PROTEIN_ID_DIGITS => 11,
	    CONTIG_CLUSTER_ID_SUBSCRIPT => 'HCC',
            WWW4EBI => "/nfs/WWW/htdocs/Users/humpub/ensembl",
	    );

sub import {
    my ($callpack) = caller(0); # Name of the calling package
    my $pack = shift; # Need to move package off @_

    # Get list of variables supplied, or else
    # all of ensConf:
    my @vars = @_ ? @_ : keys( %ensConf );
    return unless @vars;

    # Predeclare global variables in calling package
    eval "package $callpack; use vars qw("
         . join(' ', map { '$'.$_ } @vars) . ")";
    die $@ if $@;


    foreach (@vars) {
	if ( defined $ensConf{ $_ } ) {
            no strict 'refs';
	    # Exporter does a similar job to the following
	    # statement, but for function names, not
	    # scalar variables:
	    *{"${callpack}::$_"} = \$ensConf{ $_ };
	} else {
	    die "Error: ensConf: $_ not known\n";
	}
    }
}

1;
