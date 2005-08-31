#!/usr/local/bin/perl

=head1 NAME

vega_repeat_libraries.pl - set repeat_consensus.repeat_class

=head1 SYNOPSIS

vega_repeat_libraries.pl [options]

General options:
    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --dbname, db_name=NAME              use database NAME
    --host, --dbhost, --db_host=HOST    use database host HOST
    --port, --dbport, --db_port=PORT    use database port PORT
    --user, --dbuser, --db_user=USER    use database username USER
    --pass, --dbpass, --db_pass=PASS    use database passwort PASS
    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)
    -v, --verbose                       verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

Specific options:

    --repeatfile=FILE                   read repeat class definitions from FILE

=head1 DESCRIPTION

This program classifies the repeats stored in a core database into some
somewhat sensible categories. It does this through a combination of a
repeat.txt file extracted from RepeatMasker repeat libraries and through some
simple pattern matching of the repeat names.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>
Patrick Meidl <pm2@sanger.ac.uk>

Based on code by James Smith <js5@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options('repeatfile');
$support->allowed_params($support->get_common_params, 'repeatfile');

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params('repeatfile');

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

# mouse fixes
if ($support->species eq 'Mus_musculus') {
    $support->log("Making Vega mouse specific changes...\n");
    $support->log("Copying repeat_name to repeat_consensus...\n", 1);
    $dbh->do("update repeat_consensus set repeat_consensus = repeat_name where repeat_class = 'Tandem_repeat'") unless ($support->param('dry_run'));
    $support->log("Setting repeat_name to 'trf' where appropriate\n", 1);
    $dbh->do("update repeat_consensus set repeat_name = 'trf' where repeat_class = 'Tandem_repeat'") unless ($support->param('dry_run'));
    $support->log("Done.\n");
}

# clear repeat_class
$support->log("Clearing repeat_class...\n");
$dbh->do("update repeat_consensus set repeat_class = ''") unless ($support->param('dry_run'));
$support->log("Done.\n");

# read repeat classes from file
$support->log_stamped("Reading repeat classes from input file...\n");
my $fh = $support->filehandle('<', $support->param('repeatfile'));
my $C = 0;
while (<$fh>) {
    chomp;
    my ($hid, $type) = split( /\t/, $_, 2);
    $dbh->do("update repeat_consensus set repeat_class = ? where repeat_name in (?,?,?)", {} , $type, $hid, substr($hid,0,15), "$hid-int" ) unless ($support->param('dry_run'));
    $C++;
    $support->log("$C\n", 1) unless $C % 100;
}
close $fh;
$support->log_stamped("Done.\n");

# Consensifying repeat classes
$support->log_stamped("Consensifying remaining repeat classes...\n");
unless ($support->param('dry_run')) {
    $dbh->do("update repeat_consensus set repeat_class = 'Simple_repeat'  where repeat_class= '' and repeat_name like '%)n'" );
    $dbh->do("update repeat_consensus set repeat_class = 'low_complexity'  where repeat_class= '' and repeat_name like '%-rich'" );
    $dbh->do("update repeat_consensus set repeat_class = 'low_complexity'  where repeat_class= '' and repeat_name like 'poly%'" );
    $dbh->do("update repeat_consensus set repeat_class = 'LTR/ERVL'  where repeat_class= '' and repeat_name like '%ERVL%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'LTR/ERVL'  where repeat_class= '' and repeat_name like '%ERV16%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'SINE/Alu'  where repeat_class= '' and repeat_name like 'Alu%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'SINE/Alu'  where repeat_class= '' and repeat_name like '%F_AM%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'LINE/L1'  where repeat_class= '' and repeat_name like 'L1%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER2_type'  where repeat_class= '' and repeat_name like 'Tigger%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'Charlie%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/Tc2'  where repeat_class= '' and repeat_name like 'HsTC%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER2_type'  where repeat_class= '' and repeat_name like 'MER46%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER2_type'  where repeat_class= '' and repeat_name like 'MER7%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'MER91' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'MER58' " );
    $dbh->do("update repeat_consensus set repeat_class = 'DNA/MER1_type'  where repeat_class= '' and repeat_name like 'MER63' " );
    $dbh->do("update repeat_consensus set repeat_class = 'Satellite/telomeric'  where repeat_class= '' and repeat_name like 'SUBTEL_%' " );
    $dbh->do("update repeat_consensus set repeat_class = 'trf'  where repeat_class = '' and repeat_name = 'trf' " );
    $dbh->do("update repeat_consensus set repeat_class = 'dust' where repeat_class = '' and repeat_name = 'dust'" );
    $dbh->do("update repeat_consensus set repeat_class = 'novel_transposon' where repeat_class = '' and repeat_name = 'novel_transposon'");
}
$support->log_stamped("Done.\n");

# Setting repeat types
$support->log_stamped("Setting repeat types...\n");
my %mappings = (
        'Low_Comp%' => 'Low complexity regions',
        'LINE%'	=> 'Type I Transposons/LINE',
        'SINE%'	=> 'Type I Transposons/SINE',
        'DNA%'	=> 'Type II Transposons',
        'LTR%'	=> 'LTRs',
        'Other%'	=> 'Other repeats',
        'Satelli%'	=> 'Satellite repeats',
        'Simple%'	=> 'Simple repeats',
        'Other%'	=> 'Other repeats',
        'Tandem%'	=> 'Tandem repeats',
        'TRF%'	=> 'Tandem repeats',
        'dust%' => 'Dust',
        'Unknown%'	=> 'Unknown',
        '%RNA'	=> 'RNA repeats',
        'novel_transposon' => 'Novel Transposon',
);
unless ($support->param('dry_run')) {
    foreach (keys %mappings) { 
        $dbh->do(qq(update repeat_consensus set repeat_type = '$mappings{$_}' where repeat_class like '$_'));
    }

    # type all remaining repeats as unknown
    $dbh->do(qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = ''));
    $dbh->do(qq(update repeat_consensus set repeat_type = 'Unknown' where repeat_type = NULL));
}
$support->log_stamped("Done.\n");

# finish logfile
$support->finish_log;

