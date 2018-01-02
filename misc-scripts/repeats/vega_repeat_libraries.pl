#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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

    --prune				undo, i.e. delete from the database changes caused by running the script			


Specific options:

    --repeatfile=FILE                   read repeat class definitions from FILE

=head1 DESCRIPTION

This program classifies the repeats stored in a core database into some
somewhat sensible categories. It does this through a combination of a
repeat.txt file extracted from RepeatMasker repeat libraries and through some
simple pattern matching of the repeat names.


=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>
Patrick Meidl <pm2@sanger.ac.uk>

Based on code by James Smith <js5@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list http://lists.ensembl.org/mailman/listinfo/dev

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
$support->parse_extra_options('repeatfile=s', 'prune');
$support->allowed_params($support->get_common_params, 'repeatfile', 'prune');

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params('repeatfile') unless $support->param('prune');		# don't need the repeat file for pruning

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

# unless we are pruning (undo), we should make a backup copy of the repeat_consensus table
if($support->param('prune')){
	# prune (undo mode)
	# backup table must exist for this to work
	
	if(check_for_backup_table()){
		# backup table present
		if($support->user_proceed("Replace the current table 'repeat_consensus' with the backup table 'repeat_consensus_backup'?")){		
			if($dbh->do("drop table repeat_consensus")){				
				if($dbh->do("create table repeat_consensus select * from repeat_consensus_backup")){					
					$support->log("prune (undo) was successful\n");				
					$support->log_stamped("Done.\n");

					# finish logfile
					$support->finish_log;				
					exit(0);								
				}
				else {				
					$support->log_error("prune failed\n");								
				}				
			}
			else {			
				$support->log_error("prune failed\n");			
			}				
		}
		else{
		
			#user is aborting
			print "aborting...\n";
			$support->log_error("aborting...\n");				
		}	
	}
	else{	
		print "Cannot do prune, as no backup table\n";
		$support->log_error("Cannot do prune, as no backup table\n");		
	}
}
else{

	# normal run
	# check to see if the backup table 'repeat_consensus_backup' already exists		
	if(check_for_backup_table()){
		#table already exists: ask user if OK to overwrite it		
		if ($support->user_proceed("The backup table 'repeat_consensus_backup' already exists, OK to delete?")) {
                    if($dbh->do("drop table 'repeat_consensus_backup'")){
                    	$support->log("deleted previous backup table\n");
                    	make_backup_table();
                    }
                    else{
                    	$support->log_error("tried but failed to delete previous backup table\n");
                    }
                }
                else{
                	# user won't allow removing the backup table
                	print "Aborting ...\n";
                	$support->log_error("User won't allow removal of backup table ... aborting program\n");
                }
	}else{
		# table doesn't exist, therefore we can create it
		make_backup_table();
	}
}


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


sub make_backup_table{
	if($dbh->do("create table repeat_consensus_backup select * from repeat_consensus")){
		$support->log("backup table 'repeat_consensus_backup was created successfully\n");	
	}
	else{
		$support->log_error("failed to create backup table 'repeat_consensus_backup'\n");
	}
}

sub check_for_backup_table{
	# check to see if the backup table 'repeat_consensus_backup' already exists
	my @tables = $dbh->tables();
	my $found=0;

	foreach my $table(@tables){
		#print "$table\n";
		
		if($table eq '`repeat_consensus_backup`'){			
			$found=1;
			last;		
		}	
	}
	return $found;
}
