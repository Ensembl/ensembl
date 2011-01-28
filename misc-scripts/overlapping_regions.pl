#!/usr/local/ensembl/bin/perl -w

=head1 NAME

overlapping_regions.pl - calculates overlapping regions in the assembly table

=head1 SYNOPSIS

overlapping_regions.pl [arguments]

Required arguments:

  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  
  --pattern, --dbpattern=PATTERN      patch databases where name matches PATTERN
                                      Note that this is a database pattern of
                                      the form %core_41% rather than a regexp
 
Optional arguments:

  -n, --dry, --dry_run=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

Please note that where an argument expects a value, this is true for all
alternative argument styles.

=head1 DESCRIPTION

This is a script that will calculate if an assembly has overlapping regions
in same coordinate system (e.g. 2 chr_1 stored in the assembly table, one as
the current, the other as an old version). It will store a new entry in the meta
table. If there are overlapping regions, can cause problems when moving features
between coordinate systems. The pattern is a string that can be used in an
'IN' clause in SQL (e.g. "%_core_%"; if you want to patch only a single
database, use a pattern without % expansion, e.g. "homo_sapiens_core_38_36").

=head1 EXAMPLES

   perl overlapping_regions.pl -host ens-staging -user ensadmin -pass xxxx
        -pattern "%core_55_%" -port 3306 -dry_run 0

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Daniel Rios <dani@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<dev@ensembl.org>

=cut


use strict;
use warnings;
no warnings 'uninitialized';
use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use DBI;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
  'pattern|dbpattern=s'
			      );
my @params = map { $_ unless ($_ =~ /dbname/) } $support->get_common_params;
$support->allowed_params(
  @params,
  'pattern'
			 );

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'pattern'
);

# connect to database
my $dbh = $support->get_dbconnection;

# get all database names that match pattern
my ($sth, $sql);
$sql = "SHOW DATABASES LIKE '".$support->param('pattern')."'";
$sth = $dbh->prepare($sql);
$sth->execute;

# loop over databases
while (my ($dbname) = $sth->fetchrow_array) {
   $support->log_stamped("$dbname\n");
  
  if ($support->user_proceed("\nCalculate overlapping regions for $dbname?")) {
      #this query will detect if there are overlapping regions in same coord_system
      #starting in position 1. 
      $dbh->do("use $dbname;");
      $sql = qq(SELECT count(*) FROM assembly as1, seq_region s1, assembly as2, seq_region s2 
		WHERE as1.asm_seq_region_id = s1.seq_region_id AND as2.asm_seq_region_id = s2.seq_region_id 
		AND s1.coord_system_id = s2.coord_system_id AND as1.asm_start = as2.asm_start 
		AND as1.asm_end = as2.asm_end AND as1.cmp_seq_region_id = as2.cmp_seq_region_id 
		AND as1.asm_seq_region_id <> as2.asm_seq_region_id AND as1.asm_start = 1);
      my $sth1 = $dbh->prepare($sql);
      $sth1->execute;
      my $val = $sth1->fetch();
      my $overlaps;
      #let's see if we write it in the database
      if ($val->[0] > 0){
	  $support->log("\tDatabase has overlapping regions. This might make give wrong results when mapping features\n");
	  $overlaps = 'true';
      }
      else{
	  $support->log("\tDatabase looks correct\n");
	  $overlaps = 'false';
      }
      if (defined $support->param('dry_run') and !$support->param('dry_run')){
	  #let's write results to the database. First remove previous entry
	  $sql = qq(DELETE FROM meta where meta_key = 'assembly.overlapping_regions');
	  $dbh->do($sql);
	  #and now will insert it
	  $sql = qq(INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,'assembly.overlapping_regions',$overlaps));
	  $dbh->do($sql);
      }
     $support->log("Done\n"); 
  }

}
$support->finish_log;
