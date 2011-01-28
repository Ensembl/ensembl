#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;

=head1 NAME

schema_patch.pl - automagically apply schema patches to Ensembl dbs

=head1 SYNOPSIS

schema_patch.pl [arguments]

Required arguments:

  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  
  --pattern, --dbpattern=PATTERN      patch databases where name matches PATTERN
                                      Note that this is a database pattern of
                                      the form %core_41% rather than a regexp
  --schema, --dbschema=NUM            patch to schema version NUM
  --schema_type                       Schema type to patch e.g. core|variation|funcgen

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --bindir=DIR                        mysql binary directory (default:
                                      /usr/local/ensembl/bin)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry, --dry_run=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

Please note that where an argument expects a value, this is true for all
alternative argument styles.

=head1 DESCRIPTION

This is a script to facilitate patching databases to the next schema version.
It will connect to a database server and apply schema patches to all databases
where the name matches a pattern. The pattern is a string that can be used in an
'IN' clause in SQL (e.g. "%_core_%"; if you want to patch only a single
database, use a pattern without % expansion, e.g. "homo_sapiens_core_38_36").

If you only want to check which patches need to be applied, use --dry_run=1
(best done in combination with --interactive=0).

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<dev@ensembl.org>

=cut

# Should really add explicit --schema_type param to avoid applying core patches to non-core DBs
# Could also validate this against meta schema.type = core|funcgen|variation

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
$support->param('dbname', undef);
$support->parse_extra_options(
							  'pattern|dbpattern=s',
							  'schema|dbschema=s',
 							  'bindir=s',
							  'schema_type=s',
							 );
my @params = map { $_ unless ($_ =~ /dbname/) } $support->get_common_params;
$support->allowed_params(
  @params,
  'pattern',
  'schema',
  'bindir',
  'schema_type'
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

unless ($support->param('bindir')) {
  $support->param('bindir', '/usr/local/ensembl/bin');
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'pattern',
  'schema',
								'schema_type'
);

my $schema_type = $support->param('schema_type');
my %patch_dirs = (
				  'core' => "$SERVERROOT/ensembl/sql",
				  'funcgen' => "$SERVERROOT/ensembl-functgenomics/sql",
				  'variation' => "$SERVERROOT/ensembl-variation/sql",
				 );

#check schema_type is valid
if(! (defined $schema_type && exists $patch_dirs{$schema_type})){
  $support->log_error('You must specify a valid --schema_type parameter e.g. core|variation|funcgen');
}

# connect to database
my $dbh = $support->get_dbconnection;

# read patches from file
$support->log("Reading patches from file...\n");
my $patchdir = $patch_dirs{$schema_type};
my $schema = $support->param('schema');
my @patches;

opendir(DIR, $patchdir) or
  $support->log_error("Can't opendir $patchdir: $!");

while (my $file = readdir(DIR)) {
  if ($file =~ /^patch_\d+_${schema}.*\.sql$/) {
    $support->log("$file\n", 1);
    push @patches, $file;
  }
}

$support->log("Done.\n\n");

# get all database names that match pattern
my ($sth, $sql);
$sql = "SHOW DATABASES LIKE '".$support->param('pattern')."'";
$sth = $dbh->prepare($sql);
$sth->execute;

# loop over databases
while (my ($dbname) = $sth->fetchrow_array) {
  $support->log_stamped("$dbname\n");
  
  if ($support->user_proceed("\nPatch $dbname?")) {

    # check which patches have already been applied
    $sql = qq(SELECT meta_value FROM $dbname.meta WHERE meta_key = 'patch');
    my $sth1 = $dbh->prepare($sql);
    $sth1->execute;
    
    my %applied;
  
    while (my ($val) = $sth1->fetchrow_array) {
      my ($file) = split(/\|/, $val);
      $applied{$file} = 1;
    }

    # apply the missing ones
    foreach my $patch (sort @patches) {
      $support->log("$patch... ", 1);
      
      if ($applied{$patch}) {
        $support->log("already applied.\n")

      } elsif ($support->param('dry_run')) {
        $support->log("needs applying.\n")
      
      } else {
        $support->log("applying... ");

        my $cmd = $support->param('bindir')."/mysql".
          " -h ".$support->param('host').
          " -P ".$support->param('port').
          " -u ".$support->param('user').
          " -p".$support->param('pass').
          " $dbname < $patchdir/$patch";

        if (system($cmd) == 0) {
          $support->log("done.\n");
        } else {
          $support->log_warning("Error applying patch. Please check patch file.\n", 1);
        }
        
      }
    }
    
  } else {
    $support->log("Skipping on user's request.\n", 1);
  }

  $support->log("Done.\n\n");
}

# finish logfile
$support->finish_log;

