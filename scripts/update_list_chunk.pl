#!/usr/local/bin/perl -w

=head1 NAME

Update

=head1 SYNOPSIS
 
  update.pl

=head1 DESCRIPTION

This script updates a recipient database by checking its donor database

=head1 OPTIONS

    -thost     host name for database (gets put as host= in locator)

    -tport     For RDBs, what port to connect to (port= in locator)

    -tdbname   For RDBs, what name to connect to (dbname= in locator)

    -tdbuser   For RDBs, what username to connect as (dbuser= in locator)

    -tpass     For RDBs, what password to use (dbpass= in locator)

    -help      Displays script documentation with PERLDOC
    
    -nowrite   Runs entire script without writing in recipient

    -verbose   Gets all the print STDERR for testing purposes

=cut

use Bio::EnsEMBL::Analysis::UpdateManager;

use strict;
use Getopt::Long;
use vars qw(@ISA);



@ISA = qw(Bio::Root::Object);

my $tdbtype = 'rdb';
my $thost   = 'obi-wan.sanger.ac.uk';
my $tport   = '410000';
my $tdbname = 'ensembl';
my $tdbuser = 'root';
my $tpass = undef;
my $adbname = 'ens_archive';
my $arcpass = undef;
my $module = "Bio::EnsEMBL::DBSQL::Obj";
my $archive = 0;

my $freeze;
my $help;
my $nowrite;
my $verbose = 1;
my $slice;
my $usefile = 0;
#!/usr/local/bin/perl -w

=head1 NAME

Update

=head1 SYNOPSIS
 
  update.pl

=head1 DESCRIPTION

This script updates a recipient database by checking its donor database

=head1 OPTIONS

    -thost     host name for database (gets put as host= in locator)

    -tport     For RDBs, what port to connect to (port= in locator)

    -tdbname   For RDBs, what name to connect to (dbname= in locator)

    -tdbuser   For RDBs, what username to connect as (dbuser= in locator)

    -tpass     For RDBs, what password to use (dbpass= in locator)

    -help      Displays script documentation with PERLDOC
    
    -nowrite   Runs entire script without writing in recipient

    -verbose   Gets all the print STDERR for testing purposes

=cut

use Bio::EnsEMBL::Analysis::UpdateManager;

use strict;
use Getopt::Long;
use vars qw(@ISA);



@ISA = qw(Bio::Root::Object);

my $tdbtype = 'rdb';
my $thost   = 'obi-wan.sanger.ac.uk';
my $tport   = '410000';
my $tdbname = 'ensembl';
my $tdbuser = 'root';
my $tpass = undef;
my $adbname = 'ens_archive';
my $arcpass = undef;
my $module = "Bio::EnsEMBL::DBSQL::Obj";
my $archive = 0;

my $help;
my $nowrite;
my $verbose = 1;
my $slice;
my $usefile = 0;

my $from;

&GetOptions( 
	     'tdbtype:s'  => \$tdbtype,
	     'thost:s'    => \$thost,
	     'archive'    => \$archive,
	     'tport:n'    => \$tport,
	     'tdbname:s'  => \$tdbname,
	     'tdbuser=s'  => \$tdbuser,
	     'tpass:s'    => \$tpass,
	     'arcpass:s'  => \$arcpass,
             'adbname:s'  => \$adbname,
	     'module=s'  => \$module,
	     'usefile=s' => \$usefile,
	     'h|help'    => \$help,
	     'nowrite'   => \$nowrite,
	     'slice:s'   => \$slice,
	     'v|verbose' => \$verbose,
	     'from:n'    => \$from,
	     );

if ($help) {
    exec('perldoc', $0);
}

$|=1;

my $to_locator       = make_locator_string($tdbtype,$module,$thost,$tport,$tdbname,$tdbuser,$tpass);
my $tdb              = new Bio::EnsEMBL::DBLoader($to_locator);
my $tdb_update_obj      = Bio::EnsEMBL::DBSQL::Update_Obj->new($tdb);
#my $from_locator     = "Bio::EnsEMBL::TimDB::Obj";
my $from_locator     = $tdb->get_donor_locator;
my $arc_locator;
if ($archive) {
    $arc_locator = "Bio::EnsEMBL::DBArchive::Obj//host=$thost;port=$tport;dbname=$adbname;user=$tdbuser;pass=$arcpass";
}
else {
    $arc_locator = "none";
}
my $last_offset;
if($from){
    $last_offset=$from;
}else{
    $last_offset=$tdb_update_obj->get_last_update_offset;
}
my $now_offset       = $tdb_update_obj->get_now_offset;    # This should be something different

print STDERR "From time $last_offset\n";
print STDERR "To time: $now_offset\n" if ($now_offset);
print "Trying output... verbose=$verbose\n";

$| = 1;

if ($last_offset > $now_offset) {
    print "Time of last_offset update more recent than now-offset, exiting!\n";
    exit;
}

my $update_manager   = new Bio::EnsEMBL::Analysis::UpdateManager(-fromlocator => $from_locator,
								 -tolocator   => $to_locator,
								 -arclocator  => $arc_locator,
								 -fromtime    => $last_offset,
								 -totime      => $now_offset,
								 );

print "USEFILE = $usefile\n";
$update_manager->nowrite  ($nowrite);
$update_manager->verbose  ($verbose);
$update_manager->usefile ($usefile);
$update_manager->chunksize(10);
$update_manager->update;


sub make_locator_string {
    my ($type,$module,$host,$port,$dbname,$dbuser,$dbpass) = @_;

    if ($type eq "rdb") {
	return 	"$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
    } elsif ($type eq "timdb") {
	return "Bio::EnsEMBL::TimDB::Obj";
    } else {
	die "Database type [$type] not recognised\n";
    }
}



my $from;

&GetOptions( 
	     'tdbtype:s'  => \$tdbtype,
	     'thost:s'    => \$thost,
	     'archive'    => \$archive,
	     'tport:n'    => \$tport,
	     'tdbname:s'  => \$tdbname,
	     'tdbuser=s'  => \$tdbuser,
	     'tpass:s'    => \$tpass,
	     'arcpass:s'  => \$arcpass,
             'adbname:s'  => \$adbname,
	     'module=s'  => \$module,
	     'usefile=s' => \$usefile,
	     'h|help'    => \$help,
	     'freeze'    => \$freeze,
	     'nowrite'   => \$nowrite,
	     'slice:s'   => \$slice,
	     'v|verbose' => \$verbose,
	     'from:n'    => \$from,
	     );

if ($help) {
    exec('perldoc', $0);
}

$|=1;

my $to_locator       = make_locator_string($tdbtype,$module,$thost,$tport,$tdbname,$tdbuser,$tpass);
my $tdb              = new Bio::EnsEMBL::DBLoader($to_locator);
my $tdb_update_obj      = Bio::EnsEMBL::DBSQL::Update_Obj->new($tdb);
#my $from_locator     = "Bio::EnsEMBL::TimDB::Obj";
my $from_locator     = $tdb->get_donor_locator;
my $arc_locator;
if ($archive) {
    $arc_locator = "Bio::EnsEMBL::DBArchive::Obj//host=$thost;port=$tport;dbname=$adbname;user=$tdbuser;pass=$arcpass";
}
else {
    $arc_locator = "none";
}
my $last_offset;
if($from){
    $last_offset=$from;
}else{
    $last_offset=$tdb_update_obj->get_last_update_offset;
}
my $now_offset       = $tdb_update_obj->get_now_offset;    # This should be something different

print STDERR "From time $last_offset\n";
print STDERR "To time: $now_offset\n" if ($now_offset);
print "Trying output... verbose=$verbose\n";

$| = 1;

if ($last_offset > $now_offset) {
    print "Time of last_offset update more recent than now-offset, exiting!\n";
    exit;
}

my $update_manager   = new Bio::EnsEMBL::Analysis::UpdateManager(-fromlocator => $from_locator,
								 -tolocator   => $to_locator,
								 -arclocator  => $arc_locator,
								 -fromtime    => $last_offset,
								 -totime      => $now_offset,
								 );

print "USEFILE = $usefile\n";
$update_manager->nowrite  ($nowrite);
$update_manager->verbose  ($verbose);
$update_manager->usefile ($usefile);
if $freeze {
    $update_manager->freeze(2);
    $update_manager->nogene(1);
}
$update_manager->chunksize(10);
$update_manager->update;


sub make_locator_string {
    my ($type,$module,$host,$port,$dbname,$dbuser,$dbpass) = @_;

    if ($type eq "rdb") {
	return 	"$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
    } elsif ($type eq "timdb") {
	return "Bio::EnsEMBL::TimDB::Obj";
    } else {
	die "Database type [$type] not recognised\n";
    }
}


