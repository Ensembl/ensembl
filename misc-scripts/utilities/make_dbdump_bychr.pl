#This script generates a full dump of the EnsEMBL database for 
#a particular chromosome. Useful to create a small but fully functional
#EnsEMBL db

#!/usr/local/bin/perl

=head1 NAME

make_dbdump_bychr

=head1 SYNOPSIS

  1)Needs to be called within a new directory where you want
  all the files to be written
  2)with a user that is allowed to use mysqldump
  3)needs to be run on the host that runs the daemon

  make_dbdump_bychr -chr 'chr22'

=head1 DESCRIPTION

This script generates a full dump of the EnsEMBL database for 
a particular chromosome. Useful to create a small but fully functional
EnsEMBL db (e.g. laptop mini-mirror)
=cut

use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $host = "localhost";
my $port   = '';
my $dbname = 'idmap_apr01';
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $chr = 'chr22';

&GetOptions( 
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'chr:s'      => \$chr
	     );

my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
$db =  Bio::EnsEMBL::DBLoader->new($locator);

#Start from schema
my $command = "mysqldump -u $dbuser -d $dbname > table.all";
system($command);

#Then dump data from tables that are needed in full regardless of
#dump size: analysis,analysisprocess,chromosome,externalDB,meta,species

$command = ("mysqldump -u $dbuser -T . $dbname analysis analysisprocess chromosome externalDB meta species");
system ($command);

my $sth = $db->prepare("select * from static_golden_path where chr_name = '$chr' limit 1");

$sth->execute;
my @contigs;
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;

    #Relying on raw_id being 3rd element in the table!
    push (@contigs,$array[2]);
    open (FILE,">static_golden_path.txt");
    print FILE join("\t",@array)."\n";
    close (FILE);
}

#Now get partial dumps of all tables that relate to contigs
#using the @contigs array


