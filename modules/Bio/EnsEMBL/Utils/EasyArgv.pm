=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 AUTHOR

Juguang Xiao <juguang@tll.org.sg>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::EasyArgv

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::EasyArgv;

  my $db = get_ens_db_from_argv;    # this method is exported.

  use Getopt::Long;

  my $others;
  &GetOptions( 'others=s' => \$others );

=head1 DESCRIPTION

This is a lazy but easy way to get the db-related arguments. All you
need to do is to invoke get_ens_db_from_argv before using standard
Getopt. The below options will be absorbed and removed from @ARGV.

db_file, host, db_host, dbhost, user, db_user, dbuser, pass, db_pass,
dbpass, dbname, db_name.

Now you can take advantage of Perl's do method to execute a file as perl
script and get returned the last line of it. For your most accessed db
setting, you can have a file named, say, ensdb_homo_core_18.perlobj,
with the content like

  use strict;    # The ceiling line

  use Bio::EnsEMBL::DBSQL::DBAdaptor;

  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host   => 'ensembldb.ensembl.org',
    -user   => 'anonymous',
    -dbname => 'homo_sapiens_core_18_34'
  );

  $db;           # The floor line

In the your command line, you just need to write like 

  perl my_script.pl -db_file ensdb_homo_core_18.perlobj

rather than the verbose

  -host ensembldb.ensembl.org -user anonymous \
  -dbname homo_sapiens_core_18_34

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::EasyArgv;

use strict;
use vars qw($debug);
use Exporter ();
our @ISA= qw(Exporter);
our @EXPORT = qw(get_ens_db_from_argv
);
use Bio::Root::Root; # For _load_module
use Getopt::Long;

sub _debug_print;

sub get_ens_db_from_argv {
    my ($db_file, $host, $user, $pass, $dbname, $driver, $db_module);
    $host = 'localhost';
    $driver ='mysql';
    $db_module = 'Bio::EnsEMBL::SQL::DBAdaptor';
    Getopt::Long::config('pass_through');
    &GetOptions(
        'db_file=s' => \$db_file,
        'driver|dbdriver|db_driver=s' => \$driver,
        'host|dbhost|db_host=s' => \$host,
        'user|dbuser|db_user=s' => \$user,
        'pass|dbpass|db_pass=s' => \$pass,
        'dbname|db_name=s' => \$dbname,
        'db_module=s' => \$db_module
    );

    my $db;
    if(defined $db_file){
        -e $db_file or die "'$db_file' is defined but does not exist\n";
        eval { $db = do($db_file) };
        $@ and die "'$db_file' is not a perlobj file\n";
        $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')
            or die "'$db_file' is not EnsEMBL DBAdaptor\n";
        _debug_print "I get a db from file\n";
        
    }elsif(defined $host and defined $user and defined $dbname){
        Bio::Root::Root::_load_module($db_module);
        $db = $db_module->new(
            -host => $host,
            -user => $user,
            -pass => $pass,
            -dbname => $dbname,
            -driver => $driver
        );
    }else{
        die "Cannot get the db, due to the insufficient information\n";
    }
    return $db;
}

sub _debug_print {
    print STDERR @_ if $debug;
}


1;
