
=pod

=head1 NAME - EnsTestDB

=head1 SYNOPSIS

=head1 DESCRIPTION


=head1 METHODS

=cut

package MultiTestDB;

use vars qw(@ISA);
use strict;

use DBI;
use Digest::MD5;
use Data::Dumper;


#homo sapiens is used if no species is specified
my $DEFAULT_SPECIES  = 'Homo_sapiens';

#configuration file extension appended onto species name
my $FROZEN_CONF_EXT  = '.MultiTestDB.frozen.conf';

my $CONF_FILE    = 'MultiTestDB.conf';

my $DUMP_DIR = 'multi_test_dbs';




sub new {
  my( $pkg, $species ) = @_;

  my $self = bless {}, $pkg;

  unless($species) {
    $species = $DEFAULT_SPECIES;
  }

  $self->species($species);

  if($ENV{'HARNESS_ACTIVE'}) {
    #databases are loaded already, read conf hash from file
    $self->load_config($species);
  } else {
    #load the databases and generate the conf hash
    $self->load_databases($species);
  }

  #generate the db_adaptors from the $self->{'conf'} hash
  $self->create_adaptors;

  return $self;
}
 



#
# load config into $self->{'conf'} hash
#
sub load_config {
  my $self = shift;

  my $conf = $self->species . $FROZEN_CONF_EXT;
    
  eval {
    $self->{'conf'} = do $conf; #reads file into $self->{'conf'}
  };
    
  if($@) {
    die("Could not read frozen configuration file '$conf'\n");
  }
}



#
# Store $self->{'config'} hash into a file
#
sub store_config {
  my $self = shift;

  my $conf = $self->species . $FROZEN_CONF_EXT;

  local *FILE = open ">$conf" or die "Could not open config file '$conf'\n";

  my $string = Dumper($self->{'conf'});

  #strip of leading '$VAR1 = '
  $string =~ s/$[\$]VAR1\s*=//;

  #store config in file
  print FILE $string;
  
  close FILE;
}



#create a set of adaptors based on the $self->{'conf'} hash
sub create_adaptors {
  my $self = shift;

  #establish a connection to each of the databases in the configuration
  foreach my $dbtype (keys %{$self->{'conf'}}) {
    print "Connecting to $dbtype\n";
    my $db = $self->{'conf'}->{$dbtype};
    
    my $adaptor;
    
    #try to instantiate an adaptor for this database 
    eval {
      require $db->{'module'};
      $adaptor = new $db->{'module'}('dbname' => $db->{'name'},
				     'user'   => $db->{'user'},
				     'pass'   => $db->{'pass'},
				     'port'   => $db->{'port'},
				     'host'   => $db->{'host'},
				     'driver' => $db->{'driver'});
    };
	
    if ($@) {
      warn("WARNING: Could not instantiate $dbtype DBAdaptor:\n$@");
    } else {
      $self->{'db_adaptors'}->{$dbtype} = $adaptor;
    }
  }
}




sub load_databases {
  my ($self, $species) = @_;
  
  #create database from conf and from zip files 
  my $db_conf = do $CONF_FILE;

  my $port   = $db_conf->{'port'};
  my $driver = $db_conf->{'driver'};
  my $host   = $db_conf->{'host'};
  my $pass   = $db_conf->{'pass'};
  my $user   = $db_conf->{'user'};
  my $zip    = $db_conf->{'zip'};

  #create a config hash which will be frozen to a file
  $self->{'conf'} = {};

  #unzip database files
  unzip_test_dbs($zip);

  #connect to the database
  my $locator = 'DBI:$driver:host=$host;port=$port';
  my $db = DBI->connect($locator, $user, $pass, {RaiseError => 1});

  unless($db) {
    die "Can't connect to database $locator";
  }

  #create a database for each database specified
  foreach my $dbtype (keys %{$db_conf->{'databases'}}) {
    #create a unique random dbname
    my $dbname = $self->_create_db_name($species, $dbtype);

    
    unless($db->do("CREATE DATABASE $dbname")) {
      die("Could not create database [$dbname]");
    }

    #copy the general config into a dbtype specific config 
    $self->{'conf'}->{$dbtype} = {};
    %{$self->{'conf'}->{$dbtype}} = %$db_conf;

    #store the temporary database name in the dbtype specific config
    $self->{'conf'}->{$dbtype}->{'dbname'} = $dbname;

    $db->do("use $dbname");
    
    #load the database with data
    my $dir = "$DUMP_DIR/species";
    local *DIR;

    opendir(DIR, $dir) or die "could not open dump directory '$dir'";

    my @files = readdir DIR;

    local *FILE;

    #read in table creat statements from *.sql files and process them with DBI
    foreach my $sql_file (grep /\.sql$/, @files) {
      unless(-f $sql_file && -r $sql_file) {
	warn("could not read SQL file '$sql_file'\n");
	next;
      }
      
      FILE = open $sql_file;
      my @file = <FILE>;
      $db->do(join ' ', @file);
      close FILE;

      #import data from the txt files of the same name
      $sql_file  =~ /.*\/(.*)\.sql/;
      my $tablename = $1;
      my $txt_file = s/\.sql$/\.txt/;
      
      unless(-f $txt_file && -r $txt_file) {
	warn("could not read data file '$txt_file'\n");
	next;
      }

      $db->do( "load data local infile '$txt_file' into table $tablename" );
    }	     
  }
  
  closedir DIR;

  $db->disconnect;
  
  #freeze configuration in a file
  $self->store_config;
  

}


sub unzip_test_dbs {
  my ($self, $zipfile) = @_;

  if (-e $DUMP_DIR) {
    warn "Test genome dbs already unpacked\n";
    return;
  }

  unless($zipfile) {
    $self->throw("zipfile argument is required\n");
  }

  unless(-f $zipfile) {
    $self->throw("zipfile could not be found\n");
  }

  system ( "unzip $zipfile" );
}




sub get_DBAdaptor {
  my ($self, $type) = @_;

  unless($type) {
    die('type arg must be specified\n');
  }

  return $self->{'db_adaptors'}->{$type};
}


# convenience method: by calling it, you get the name of the database,
# which  you can cut-n-paste into another window for doing some mysql
# stuff interactively
sub pause {
  my ($self) = @_;
  
  print STDERR "pausing to inspect databases\n";
  foreach my $dbtype (keys %{$self->{'db_adaptors'}}) {
    my $db_adaptor = $self->{'db_adaptors'}->{$dbtype};
    print STDERR " [$dbtype]\n";
    print STDERR "    name=[".$db_adaptor->dbname."]\n";
    print STDERR "    port=[".$db_adaptor->port."]\n";
    print STDERR "    host=[".$db_adaptor->host."]\n";
    print STDERR "    user=[".$db_adaptor->user."]\n";
  }
  print STDERR "press ^D to continue\n";
  `cat `;
}



sub species {
  my ($self, $species) = @_;

  if($species) {
    $self->{'species'} = $species;
  }

  return $self->{'species'};
}



sub _create_db_name {
    my( $self, $species, $dbtype ) = @_;

    my $rand = &Digest::MD5::md5_hex(rand());
    my $db_name = "_test_db_${species}_${dbtype}_${rand}";

    return $db_name;
}




sub do_sql_file {
  my( $self, @files ) = @_;
  local *SQL;
  my $i = 0;
  my $dbh = $self->db_handle;
  
  my $comment_strip_warned=0;

  foreach my $file (@files) {
    my $sql = '';
    open SQL, $file or die "Can't read SQL file '$file' : $!";
    while (<SQL>) {
      # careful with stripping out comments; quoted text
      # (e.g. aligments) may contain them. Just warn (once) and ignore
      if (    /'[^']*#[^']*'/ 
	|| /'[^']*--[^']*'/ ) {
	  if ( $comment_strip_warned++ ) { 
	    # already warned
	  } else {
	    warn "#################################\n";
	    warn "# found comment strings inside quoted string;" .
	         "not stripping, too complicated: $_\n";
	    warn "# (continuing, assuming all these they are simply " .
	      "valid quoted strings)\n";
	    warn "#################################\n";
	  }
	} else {
	  s/(#|--).*//;       # Remove comments
	}
        next unless /\S/;   # Skip lines which are all space
        $sql .= $_;
        $sql .= ' ';
      }
  close SQL;
        
    #Modified split statement, only semicolumns before end of line,
    #so we can have them inside a string in the statement
    #\s*\n, takes in account the case when there is space before the new line
    foreach my $s (grep /\S/, split /;[ \t]*\n/, $sql) {
      $s =~ s/\;\s*$//g;
      $self->validate_sql($s);
      $dbh->do($s);
      $i++
    }
  }
  return $i;
}                                       # do_sql_file

sub validate_sql {
  my ($self, $statement) = @_;
  if ($statement =~ /insert/i) {
    $statement =~ s/\n/ /g; #remove newlines
    die ("INSERT should use explicit column names " .
	 "(-c switch in mysqldump)\n$statement\n")
      unless ($statement =~ /insert.+into.*\(.+\).+values.*\(.+\)/i);
  }
}



sub DESTROY {
    my( $self ) = @_;

}



1;

