#!/usr/local/bin/perl

# script to manage db_history table in mysql which records
# history of database

use strict;
use Getopt::Long;
use DBI;
use Sys::Hostname;

# hard wired
my $driver="mysql";
my $port=3306;
my $db_history="db_history";
my $create_user;       # username for read/write access
my $create_password;   # password for read/write access
my $user;              # username for readonly access
my $password;          # password for readonly access

# derived
# name of user running db_history.pl
my ($uname)=getpwuid($<);

# user restrictions
my $host;
my $cluster;
my $db;
my $label;
my $full_history;

# count action
my $table;

# insert values into database
my $set_owner;
my $set_user;
my $set_title;
my $set_note;
my $set_label;

my $opt_C;

# help
my $help;
my $phelp;

# options not listed in help: mainly testing
my $all_databases;
my $opt_v;

$Getopt::Long::ignorecase=0;

&GetOptions(

	    'host|H:s'=>\$host,
	    'cluster|r:s'=>\$cluster,
	    'db|d:s'=>\$db,
	    'label|l:s'=>\$label,

	    'full_history|F'=>\$full_history,

	    'table|t:s'=>\$table,

	    'set_owner|O:s'=>\$set_owner,
	    'set_user|U:s'=>\$set_user,
	    'set_title|T:s'=>\$set_title,
	    'set_note|N:s'=>\$set_note,
	    'set_label|L:s'=>\$set_label,
	    'create|C'=>\$opt_C,

	    'help'=>\$phelp,
	    'h'=>\$help,

	    'user|u:s'=>\$user,
	    'all_databases|A'=>\$all_databases,
	    'v'=>\$opt_v,
	    'passwd|password|p:s'=>\$password,
	    'create_password|P:s'=>\$create_password,
	    'create_user|c:s'=>\$create_user,

	     );

# help
if($phelp){
    exec('perldoc', $0);
    exit 0;
}
if($help){
    print<<ENDOFTEXT;
db_history.pl

  -H|host         host    restrict to this host
  -r|cluster      name    restrict to databases on this cluster
  -d|db           dbname  restrict to this database
  -l|label        label   restrict to databases with this label set

  -F|full_history         show full history of each database

  -t|table     string     count entries in this table across selected databases

  -O|set_owner id         identify the owner/creator of the database
  -U|set_user  id         identify a user who needs this database
  -T|set_title string     short description of database
  -N|set_note  string     append a note about the database
  -L|set_label string     add a label to a database, e.g identifying it as part of a set

  -C|create               create table in database for history

  -h                      this help
  -help                   perldoc help

ENDOFTEXT
    exit 0;
}

# databases to ignore
my %skip_db=('test'=>1,
	     'mysql'=>1,
	     );

# set hosts to process
my @hosts;
my @ecs1=('ecs1a',
	  'ecs1b',
	  'ecs1c',
	  'ecs1d',
	  'ecs1e',
	  'ecs1f',
	  );
my @ecs2=('ecs2a',
	  'ecs2b',
	  'ecs2c',
	  'ecs2d',
	  'ecs2e',
	  'ecs2f',
	  );

if ($host) {
    @hosts=($host);
}
elsif($cluster) {
    if ($cluster eq 'ecs1') {
	@hosts = @ecs1;
    }
    elsif ($cluster eq 'ecs2') {
	@hosts = @ecs2;
    }
    else {
	die "Unknown cluster $cluster";
    }
}
unless(@hosts){@hosts=(@ecs1, @ecs2);}

# mode 1: write data
if($db && $host){
    if($set_owner || $set_user || $set_title || $set_note || $set_label){
	&_db_history_add($host,$db,
			 $set_owner,$set_user,
			 $set_title,$set_note,$set_label);
	exit 0;
    }

    # mode 2: detailed report about a particular database
    # if databast and host is given
    my @history;
    my $err=&_db_history($host,$db,\@history);

    if($err<-1){
	print "$db not found on $host\n";
	exit 0;
    }elsif($err<0){
	print "No history for $db on $host\n";
	exit 0;
    }

    # report full history
    my $info;
    &_summarize_history(\@history,\$info);
    print "$info\n";

    exit 0;
}

# mode 3: summary of all databases (count of table if -table)
# applying restrictions by label and or host

# get info about each database
my $description;
if($table){
    $description="Rows in Table \'$table\'";
}elsif($full_history){
    $description="";
}else{
    $description="TITLE";
}
printf "\n%-10s %-30s %s\n",'HOST','DATABASE',$description;
print "-" x 70 ."\n";

foreach my $host (@hosts){
    print "checking $host\n" if $opt_v;

    # get list of databases
    my @databases;
    if($db){
	@databases=($db);
    }else{
	&_db_list($host,\@databases);
	print "found ".scalar(@databases)." databases on $host\n" if $opt_v;
    }

    foreach my $database (@databases){

	# unless explicitly told to show all databases, skip ones
	# in skip list
	next if(!$all_databases && $skip_db{$database});

	if($table){
	    my $count;
	    my $err=&_db_count($host,$database,$table,\$count);
	    # no database
	    next if $err<0;
	    # no table
	    next if $count<0;

	    $count=0 if !$count;
	    printf "%-10s %-30s %10d\n",$host,$database,$count;
	}elsif($full_history){
	    my @history;
	    my $err=&_db_history($host,$database,\@history);
	    my $info;
	    my $err2=&_summarize_history(\@history,\$info);
	    if($err2<-2){
		# label not found - skip
		next;
	    }else{
		printf "%-10s %-30s\n",$host,$database;
		if($err<0){
		    print " No history\n\n";
		}else{
		    print "$info\n";
		}
	    }
	}else{
	    my $title;
	    my $err=&_db_history_title($host,$database,\$title);
	    # skip if database doesn't exist or missing label; continue if table exists
	    next if $err<-1;

	    printf "%-10s %-30s \'%s\'\n",$host,$database,$title;
	}
    }
}

# report summary

# get list of db's for a host
sub _db_list{
    my($host,$radatabases)=@_;
    #$host='';
    eval{
	my $dsn = "DBI:$driver:database=$db;host=$host;port=$port";
	my $dbh = DBI->connect($dsn, $user, $password);
	@$radatabases = $dbh->func('_ListDBs');
	$dbh->disconnect();
    };
    if($@){
	print "failed to fetch list of databases on $host\n$@\n" if $opt_v;
    }
}

sub _summarize_history{
    my($rahistory,$rinfo)=@_;
    my($notes,$title,$owner,@users,@labels,%labels);
    foreach my $row (@$rahistory){
	my($type,$user,$date,$isc,$text)=@$row;
	if($type eq 'owner' && $isc==1){
	    $owner=$text;
	}elsif($type eq 'title' && $isc==1){
	    $title=$text;
	}elsif($type eq 'note'){
	    my $ymd=&_time2date($date);
	    $notes.=sprintf("   %s %-10s %s\n",$ymd,$user,$text);
	}elsif($type eq 'user'){
	    push(@users,$text);
	}elsif($type eq 'label'){
	    push(@labels,$text);
	    $labels{$text}=1;
	}
    }
    if($label && !$labels{$label}){
	return -3;
    }
    my $info;
    if(!$owner && !$title && !$notes && scalar(@users)==0 && scalar(@labels)==0){
	$info="no history data";
    }else{
	$info=" TITLE: \'$title\'\n OWNER: \'$owner\', USERS: \'".join('\',\'',@users)."\'\n";
	if(scalar(@labels)){
	    $info.=" LABELS: \'".join('\',\'',@labels)."\'\n";
	}
	$info.=$notes;
    }
    $$rinfo=$info;
    return;
}

sub _time2date{
    my($date)=@_;
    if($date=~/(\d+\-\d+\-\d+)/){
	return $1;
    }else{
	return 'XXXX-XX-XX';
    }
}

sub _db_history_add{
    my($host,$database,
       $set_owner,$set_user,
       $set_title,$set_note,$set_label)=@_;

    # create history table if it doesn't exist
    &_db_history_create($host,$database) if $opt_C;

    my $dbh;
    if(my $err=&_db_connect(\$dbh,$host,$database,$create_user,$create_password)){return $err};
    &_db_history_insert($dbh,$host,$database,'owner',$set_owner) if $set_owner;
    &_db_history_insert($dbh,$host,$database,'user',$set_user) if $set_user;
    &_db_history_insert($dbh,$host,$database,'title',$set_title) if $set_title;
    &_db_history_insert($dbh,$host,$database,'note',$set_note) if $set_note;
    &_db_history_insert($dbh,$host,$database,'label',$set_label) if $set_label;
    $dbh->disconnect();
}

# insert an entry into table
sub _db_history_insert{
    my($dbh,$host,$database,$type,$string)=@_;
    print "adding $type of value $string\n" if $opt_v;
    eval{
	# mark all existing rows as not current
	my $query="update $db_history set ";
	$query.="is_current=0 where ";
	$query.="type=".$dbh->quote($type);
	$dbh->do($query);

	# add new row, marking as current
	my $query="insert into $db_history set ";
	$query.="type=".$dbh->quote($type).", ";
	$query.="user=".$dbh->quote($uname).", ";
	$query.="date=now(), is_current=1, text=".$dbh->quote($string);
	$dbh->do($query);
    };
    if($@){
	my $old_err=$@;
	
	# if failed, maybe because table didn't exist - check
	my @history;
	my $err=&_db_history($host,$db,\@history);
	if($err){
	    print "$db on $host does not have a history table - create using '-C'\n";
	    $dbh->disconnect();
	    exit 0;
	}
	print "failed to change $type\n$@\n" if $opt_v;
	return -1;
    }
}

# create table
sub _db_history_create{
    my($host,$database)=@_;
    my $dbh;
    my $local_host=hostname();
    my $connect_host;
    if($local_host eq $host){
	$connect_host="";
    }else{
	$connect_host=$host;
    }
    if(my $err=&_db_connect(\$dbh,$connect_host,$database,$create_user,$create_password)){return $err};
    eval{
	my $query="create table if not exists $db_history".
	    "(type enum ('owner','user','title','note','label') not null,".
		"user varchar(10), date datetime default '0000-00-00 00:00:00' not null, ".
		    "is_current tinyint(1), text varchar(80) )";
	$dbh->do($query);
    };
    if($@){
	print "failed to create table $db_history in $database on $host\n$@\n" if $opt_v;
	return -1;
    }
    $dbh->disconnect();
}

# 2 possible conditions
# - no database with this name
# - no table with this name
sub _db_history{
    my($host,$database,$rahistory)=@_;
    my $dbh;
    if(my $err=&_db_connect(\$dbh,$host,$database,$user,$password)){return $err};

    eval{
	my $sth = $dbh->prepare("SELECT * FROM $db_history order by date");
	$sth->execute;
	while(my @array=$sth->fetchrow_array()){
	    push(@$rahistory,[@array]);
	}
	$sth->finish;
    };
    $dbh->disconnect();
    if($@){
	print "failed to fetch entries in $db_history from $database on $host\n$@\n" 
	    if $opt_v;
	return -1;
    }
}

sub _db_history_title{
    my($host,$database,$rtitle)=@_;
    my $dbh;
    if(my $err=&_db_connect(\$dbh,$host,$database,$user,$password)){return $err};

    # look for a title
    eval{
	my $sth = $dbh->prepare("SELECT text FROM $db_history ".
				"where type='title' and is_current=1");
	$sth->execute;
	($$rtitle)=$sth->fetchrow_array();
	if($$rtitle eq ""){$$rtitle="no title";}
	$sth->finish;
    };

    # no title - failed due to missing table
    if($@){
	$dbh->disconnect();
	if($label){
	    # BUG - assumes that if no title, won't have label
	    print "failed to fetch title - skip since selection by label\n$@\n" 
		if $opt_v;
	    return -3;
	}else{
	    print "failed to fetch title\n$@\n" 
		if $opt_v;
	    $$rtitle="no title";
	    return -1;
	}
    }

    # if got a title and label, check label
    if($label){
	my $count;
	eval{
	    my $sth = $dbh->prepare("SELECT count(*) FROM $db_history ".
				    "where type='label' and text='$label'");
	    $sth->execute;
	    $count=$sth->fetchrow_array();
	    $sth->finish;
	};
	if($@ || $count==0){
	    print "$database does not have label $label\n$@\n" 
		if $opt_v;
	    $dbh->disconnect();
	    return -3;
	}
    }
    $dbh->disconnect();
}

# return row count for particular table
sub _db_count{
    my($host,$database,$table,$rcount)=@_;
    my $dbh;
    if(my $err=&_db_connect(\$dbh,$host,$database,$user,$password)){return $err};

    # try to access
    eval{
	my $sth = $dbh->prepare("SELECT count(*) FROM $table");
	$sth->execute;
	$$rcount=$sth->fetchrow_array();
	$sth->finish;
    };
    $dbh->disconnect();
    if($@){
	print "failed to count rows in $table from $database on $host\n$@\n" if $opt_v;
	$$rcount=-1;
    }
}

# connect to db with error handling
sub _db_connect{
    my($rdbh,$host,$database,$user,$password)=@_;
    my $dsn = "DBI:$driver:database=$database;host=$host;port=$port";

    # try to connect to database
    eval{
	$$rdbh = DBI->connect($dsn, $user, $password,
			    { RaiseError => 1, PrintError => 0 });
    };
    if($@){
	print "$database not on $host\n$@\n" if $opt_v;
	return -2;
    }
}

__END__

=pod

=head1 db_history.pl

=head1 DESCRIPTION

Tool to annotate mysql databases with a descriptive title, labels,
users and other notes.  Allows this information to be retrivied across
all mysql databases on acari cluster.

=head1 EXAMPLES

Get titles of all databases:

    db_history.pl


Set of databases can be restricted by host (-H host), by database name
(-d database) and by label [user added] (-l ensembl100).

e.g. Find location of all copies of database 'ensembl100':

    db_history.pl -d new_genes


Various types of annotation can be added to a database using (-O
owner, -U user, -T title, -N notes, -L label).  The proposed usage for these field are:

title = one line title that everyone can understand
owner = person responsible for database
user  = others who have an interest in the database not being deleted(!)
notes = single line coments adding to the status of the database (timestamped)
label = subsets of databases can be selected by label, so this should allow databases
        required for a particular release to be labelled, e.g. 'ensembl100'.

e.g. Set the title of the new_genes database on host ensembl1, the owner, an interested user and
a note about it being incomplete.

    db_history.pl -d new_genes -H hostx -O jack -U jill 
        -T 'database of extra genes from RTpcr experiments'
        -N 'only genes on chr1 have been loaded'


Full annotation of individual databases can be viewed (notes are viewed in historical order):

    db_history.pl -d new_genes -H ens1d

Full history of a list of database can be shown using the -F flag

    db_history.pl -d new_genes -F


As well as return user added notes about the databases, its also possible to compare
the number of rows in a table across databases:

e.g. Count number of rows in gene table across all databases:

    db_history.pl -t gene




Count

=head1 FLAGS

=over 4

=item -h

Displays short help

=item -help

Displays this help message

=back

=head1 VERSION HISTORY

=over 4

=item 09-JUL-2001

B<th> released first version

=back

=head1 BUGS

=head1 AUTHOR

B<Tim Hubbard> Email th@sanger.ac.uk

=cut
