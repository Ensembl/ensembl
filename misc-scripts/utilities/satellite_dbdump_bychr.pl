#!/usr/local/bin/perl
# 
# $Id$
# 
#

=head1 NAME

 sattelite_dbdump_bychr

=head1 SYNOPSIS

  This script generates a dump of an EnsEMBL satellite database for
  particular chromosome. Useful to create a small but fully functional
  ensembl installation, e.g. for a laptop. It needs access to an ensembl-lite
  database (for golden path etc.)

  (1) Needs to be called within a new directory where you want
     all the files to be written

  (2) with a user that is allowed to use mysqldump

  (3) needs to be run on the host that runs the daemon

  (4) Usage: 

       satellite_dbdump_bychr  -<dbtype> <dbinstance>
  
     e.g
  
       satellite_dbdump_bychr  -disease homo_sapiens_disease_110

     Known types are: embl est expression family maps snp

=head1 DESCRIPTION

This script generates a full dump of one or several EnsEMBL sattelite
database for a particular chromosome. Useful to create a small but fully
functional EnsEMBL db (e.g. laptop mini-mirror) 

Based on make_dbdumpk_bychr (which does the core database). 

=cut

;

use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $workdir = `pwd`; chomp($workdir);
my $host = "localhost";
my $port   = '';
my $litedb = ''; # 'homo_sapiens_lite_110'; # force user to provide it
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $chr = 'chr21';                      # smaller than chr22
# my $lim;
my $mysql = 'mysql'; 
my $mysqldump = 'mysqldump'; # in $PATH we trust
#                  /mysql/current/bin/mysqldump

# satellites:
my $famdb;
my $diseaseb;
# end of satellites

&GetOptions( 
            'port:n'     => \$port,
            'litedb:s'   => \$litedb,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'module:s'   => \$module,
            'chr:s'      => \$chr,
            'workdir:s'  => \$workdir,
            'limit:n'    => \$lim,
            'family:s' => \$famdb,
            'disease:s' => \$diseasedb,
	     );

die "need a litedb; use -litedb something " unless $litedb;
die "chromosome names should start with 'chr'" unless $chr =~ /^chr/;
my $pass_arg=""; $pass_arg="-p$dbpass" if $dbpass;

my $limit;
if ($lim) {
    $limit = "limit $lim";
}

if ($famdb) {
    my $dumpdir = "$workdir/$famdb";
    dump_schema($famdb, $dumpdir, 'family.sql');

    my $sql;
    $sql = "
SELECT distinct f.* 
FROM $famdb.family f, $litedb.gene g
WHERE g.chr_name = '$chr'
  and g.family = f.id
  $limit
";
    dump_data($litedb, $sql, $dumpdir, 'family.dat' );

    $sql = "
SELECT fm.* 
FROM $famdb.family_members fm, $famdb.family f, $litedb.gene g
WHERE g.chr_name = '$chr'
  and g.family = f.id
  and f.internal_id  = fm.family
  $limit
";
    dump_data($litedb, $sql, $dumpdir, 'family_members.dat' );
}

if ($diseasedb) {
    my $dumpdir = "$workdir/$diseasedb";
    dump_schema($diseasedb, $dumpdir, 'disease.sql');

# may need an ALTER TABLE gene ADD KEY(gene_symbol);
    my $sql;
    $sql = "
SELECT dg.*
FROM  $diseasedb.gene dg, 
      $litedb.gene lg, 
      $litedb.gene_xref lgx
WHERE lg.chr_name = '$chr' 
  AND lg.gene = lgx.gene 
  AND lgx.display_id = dg.gene_symbol
";
    dump_data($litedb, $sql, $dumpdir, 'gene.dat' );

    $sql = "
SELECT dd.*
FROM  $diseasedb.gene dg, 
      $diseasedb.disease dd,
      $litedb.gene lg, 
      $litedb.gene_xref lgx
WHERE lg.chr_name = '$chr' 
  AND lg.gene = lgx.gene 
  AND lgx.display_id = dg.gene_symbol
  AND dd.id = dg.id;
";
    dump_data($litedb, $sql, $dumpdir, 'disease.dat' );

# here's the sql to restrict the disease_index_*list, but they're so small
# it's really not worth the trouble. Left here in case anyone is interested
#     $sql = "
# SELECT ddl.*
# FROM  $diseasedb.gene dg, 
#       $diseasedb.disease_index_doclist ddl,
#       $litedb.gene lg, 
#       $litedb.gene_xref lgx
# WHERE lg.chr_name = '$chr' 
#   AND lg.gene = lgx.gene 
#   AND lgx.display_id = dg.gene_symbol
#   AND ddl.id  = dg.id
# ";

    foreach my $w ( qw(doc stop vector word) ) {
        my $table = "disease_index_${w}list";
        $sql = "select * from $diseasedb.$table";
        dump_data($litedb, $sql, $dumpdir, "$table.dat" );
    }
}

sub dump_schema {
    my ($dbinstance, $destdir, $destfile) = @_;

    unless (-d $destdir) {
        mkdir $destdir, 0755 || die "mkdir $destdir: $!";
    }

    my $d = "$destdir/$destfile";

    warn "Dumping database schema of $dbinstance to $d\n";
    die "$d exists" if -s $d ;
    $command = "$mysqldump -u $dbuser $pass_arg -d $dbinstance > $d ";
    if ( system($command) ) {
        die "Error: ``$command'' ended with exit status $?";
    }
}

sub dump_data {
    my($db, $sql, $destdir, $destfile) = @_;

    unless (-d $destdir) {
        mkdir $destdir, 0755 || die "mkdir $destdir: $!";
    }
    
    $sql =~ s/\s+/ /g;
    
    my $cmd = "echo \"$sql\" | $mysql -q --batch -u $dbuser -p$dbpass $db > $destdir/$destfile";
    warn "dumping: $cmd\n";

    if ( system($cmd) ) { 
        die "``$cmd'' exited with exit-status $?";
    }
}

## stuff below is not used (yet), since everything is done by plain SQL

## This comes from family-input.pl, and should at one point be put somewhere
## more central (the ones in EnsEMBL load modules etc. that are not relevant)
## Takes string that looks like
## "database=foo;host=bar;user=jsmith;passwd=secret", connects to mysql
## and return the handle
sub db_connect { 
    my ($dbcs) = @_;

    my %keyvals= split('[=;]', $dbcs);
    my $user=$keyvals{'user'};
    my $paw=$keyvals{'pass'};
#    $dbcs =~ s/user=[^;]+;?//g;
#    $dbcs =~ s/password=[^;]+;?//g;
# (mysql doesn't seem to mind the extra user/passwd values, leave them)

    my $dsn = "DBI:mysql:$dbcs";

    my $dbh=DBI->connect($dsn, $user, $paw) ||
      die "couldn't connect using dsn $dsn, user $user, password $paw:" 
         . $DBI::errstr;
    $dbh->{RaiseError}++;
    $dbh;
}                                       # db_connect

sub unique {
    
    my @unique;
    my %seen = ();
    foreach my $item (@_) {
	push(@unique,$item) unless $seen{$item}++;
    }
    return @unique;
}

sub get_inlist {
    my $string_flag = shift (@_);
    my $string;
    foreach my $element (@_) {
	if ($string_flag) {
	    $string .= "'$element',";
	}
	else {
	    $string .= "$element,";
	}
    }
    $string =~ s/,$//;
    return "($string)";
} 
