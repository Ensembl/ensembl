#!/usr/local/bin/perl -w
# $Header$
# 
# Script to create a tab-separated file (to be loaded into the 'gene_description' table)
# from protein database file(s) either SWISSPROT/SPTREMBL or RefSeq and
# the protein mapping information stored in the currently used Ensembl core database.
#
# When file is created and checked for "junk" data, the same script is used to load the
# corresponding data to the Ensembl core database in 'gene_description' table

# EXIT STATUS
# 0 all is fine
# 1 Problem in parsing protein database file
# 2 Problem when loading data
# 3 Gene_description_file does not exist
# 4 Problem with an Acc without description

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

$| = 1;

# db_name used in the Ensembl core database (to be checked in table externalDB) and
# used in the protein mapping.
my %protein_db_naming = ('swissprot' => "SWISSPROT",
			 'sptrembl' => "SPTREMBL",
			 'refseq' => "RefSeq");

my @databases_used_for_guessing_descriptions;

## list of regexps tested, in order of increasing preference, when
## deciding which description to use (used by sub compare_desc):

my @word_order = qw(unknown hypothetical putative novel probable [0-9]{3} kDa fragment cdna protein);

## Getting parameters/options and usage

my $usage = "
 Usage: $0 -host host \\
           -port port \\
           -dbuser dbuser \\
           -debug \\
           -dbname dbname 
           -consortium dbname in external_db table (e.g. ZFIN_ID for zebrafish)
           -r regexps_file
           fileprotein_database_file(s) > gene-descriptions.tab

 where protein_database_file could SWISSPROT/SPTREMBL entries in SWISSPROT-like format
                                   and/or
                                   RefSeq entries in gnp-like format.
                                   and/or
                                   \"consortium\" description file, which format should
                                   match this regexp /^(\\S+)\\t(.*)\$/, making sure that
                                   the mapping includes populating identity_xref table.

 OR to load the data from gene-descriptions.tab file to \'gene-description\' table

        $0 -host host \\
           -port port \\
           -dbuser dbuser \\
           -dbname dbname \\
           -dbpass password \\
           -load gene-descriptions.tab
\n";

my $help = 0;
my ($host, $port, $dbname, $dbuser, $dbpass, $gene_description_file);
my $regexp_file;
my $debug = 0;
my $consortium = "";

unless (GetOptions('help' => \$help,
		   'host=s' => \$host,
		   'port=i' => \$port,
		   'dbname=s' => \$dbname,
		   'dbuser=s' => \$dbuser,
		   'dbpass=s' => \$dbpass,
		   'r=s' => \$regexp_file,
		   'consortium' => \$consortium,
		   'load=s' => \$gene_description_file,
		   'debug' => \$debug)) {
  die $usage;
}

if ($help) {
  print $usage;
  exit 0;
}

$protein_db_naming{$consortium} = $consortium;

## Some checks before starting

if (! defined $host) {
  warn "\n Must specify a host with -h\n";
  die $usage;
} elsif (! defined $dbname) {
  warn "\n Must specify a database with -d\n";
  die $usage;
} elsif (! defined $dbuser) {
  warn "\n Must specify a database username with -u\n";
  die $usage;
}

if (defined $gene_description_file || defined $dbpass) {

  if (! defined $dbpass) {
    warn "\n Must specify a database password with -p to load your data\n";
    die $usage;
  } elsif (! defined $gene_description_file) {
    warn "
 Assuming you want to load your data in the database...
 If so you must specify a path for your gene-descriptions.tab file with -load.
 If not -p is not needed\n";
    die $usage;
  }

} elsif (scalar @ARGV ==0) {
  warn "\n Must specify an protein database input file\n";
  die $usage;
}

## Loading data to database if requested and exit

if (defined $gene_description_file) {
  if (load_data($gene_description_file,$host,$port,$dbname,$dbuser,$dbpass) == 1) {
    exit 0;
  } else {
    warn "Problem when loading data
exit 2";
    exit 2;
  }
}

## Parsing of protein database flat files(Swissprot/SPTrEMBL/RefSeq) used in protein mapping

my %sp_desc;

while (defined (my $protein_database_file = shift @ARGV)) {
  unless (parse_protein_database($protein_database_file,\%protein_db_naming,\%sp_desc) == 1) {
    warn "Problem in parsing protein database file : $protein_database_file
exit 1";
    exit 1;
  }
}

## Connecting to the Ensembl database where protein mapping is stored

print STDERR "Connecting to $dbname $host......";

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor (-host => $host,
                                             -port => $port,
					     -dbname => $dbname,
					     -user => $dbuser);
print STDERR "Done\n";

foreach my $db (@databases_used_for_guessing_descriptions) {
  $db = "'".$db."'";
}
my $db_query = join ",",@databases_used_for_guessing_descriptions;

print STDERR "Preparing query...";


my $sth = $db->prepare("
SELECT
 tsl.translation_id,tsc.gene_id,xdb.db_name,x.dbprimary_acc,ix.query_identity,ix.target_identity
FROM
 translation tsl,transcript tsc,object_xref ox,xref x,external_db xdb,identity_xref ix
WHERE
 tsl.transcript_id = tsc.transcript_id AND
 tsl.translation_id = ox.ensembl_id AND
 ox.xref_id = x.xref_id AND
 x.external_db_id = xdb.external_db_id AND
 xdb.db_name in ($db_query) AND
 ox.object_xref_id = ix.object_xref_id
ORDER BY
 tsc.gene_id asc,
 xdb.db_name desc,
 x.dbprimary_acc asc");

print STDERR "Done\n";
print STDERR "Executing query...";

unless ($sth->execute()) {
  $db->throw("Failed execution of a select query");
}
print STDERR "Done\n";
## Assigning the descriptions to Ensembl genes.

my %gene_desc;

while (my ($ensp, $ensg, $db, $acc, $qy_percid, $tg_percid) = $sth->fetchrow_array()) {

  my $desc;
  if (defined $gene_desc{$ensg}) {
    my ($prevdb, $prev_desc, $prev_acc, $prev_qy_percid, $prev_tg_percid)  = @{$gene_desc{$ensg}}; 

    if ($prevdb eq $db &&
	($db eq $protein_db_naming{$consortium} ||
	 $db eq $protein_db_naming{'swissprot'} ||
	 $db eq $protein_db_naming{'refseq'})) {
      if ($prev_qy_percid < $qy_percid) {
	$desc = $sp_desc{"$db:$acc"};
	$gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # taking better SWISSPROT/RefSeq/Consortium desc.
	next;
      } elsif ($prev_qy_percid == $qy_percid &&
	       $prev_tg_percid < $tg_percid) {
	$desc = $sp_desc{"$db:$acc"};
	$gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # taking better SWISSPROT/RefSeq/Consortium desc.
	next;
      } else {
	next;                 # nothing to change
      } 
    }

    if ($db eq $protein_db_naming{$consortium}) { 
      $desc = $sp_desc{"$db:$acc"};
      $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # kick out the SWISSPROT/RefSeq/SPTREMBL desc -> priority to consortium desc
      next;
    }

    if ($db eq $protein_db_naming{'swissprot'} &&
	$prevdb ne $protein_db_naming{$consortium}) { 
      $desc = $sp_desc{"$db:$acc"};
      $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # kick out the RefSeq/SPTREMBL desc -> priority to SWISSPROT desc
      next;
    }

    if ($db eq $protein_db_naming{'refseq'} &&
	$prevdb ne $protein_db_naming{'swissprot'} &&
	$prevdb ne $protein_db_naming{$consortium}) { 
      $desc = $sp_desc{"$db:$acc"};
      $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # kick out the SPTREMBL desc -> priority to RefSeq desc
      next;
    }
    
    if ($prevdb eq $db &&
	$db eq $protein_db_naming{'sptrembl'}) {   
      $desc = $sp_desc{"$db:$acc"}; 
      if (compare_desc($prev_desc,$desc) < 0) {
	# new desc is better
	# warn "new better: $desc (old was: $prev_desc)\n";
	$gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ];
	next;
      } else {
	# warn "old better: $prev_desc (new is: $desc)\n";
	next;
      }
#      die "should not reach this point: only know SWISSPROT and SPTREMBL";
    }
  } else {
    $desc = $sp_desc{"$db:$acc"};
    $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ];
  }
}

## Load regexps

open (REGEXP, $regexp_file) || die "Can't open REGEX file: $regexp_file\n";
my @regexps;

while (my $line = <REGEXP>) {
  next if ($line =~ /^\#.*/ || $line =~ /^$/);
  chomp $line;
  push @regexps, $line;
}

close REGEXP;

##  Dump Ensembl gene descriptions assignations to stdout

foreach my $ensg ( keys %gene_desc )  { 

  my ($db,$desc,$acc,$qy_percid,$tg_percid) = @{$gene_desc{$ensg}};

  ### final cleanup 
  ### get rid of the useless descriptions matching regexps in the regexp_file provided

  if (defined $desc) {
    $_ = $desc;

    foreach my $regexp (@regexps) {
      s/$regexp//ig;
    }
    
    unless ($_ eq "") {
      if ($debug) {
	print STDOUT "$ensg\t".$desc." [Source:$db;Acc:$acc]\t>>$_<<\n"; 
      } else {
	print STDOUT "$ensg\t".$desc." [Source:$db;Acc:$acc]\n"; 
      }
    } else {
      warn "throwing away: $desc\n";
    }
  } else {
    warn "Description for Acc $acc in $db file not defined.
Check in Ensembl core database if all Acc included in protein mapping are in the protein database files
exit 4";
    exit 4;
  }
}

## following taken from ensembl-external/scripts/family-input.pl

sub compare_desc {
  my ($a, $b) = @_; 

  my ($am, $bm);
  foreach my $w (@word_order) {
    $am = ($a =~ /$w/i)?1:0;
    $bm = ($b =~ /$w/i)?1:0;
    
    if ($am  != $bm ) {             # ie, one matches, other doesn't
      if ( $am == 1 ) {           # first one worse than second 
	return -1;
      } else { 
	return 1; 
      }
    }
  }
  # still look same; base result on length: longer is better
  return length($a) <=> length($b);
} # compare_desc

sub load_data {
  my ($gene_description_file,$host,$port,$dbname,$dbuser,$dbpass) = @_;

  if (-e $gene_description_file) {
    symlink $gene_description_file,"$ENV{PWD}/gene_description_file.symlink";
  } else {
    warn "File $gene_description_file does not exist
exit 3";
    exit 3;
  }
 
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor (-host => $host,
                                               -port => $port,
					       -dbname => $dbname,
					       -user => $dbuser,
					       -pass => $dbpass);

  my $sth = $db->prepare("load data infile '$ENV{PWD}/gene_description_file.symlink' into table gene_description");

  unless ($sth->execute()) {
    unlink "$ENV{PWD}/gene_description_file.symlink";
    $db->throw("Failed execution of loading data");
  } else {
    unlink "$ENV{PWD}/gene_description_file.symlink";
    print STDERR "Data in $gene_description_file file successfully loaded in \'gene_description\' table in ",$dbname,"\@",$host,"\n";
    return 1;
  }
} # load_data

sub parse_protein_database {
  my ($protein_database_file,$protein_db_naming_href,$sp_desc_href) = @_;

  my $db;
  my $acc; 
  my $desc = "";

  open PROTDBF,$protein_database_file || die "$protein_database_file:$!";
  
  while (defined (my $line = <PROTDBF>)) {
    next if ($line =~ /^$/);
    if ($line =~ /^ID\s{3}.*$/o) {
      if ($line =~ /^ID\s{3}.*\s+PRELIMINARY;\s+.*$/o) {
	$db = $protein_db_naming_href->{'sptrembl'};
	push @databases_used_for_guessing_descriptions, $db unless (grep /^$db$/, @databases_used_for_guessing_descriptions);
      } elsif ($line =~ /^ID\s{3}.*\s+STANDARD;\s+.*$/o) { 
	$db = $protein_db_naming_href->{'swissprot'};
	push @databases_used_for_guessing_descriptions, $db unless (grep /^$db$/, @databases_used_for_guessing_descriptions);
      } elsif ($line =~ /^(\S+)\t(.*)$/) {
	$db = $consortium;
	($acc,$desc) = ($1,$2);
	push @databases_used_for_guessing_descriptions, $db unless (grep /^$db$/, @databases_used_for_guessing_descriptions);
	$sp_desc{"$db:$acc"} = $desc; 
	$db = undef;
	$acc = undef; 
	$desc = "";
	
      } else {
	chomp $line;
#took away the warning as many lines in SP have ID  AC and nothing else

#	warn "\nCould not recognize line : \"$line\"\nCheck the input file format in $protein_database_file.\n\n";
#	return 0;
      }
    }
    
    if ($line =~ /^LOCUS\s+\S+\s+.*$/o) {
      $db = $protein_db_naming_href->{'refseq'};
      push @databases_used_for_guessing_descriptions, $db  unless (grep /^$db$/, @databases_used_for_guessing_descriptions);
    }
    if ($db eq $protein_db_naming_href->{'swissprot'} ||
	$db eq $protein_db_naming_href->{'sptrembl'}) {
      if ($line =~ /^AC\s{3}(\S+);.*$/o &&
	  ! defined $acc)  {
	$acc = $1;
      }
      
      if ($line =~ /^DE\s{3}(\S.*\S)$/o) {
	$desc .= " " unless ($desc eq "");
        $desc .= $1;
      }
      
      if ($line =~ /^\/\/$/o) {
	$desc =~ s/\{.*\}//g; 
	$sp_desc{"$db:$acc"} = $desc; 
	$db = undef;
	$acc = undef; 
	$desc = ""; 
      }
      
    } elsif ($db eq $protein_db_naming_href->{'refseq'}) {
      
      if ($line =~ /^DEFINITION\s+(\S.*\S)$/o) {
        $desc .= $1;
	while (defined ($line = <PROTDBF>)) {
	  last if ($line =~ /^ACCESSION\s+(\S+)\s*.*$/o);
	  if ($line =~ /^\s+(\S.*\S)$/o) {
	    $desc .= " ". $1;
	  }
	}
      }
      
      if ($line =~ /^DBSOURCE\s+(\S+)\s*.*$/o)  {
# Accession Number got in this line here during protein mapping
# not in  ACCESSION line.
	if ($line =~ /^DBSOURCE\s+REFSEQ:\s+accession\s+(\S+)\.\d+$/o)  {
	  $acc = $1;
	} else {
	  chomp $line;
	  warn "Problem in parsing DBSOURCE line: $line
Expecting to match this regexp /^DBSOURCE\\s+REFSEQ:\\s+accession\\s+(\\S+)\\.\\d+\$/o\n";
	  return 0;
	}
      }
      
      if ($line =~ /^\/\/$/o) {
	$desc =~ s/\s*\[.*\]//g; 
	$sp_desc{"$db:$acc"} = $desc; 
	$db = undef;
	$acc = undef; 
	$desc = ""; 
      }
    }

  }

  close PROTDBF;
  
  return 1;

} # parse_protein_database

exit 0;
