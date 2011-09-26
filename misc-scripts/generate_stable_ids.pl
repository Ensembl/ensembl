# Generate stable IDs for genes/transcripts/translations/exons that have none
# Start from current max stable ID + 1

use strict;
use DBI;
use Getopt::Long;
my $port = 3306 ; 
my ($host, $dbname, $user, $pass, @types, $start, $verbose );

GetOptions('dbuser|user=s'       => \$user,
	   'dbpass|pass=s'       => \$pass,
	   'dbhost|host=s'       => \$host,
	   'dbport|port=i'       => \$port,
	   'dbname=s'     => \$dbname,
	   'types=s'      => \@types,
	   'start=s'      => \$start,   # USE ENS000001 or ENS for human, ENSMUS00001 or ENSMUS for mouse etc 
                                        # don't add G/T/E/P for specific types !!!
	   'help'         => sub { usage(); exit(0); },
           'verbose!'    => \$verbose ,
            );

@types = ('gene','transcript','translation','exon') if (!@types);
@types = split(/,/,join(',',@types));



if (!$user || !$host || !$dbname || !@types) {
  usage();
  exit(1);
}

my $dbi = DBI->connect( "DBI:mysql:host=$host:port=$port;database=$dbname", $user, $pass,
			{'RaiseError' => 1}) || die "Can't connect to database\n";



foreach my $type (@types) {

  my $table = $type . "_stable_id";
  my $sth;

  # get starting stable ID, either specified or current max
  my $new_stable_id;

  if ($start) {
    $new_stable_id = $start;
  } else {
    $new_stable_id = get_highest_stable_id($dbi, $type ) ; 
  } 
  
  print "Highest, pruned $type\_stable_id found : $new_stable_id \n"  if $verbose ; 

  # get timestamp so all new stable IDs have the same created/modified dates
  $sth = $dbi->prepare("SELECT NOW()");
  $sth->execute();
  my $ts;
  if (my @row = $sth->fetchrow_array()) {
    $ts= $row[0];
  } else {
    die "Can't get timestamp\n";
  }

  # get a list of objects that don't currently have stable IDs assigned
  # and assign new ones, incrementing & preserving formatting as we go
  my $sql = "SELECT $type.${type}_id FROM $type LEFT JOIN $table sid ON $type.${type}_id=sid.${type}_id WHERE sid.stable_id IS NULL";
  $sth = $dbi->prepare($sql);
  $sth->execute();

  while (my @row = $sth->fetchrow_array()) { 
    ($new_stable_id,my $nis) = @{increment_stable_id($new_stable_id,$type)};  
    print "INSERT INTO $table VALUES($row[0],\'$nis\',1,\'$ts\',\'$ts\');\n";
  }

}


# --------------------------------------------------------------------------------

sub increment_stable_id {

  my $stable_id = shift;
  my $type = shift ;   
  
  my ($prefix, $suffix) ;  

  # check stable_id format ...  
  if ( $stable_id =~m/([a-zA-Z]+)([0-9]+)/ ){
    ($prefix,$suffix) = $stable_id =~ /([a-zA-Z]+)([0-9]+)/;   
  } elsif( $stable_id =~m/([a-zA-Z]+)/){
    $prefix = $stable_id ; 
  }else { 
    die "unrecongnized stable_id format - should match ([a-zA-Z]+)([0-9]+) or ([a-zA-Z]+) !!\n"; 
  } 
  my $new_sid; 

  if ($type=~m/gene/){ 
     $new_sid=$prefix."G"; 
  } elsif ($type=~m/transcript/){ 
     $new_sid=$prefix."T";
  } elsif ($type=~m/translation/){ 
     $new_sid=$prefix."P";
  } elsif ($type=~m/exon/){ 
     $new_sid=$prefix."E";
  }
  my $new_stable_id = sprintf "%s%011d", $new_sid , $suffix+1 ;   


  my $old = sprintf "%s%011d", $prefix, $suffix+1 ;    
  #return [$old, $new_stable_id] ; 
  my @tmp;  
  
  push @tmp, $old ; 
  push @tmp, $new_stable_id ;  
 
  return \@tmp ; 

}

# -------------------------------------------------------------------------------- 
sub get_max_stable_id_from_gene_archive { 
  my ($dbi, $type) = @_;   

  # try to get from relevant archive
  my $sth = $dbi->prepare("SELECT MAX($type) FROM gene_archive");
  $sth->execute();  

  my $rs ;  
  if (my @row = $sth->fetchrow_array ) { 
    $rs = $row[0];
  }  
  if (length($rs)> 0 ) { 
   return $rs ; 
  } else { 
   print STDERR "no entry for $type found in gene_archive table - returning undef\n" ;  
   return undef ; 
  }  

} 
# -------------------------------------------------------------------------------- 

sub get_highest_stable_id {
  my ($dbi, $type) = @_;

  my $sid = $type . "_stable_id";

  my ($highest_from_current, $highest_from_archive);

  # get highest stable ID from the relevant table 
  
  my $sth = $dbi->prepare("SELECT MAX(stable_id) FROM $sid");
  $sth->execute(); 

  if (my @row = $sth->fetchrow_array()) {
    $highest_from_current = $row[0];
  } else {
    die "Can't get max $type stable ID from $sid\n";
  }

  if (length($highest_from_current) == 0 ) {   
    print STDERR " Warning ! length of stable_id for $type is zero \n" ;  
  } 
  
  if ($type eq "exon"){   
     # Archive doesn't store information about exon_stable_ids so try without archive first ... 
     
     if ( length($highest_from_current) == 0 ) { 
       print ("\nWARNING:\n No exon_stable_id for exon found ( exon_stable_id_table empty)\n". 
           " I got no prefix to generate new stable_ids for type $type!!! - i try to use gene_archive now\n") ;   

       my $max =  get_max_stable_id_from_gene_archive($dbi, "gene_stable_id") ;    

       my $prefix ;   

       if ( length($max) >  0 ){  
         ($prefix, my $suffix) = $max =~ /([a-zA-Z]+)([0-9]+)/;  
         $prefix=~s/G$//g;  
       }else {  
           die ("ERROR : No entries in TABLE  exon_stable_id and TABLE gene_archive found\n" .
                " Don't know which species prefix to use for species \n") ; 
  
         $highest_from_current = sprintf "%s%011d", $prefix, 0; 
       }
     }
      return $highest_from_current  ; 
  }

  # and from relevant archive 
  
  $highest_from_archive =  get_max_stable_id_from_gene_archive($dbi, $type."_stable_id") ; 
  
  my $max = 
   ($highest_from_current ge $highest_from_archive) ? $highest_from_current : $highest_from_archive; 

  if (length($max)==0){
     die ("ERROR : no stable_id in TABLE gene_archive or found in $type\_stable_id - tables\n") ; 
  }

  # assuming that this is a correctly formatted stable id -> remove the G / T / P / E for exon etc. 
 
  my ($prefix,$suffix) = $max =~ /([a-zA-Z]+)([0-9]+)/;  
  if ($type =~m/exon/){  
    $prefix=~s/E$//; 
  }elsif ($type =~m/gene/){ 
    $prefix=~s/G$//; 
  }elsif ($type =~m/transcript/){ 
    $prefix=~s/T$//; 
  }elsif ($type =~m/translation/){ 
    $prefix=~s/P$//; 
  }
  return $prefix. $suffix; 
}


sub usage {

  print << "EOF";

  USAGE :  

  generate_stable_ids.pl -dbuser|user {user} 
                         -dbpass|pass {password} 
                         -dbhost|host {host}
                         -dbport|port {port} 
                         -dbname {database} 
                         -types {gene,exon,transcript,translation} 
                         -start {first stable ID}

  Argument to -types is a comma-separated list of types of stable IDs to be produced.

  If the -types argument is ommitted, stable IDs are generated for all types (gene,transcript,translation,exon).

  Assigns stable IDs to objects that currently have none. Starting stable ID is found by incrementing the highest 
  current stable ID for that type *or* by using -start argument. The stable_ids are written to STDOUT. If no -start
  option is used the script tries to find the max. given stable_id for each object by looking up the <OBJ>_stable_id
  tables in the database and the gene_archive table ( only  for gene,translation and transcript, not for exon-stable-ids!)


  Note : 

  -start option requires to not submit an initial stable_id without any Gene/Transcript/Exon/Translation ending, like 
   ENSMUS000001 ( not ENSMUSG0001 than you end up with stable-ids like ENSMUSGG001 ENSMUSGT0001...) !

  Again,
  the parameter to -start should be the stable ID you wish to start from without the Gene/Transcript/Translation/Exon identifier

  Examples :  

    - to generate only Exon-stable-ids starting with 223 for Mouse, use 

                    -start ENSMUS222 -types exon 


    - to generate Exon-and Gene-stable-ids starting with 223 for Mouse, use 

                    -start ENSMUS222 -types exon,gene

    - to generate exon,transcript,translation and gene-stable-ids for Human starting, which all start with ID 666 use 

                    -start ENS665     

    - to generate a whole new set of stable_ids ( exon,transcript, translation, gene ) starting with 1 for an organism with 
      prefix ENSNEW you can use one of the following options : 
             
                    -start ENSNEW0          <or> 
                    -start ENSNEW0          <or> 
                    -start ENSNEW00000000
             


  Produces SQL which can be run against the target database.

EOF

}

# --------------------------------------------------------------------------------
