#Contact: Emmanuel Mongin (mongin@ebi.ac.uk)

use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf; # configuration options


# global vars

my $org_list = $conf{'organism_list'};
my $refseq_gnp = $conf{'refseq_gnp'};
my $xmap       = $conf{'x_map_out'};
my $mapping    = $conf{'pmatch_out'};
my $dbname     = $conf{'db'};
my $host       = $conf{'host'};
my $user       = $conf{'dbuser'};
my $pass       = $conf{'password'};
my $port       = $conf{'port'};
my $organism   = $conf{'organism'};
my $query_pep  = $conf{'query'};
my $refseq_pred = $conf{'refseq_pred_gnp'};
my $help;


&GetOptions(
	    'help' => \$help,
	    );

if ($help) {
    print STDERR $conf{'help'}."\n";
    exit();
}

#Organism specific options
#Drosophila
my $dros_ext_annot = $conf{'dros_ext_annot'};

#Elegans
my $cefile         = $conf{'eleg_nom'};
my $type = $conf{'elegans_pseudo'};

#working option but obsolete
my $check      = $conf{'check'};

my %map;
my %cemap;
my %ref_map;
my %sp2embl;
my %ens2embl;
my %embl2sp;
my %errorflag;
my %ref_map_pred;

my $help = $conf{'help'};

print STDERR "$help\n";

#Checks

my %check;
my $seenorg = 0;

#Check if the organism is correct
foreach my $or (@{$org_list}) {
    if ($or eq $organism) {
	$seenorg = 1;
    }
}

if ($seenorg == 0) {
    print STDERR "Either the organism name you are using ($organism) is not define or is not allowed\n";
    print STDERR "Here is a list of authorised organisms:\n";
    foreach my $or (@{$org_list}) {
	print STDERR "$or\n";
    }

    exit();
}

#Organism specific checks
if($organism eq "human") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    
    foreach my $k (keys %check) {
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "mouse") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "elegans") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'elegans_nom'} = $conf{'elegans_nom'};
    $check{'elegans_pseudo'} = $conf{'elegans_pseudo'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "anopheles") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "drosophila") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};
    $check{'dros_ext_annot'} = $conf{'dros_ext_annot'};
    
    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "rat") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "zebrafish") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'pmatch_out'} = $conf{'pmatch_out'};
    $check{'db'}     = $conf{'db'};
    $check{'host'}       = $conf{'host'};
    $check{'dbuser'}       = $conf{'dbuser'};
    $check{'password'} = $conf{'password'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

#End of checks

print STDERR "Connecting to the database... $dbname:$host\n";
print STDERR "dealing with organism ".$organism."\n";

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $user,
        -dbname => $dbname,
        -host   => $host,
	-pass   => $pass,
	-port   => $port,
        -driver => 'mysql',
	);



my $adaptor = $db->get_DBEntryAdaptor();

if (($organism eq "human") || ($organism eq "mouse") || ($organism eq "rat") || ($organism eq "zebrafish")) {
    print STDERR "Reading Refseq file\n";
    open (REFSEQ,"$refseq_gnp") || die "Can't open REFSEQ $refseq_gnp\n";
#Read the file by genbank entries (separated by //) 
    $/ = "\/\/\n";
    while (<REFSEQ>) {
#This subroutine store for each NP (refseq protein accession number) its corresponding NM (DNA accession number)
      my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
      my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;
      
      $ref_map{$prot_ac} = $dna_ac;
    }
    #Put back the default (new line) for reading file
    $/ = "\n"; 
}
close(REFSEQ);

open (XMAP,"$xmap") || die "Can't open XMAP $xmap\n";

print STDERR "Reading X_map ($xmap)\n";

while (<XMAP>) {
    #print STDERR;
    chomp;
    my ($targetid,$targetdb,$xac,$xdb,$xid,$xsyn,$status) = split (/\t/,$_);
    #print STDERR "tag ".$xac."\n";
    if ($check eq "yes") {
#Get the all of the EMBL accessions for a given SP
	if (($targetdb eq "SPTR") && ($xdb eq "EMBL")) {
	    push(@{$sp2embl{$targetid}},$xac);
	}
    }

    if ($targetid =~ /^NP_\d+/) {
	
	    ($targetid) = $targetid =~ /^(NP_\d+)/;
	    $targetid = $ref_map{$targetid};
	}


    if ($xac =~ /^NP_\d+/) {
	
	    ($xac) = $xac =~ /^(NP_\d+)/;
	    $xac = $ref_map{$xac};
	}

    if ($xid =~ /^NP_\d+/) {
	
	($xid) = $xid =~ /^(NP_\d+)/;
	$xid = $ref_map{$xid};
    }


    if ($targetid =~ /^XP_\d+/) {
	
	    ($targetid) = $targetid =~ /^(XP_\d+)/;
	    $targetid = $ref_map_pred{$targetid};
	}


    if ($xac =~ /^XP_\d+/) {
	
	    ($xac) = $xac =~ /^(XP_\d+)/;
	    $xac = $ref_map_pred{$xac};
	}

    if ($xid =~ /^XP_\d+/) {
	
	($xid) = $xid =~ /^(XP_\d+)/;
	$xid = $ref_map_pred{$xid};
    }


    my $p= Desc->new;
    $p->targetDB($targetdb);
    $p->xAC($xac);
    $p->xDB($xdb);
    $p->xID($xid);
    $p->xSYN($xsyn);
    $p->stat($status);

    push(@{$map{$targetid}},$p);
}

close (XMAP);

if ($check eq "yes") {
    open (QUERY,"$query_pep") || die "Can't open QUERY PEP $query_pep\n";
    while (<QUERY>) {
	if ($_ =~ /^>\S+\s*\S+\s* Clone:\S+/) {
	    my ($pepac,$cloneac) = $_ =~ /^>(\S+)\s*\S+\s* Clone:(\S+)/; 
	    $ens2embl{$pepac} = $cloneac;
	    $embl2sp{$cloneac} = $pepac;
	}
    }
    close (QUERY);
}

open (MAP,"$mapping") || die "Can't open MAPPING $mapping\n";

print STDERR "Reading mapping output\n";

MAPPING: while (<MAP>) {
    my $target;
    my $analysis = Bio::EnsEMBL::Analysis->new(
					       -db           => 'SPTR/ANOSUB',
					       -program      => 'exonerate',
					       -program_file => 'exonerate',
					       -gff_source   => 'mapping',
					       -gff_feature  => 'mapping',
					       -module       => 'NULL',
					       -logic_name   => 'mapping'
					       );

    chomp;
#    my ($queryid,$tid,$tag,$queryperc,$targetperc) = split (/\t/,$_);
#    my ($tid,$queryid,$tag,$targetperc,$queryperc) = split (/\t/,$_);
    
     my ($queryid,$queryperc,$qalignstart,$qalignend,$tid,$targetperc,$talignstart,$talignend,$score,$cigarline) = split (/\,/,$_);

    my $m = $tid; 
    #print STDERR "$queryid,$tid,$tag,$queryperc,$targetperc\n";

    if ($tid =~ /^NP_\d+/) {
	
	($tid) = $tid =~ /^(NP_\d+)/;
	$tid = $ref_map{$tid};
    }

    if ($tid =~ /^XP_\d+/) {
	
	($tid) = $tid =~ /^(XP_\d+)/;
	$tid = $ref_map_pred{$tid};
    }
    
    if (($tid =~ /^(\w+-\d+)/)&&($organism ne "anopheles")) {
	($tid) = $tid =~ /^(\w+)-\d+/;
    }
 
    #print STDERR "TID: $tid\n";
   
    if ((defined $tid) && (defined $map{$tid})) {
	
	
	my @array = @{$map{$tid}};
	
	
	foreach my $a(@array) {
#If the target sequence is either an SPTR or RefSeq accession number, we have some information concerning the percentage of identity (that the sequences we directly used for the pmatch mapping) 
	    
	    

	    if (($a->xDB eq "SPTREMBL") || ($a->xDB eq "SWISSPROT") || ($a->xDB eq "RefSeq") || ($a->xDB eq "ANOSUB") || $a->xDB eq 'BRIGGSAE_HYBRID') {

		my $dbentry = Bio::EnsEMBL::IdentityXref->new
		    ( -adaptor => $adaptor,
		      -primary_id => $a->xAC,
		      -display_id => $a->xID,
		      -version => 1,
		      -release => 1,
		      -dbname => $a->xDB);

		$dbentry->status($a->stat);

		if (($check eq "yes") && (($a->xDB eq "SPTREMBL") || ($a->xDB eq "SWISSPROT"))) {

		    if (($sp2embl{$a->xAC}) && ($ens2embl{$queryid})) {

			if (grep($ens2embl{$queryid}=~ /$_/,@{$sp2embl{$a->xAC}})) {
		       

			}
			else {
			    foreach my $b(@{$sp2embl{$a->xAC}}) {
				if ($embl2sp{$b}) {
				    print "DODGY: ".$a->xAC."\n";
				    next MAPPING;
				}
			    }
			}
		    }
		}


		$dbentry->query_identity($queryperc);
		$dbentry->target_identity($targetperc);
		$dbentry->score($score);
		$dbentry->cigar_line($cigarline);
		$dbentry->hit_start($talignstart);
		$dbentry->hit_end($talignend);
		$dbentry->translation_start($qalignstart);
		$dbentry->translation_end($qalignend);
		$dbentry->analysis($analysis);
		    
		my @synonyms = split (/;/,$a->xSYN);
		
		
		foreach my $syn (@synonyms) {
		    if ($syn =~ /\S+/) {
			$dbentry->add_synonym($syn);
		    }
			}
		if($queryid == 0){
		  die "have no translation_id $!";
		}
		#print STDERR "storing ".$dbentry->dbname." ".$dbentry->primary_id." with ".$queryid."\n";
		$adaptor->store($dbentry,$queryid,"Translation");
	    }
	    
	    elsif ($a->xDB eq "GO") {
		print STDERR "HERE\n";

		my $dbentry = Bio::EnsEMBL::GoXref->new 
		     ( -adaptor => $adaptor,
		      -primary_id => $a->xAC,
		      -display_id => $a->xAC,
		      -version => 1,
		      -release => 1,
		      -dbname => $a->xDB);

		$dbentry->status($a->stat);
		
		$dbentry->add_linkage_type($a->xID);
		$adaptor->store($dbentry,$queryid,"Translation");
	    }

	    
	    else {
		my $dbentry = Bio::EnsEMBL::DBEntry->new
		    ( -adaptor => $adaptor,
		      -primary_id => $a->xAC,
		      -display_id => $a->xID,
		      -version => 1,
		      -release => 1,
		      -dbname => $a->xDB );
		$dbentry->status($a->stat);
		
		
		my @synonyms = split (/;/,$a->xSYN);
		
			
		foreach my $syn (@synonyms) {
		    if ($syn =~ /\S+/) {
			$dbentry->add_synonym($syn);
		    }
		}
		if($queryid == 0){
		  die "have no translation_id $!";
		}
		#print STDERR "storing ".$dbentry->dbname." ".$dbentry->primary_id." with ".$queryid."\n";
		$adaptor->store($dbentry,$queryid,"Translation");
		    
	    }
	}
    }
    
	
    else  {
	 print STDERR " $tid not defined in x_map...hum, not good\n";
    }  
}


if ($organism eq "anopheles") {
    open (ANOANNOT,"/acari/work4/mongin/anopheles_mai/mapping/Primary/celera_mapping.txt") || die "Can't open Anno file\n";
    
    while (<ANOANNOT>) {
	chomp;

	my ($stable,$id,$extdb) = split;
	
	my $query = "select translation_id from translation_stable_id where stable_id = '$stable'";
	my $sth = $db->prepare($query);
	$sth->execute();
	
	while (my $trans_id = $sth->fetchrow) {
	    if ($trans_id) {
		
		my $dbentry = Bio::EnsEMBL::DBEntry->new
		    ( -adaptor => $adaptor,
		      -primary_id => $id,
		      -display_id => $id,
		      -version => 1,
		      -release => 1,
		      -dbname => $extdb );
		$dbentry->status("XREF");
		$adaptor->store($dbentry,$trans_id,"Translation");
	    }
	}
	
    }
}

if ($organism eq "drosophila") {
    open (DROSANNOT,"$dros_ext_annot") || die "Can't open Dros file $dros_ext_annot\n";
    
    my %seen;
    while (<DROSANNOT>) {
	chomp;
	
	my ($cg,$extdb,$ac) = split;
	
	if ($extdb eq "Name") {
	    $extdb = "flybase_symbol";
	}
	
	if ((($extdb eq "flybase") || ($extdb eq "FlyBase")) && ($ac =~ /FBgn/)) {
	    $extdb = "flybase_gene";
	}

	if ($ac =~ /FBan/) {
	    next;
	}

	#Here we get CG accession numbers ($cg) corresponding to gene_stable_id in Ensembl, get all of the translation internal ids for the given entry

	my $query = "select t.translation_id from transcript t, gene_stable_id g where g.gene_id = t.gene_id and g.stable_id = '$cg'";
	my $sth = $db->prepare($query);
	$sth->execute();

	while (my $trans_id = $sth->fetchrow) {
	
	    if ($seen{$cg} != 1) {
		my $extdb = "drosophila_gene_id"; 
		my $dbentry = Bio::EnsEMBL::DBEntry->new
		    ( -adaptor => $adaptor,
		      -primary_id => $cg,
		      -display_id => $cg,
		      -version => 1,
		      -release => 1,
		      -dbname => $extdb );
		$dbentry->status("KNOWNXREF");
		if($trans_id == 0){
		    die "have no translation_id $!";
		}
		$adaptor->store($dbentry,$trans_id,"Translation");
		$seen{$trans_id} = 1;
	    }

	    if ($trans_id) {
	    #Create a new dbentry object
		my $dbentry = Bio::EnsEMBL::DBEntry->new
		    ( -adaptor => $adaptor,
		      -primary_id => $ac,
		      -display_id => $ac,
		      -version => 1,
		      -release => 1,
		      -dbname => $extdb );
		$dbentry->status("KNOWNXREF");
		if($trans_id == 0){
		  die "have no translation_id $!";
		}
		$adaptor->store($dbentry,$trans_id,"Translation");
	    }
	}

    }

}



if ($organism eq "elegans") {
  #print STDERR " parsing wormbase information\n";
  open (IN,"$cefile") || die "can't open file\n";
  #print STDERR "have opened ".$cefile."\n";
  while (<IN>) {
    #print STDERR "have open file\n";
    #print STDERR;
    chomp;
    my @array = split;
    #print STDERR "mapping ".$array[0]." to ".$array[1]."\n";
    $cemap{$array[0]} = $array[1];
  }
  close (IN);
  
  my $adaptor = $db->get_DBEntryAdaptor();
  
  my $db1 = "wormbase_gene";
  my $db2 = "wormbase_transcript";
  my $db3 = "wormpep_id";
  my $db4 = "wormbase_pseudogene";
  my $query = "select t.translation_id, ts.stable_id, gs.stable_id, g.type from transcript t, gene_stable_id gs, transcript_stable_id ts, gene g where t.gene_id = gs.gene_id and t.transcript_id = ts.transcript_id and t.gene_id = g.gene_id";
    my $sth = $db->prepare($query);
    $sth->execute();
    while (my @res = $sth->fetchrow) {
	my $transl_dbid = $res[0];
	my $transc_stable_id = $res[1];
	my $gene_stable_id = $res[2];
	my $gene_type = $res[3];
	if($gene_type ne $type){
	my $dbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $gene_stable_id,
	      -display_id => $gene_stable_id,
	      -version => 1,
	      -release => 1,
	      -dbname => $db1);
	$dbentry->status("KNOWNXREF");
	if($transl_dbid == 0){
	  die "have no translation_id  for $transc_stable_id $!";
	}
	#print STDERR "storing ".$dbentry->dbname." ".$dbentry->primary_id." with ".$transl_dbid."\n";
	$adaptor->store($dbentry,$transl_dbid,"Translation");
	
	my $transdbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $transc_stable_id,
	      -display_id => $transc_stable_id,
	      -version => 1,
	      -release => 1,
	      -dbname => $db2);
	$transdbentry->status("KNOWNXREF");
	if($transl_dbid == 0){
	  die "have no translation_id $!";
	}
	#print STDERR "storing ".$dbentry->dbname." ".$dbentry->primary_id." with ".$transl_dbid."\n";	
	$adaptor->store($transdbentry,$transl_dbid,"Translation");

	my $ce = $cemap{$transc_stable_id};
	
	if ($ce) {
	  #print STDERR "$db3\t$ce\n";
	  my $ceentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $ce,
	      -display_id => $ce,
	      -version => 1,
	      -release => 1,
	      -dbname => $db3);
	  $ceentry->status("KNOWNXREF");
	  if($transl_dbid == 0){
	    die "have no translation_id $!";
	  }
	  #print STDERR "storing ".$dbentry->dbname." ".$dbentry->primary_id." with ".$transl_dbid."\n";
	  $adaptor->store($ceentry,$transl_dbid,"Translation");
	}
      }else{
	my $transdbentry = Bio::EnsEMBL::DBEntry->new
	    ( -adaptor => $adaptor,
	      -primary_id => $transc_stable_id,
	      -display_id => $transc_stable_id,
	      -version => 1,
	      -release => 1,
	      -dbname => $db4);
	$transdbentry->status("PSEUDO");
	$adaptor->store($transdbentry,$transl_dbid,"Transcript");
      }
    }
  }

sub usage {
    
    print STDERR <<HELP

Usage: maps2db.pl 
One of the element of the configuration file has not been properly loaded
for the organism $organism
Please fill in properly your configuration file

Here is your set up:
HELP
;

 foreach my $k (keys %check) {
	print STDERR "$k:\t$check{$k}\n";
    }



  exit();
}

###############
#Some OO stuff#
###############

package Desc;

=head2 new

 Title   : new
 Usage   : $obj->new($newval)
 Function: 
 Returns : value of new
 Args    : newvalue (optional)


=cut

sub new{
 my $class= shift;
    my $self = {};
    bless $self,$class;
    return $self;
}

=head2 targetDB

 Title   : targetDB
 Usage   : $obj->targetDB($newval)
 Function: 
 Returns : value of targetDB
 Args    : newvalue (optional)


=cut

sub targetDB{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'targetDB'} = $value;
    }
    return $obj->{'targetDB'};

}


=head2 xAC

 Title   : xAC
 Usage   : $obj->xAC($newval)
 Function: 
 Returns : value of xAC
 Args    : newvalue (optional)


=cut

sub xAC{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'xAC'} = $value;
    }
    return $obj->{'xAC'};

}

=head2 xDB

 Title   : xDB
 Usage   : $obj->xDB($newval)
 Function: 
 Returns : value of xDB
 Args    : newvalue (optional)


=cut

sub xDB{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'xDB'} = $value;
    }
    return $obj->{'xDB'};

}

=head2 xID

 Title   : xID
 Usage   : $obj->xID($newval)
 Function: 
 Returns : value of xID
 Args    : newvalue (optional)


=cut

sub xID{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'xID'} = $value;
    }
    return $obj->{'xID'};

}

=head2 xSYN

 Title   : xSYN
 Usage   : $obj->xSYN($newval)
 Function: 
 Returns : value of xSYN
 Args    : newvalue (optional)


=cut

sub xSYN{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'xSYN'} = $value;
    }
    return $obj->{'xSYN'};

}

=head2 stat

 Title   : stat
 Usage   : $obj->stat($newval)
 Function: 
 Returns : value of stat
 Args    : newvalue (optional)


=cut

sub stat{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'stat'} = $value;
    }
    return $obj->{'stat'};

}


