#Contact: Emmanuel Mongin (mongin@ebi.ac.uk)

use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::SeqIO;

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf; # configuration options

# global vars

my $refseq_gnp = $conf{'refseq_gnp'};
my $xmap       = $conf{'x_map_out'};
my $map        = $conf{'pmatch_out'};
my $dbname     = $conf{'db'};
my $host       = $conf{'host'};
my $user       = $conf{'dbuser'};
my $pass       = $conf{'password'};
my $organism   = $conf{'organism'};

my %map;
my %ref_map;

if ((!defined $organism) || (!defined $xmap) || (!defined $map) || (!defined $dbname) || (!defined $host)) {
    die "\nSome basic options have not been set up, have a look at mapping_conf\nCurrent set up (required options):\norganism: $organism\nx_map: $xmap\npmatch_out: $map\ndb: $dbname\nhost: $host\n\n";
}

print STDERR "Connecting to the database...\n";


my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => $user,
        -dbname => $dbname,
        -host   => $host,
	-pass   => $pass,			     
        -driver => 'mysql',
	);

my $adaptor = $db->get_DBEntryAdaptor();

if ($organism eq "human") {
    open (REFSEQ,"$refseq_gnp") || die "Can't open $refseq_gnp\n";
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


open (XMAP,"$xmap") || die "Can't open $xmap\n";

while (<XMAP>) {
    chomp;
    my ($targetid,$targetdb,$xac,$xdb,$xid,$xsyn) = split (/\t/,$_);

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

    my $p= Desc->new;
    $p->targetDB($targetdb);
    $p->xAC($xac);
    $p->xDB($xdb);
    $p->xID($xid);
    $p->xSYN($xsyn);

    push(@{$map{$targetid}},$p);
}

close (XMAP);



open (MAP,"$map") || die "Can't open $map\n";

while (<MAP>) {
    my $target;
    chomp;
    my ($queryid,$tid,$tag,$queryperc,$targetperc) = split (/\t/,$_);
    
    my $m = $tid; 

    if ($tid =~ /^NP_\d+/) {
	
	($tid) = $tid =~ /^(NP_\d+)/;
	$tid = $ref_map{$tid};
    }
    
    if ($tid =~ /^(\w+-\d+)/) {
	($tid) = $tid =~ /^(\w+)-\d+/;
    }
    

    if ((defined $tid) && (defined $map{$tid})) {
	
	
	my @array = @{$map{$tid}};
	
	foreach my $a(@array) {
#If the target sequence is either an SPTR or RefSeq accession number, we have some information concerning the percentage of identity (that the sequences we directly used for the pmatch mapping) 
	    if (($a->xDB eq "SPTREMBL") || ($a->xDB eq "SWISS-PROT") || ($a->xDB eq "RefSeq")) {
		my $dbentry = Bio::EnsEMBL::IdentityXref->new
		    ( -adaptor => $adaptor,
		      -primary_id => $a->xAC,
		      -display_id => $a->xID,
		      -version => 1,
		      -release => 1,
		      -dbname => $a->xDB);

		print STDERR $a->xAC,"\t",$a->xID,"\t",$a->xDB,"\n";

		$dbentry->query_identity($queryperc);
		$dbentry->target_identity($targetperc);
		
		my @synonyms = split (/;/,$a->xSYN);
		
		
		foreach my $syn (@synonyms) {
		    if ($syn =~ /\S+/) {
			$dbentry->add_synonym($syn);
		    }
			}

		print STDERR "Calling store on $queryid\n";
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
		
		
		
		my @synonyms = split (/;/,$a->xSYN);
		
			
		foreach my $syn (@synonyms) {
		    if ($syn =~ /\S+/) {
			$dbentry->add_synonym($syn);
		    }
		}
		print STDERR "Calling store1 on $queryid\n";
		$adaptor->store($dbentry,$queryid,"Translation");
		    
	    }
	}
    }
    
	
    else  {
	#print STDERR " $tid not defined in x_map...hum, not good\n";
    }  
}


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


