use strict;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBLoader;
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
my $xmap       = $conf{'x_map'};
my $map        = $conf{'human_map'};
my $dbname     = $conf{'db'};
my $host       = $conf{'host'};

my %map;
my %ref_map;

print STDERR "Connecting to the database...\n";
my $enslocator = "Bio::EnsEMBL::DBSQL::DBAdaptor/host=$host;dbname=$dbname;user=ensadmin;perlonlyfeatures=1";
print STDERR "Bio::EnsEMBL::DBSQL::DBAdaptor/host=$host;dbname=$dbname;user=ensadmin;perlonlyfeatures=1\n";
my $db =  Bio::EnsEMBL::DBLoader->new($enslocator);

my $adaptor = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($db);

open (REFSEQ,"$refseq_gnp") || die "Can't open $refseq_gnp\n";
#Read the file by genbank entries (separated by //) 
$/ = "\/\/\n";
while (<REFSEQ>) {
#This subroutine store for each NP (refseq protein accession number) its corresponding NM (DNA accession number)
    my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
    my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;

    #print STDERR "$prot_ac\n";
    $ref_map{$prot_ac} = $dna_ac;
}
#Put back the default (new line) for reading file
$/ = "\n"; 

open (XMAP,"$xmap") || die "Can't open $xmap\n";

while (<XMAP>) {
    chomp;
    my ($targetid,$targetdb,$xac,$xdb,$xid,$xsyn) = split (/\t/,$_);
    
    #print STDERR "$targetid,$targetdb,$xac,$xdb,$xid,$xsyn\n";

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
    #print STDERR "Loading data in the database\n";
    chomp;
    my ($queryid,$tid,$tag,$queryperc,$targetperc) = split (/\t/,$_);
    #print STDERR "$queryid,$tid,$tag,$queryperc,$targetperc\n";

    
    #print STDERR "TARGETID0: $tid\n";
    if ($tid ne "orphan") {
	if ($tid =~ /^NP_\d+/) {
	    
	    ($tid) = $tid =~ /^(NP_\d+)/;
	    #print STDERR "$tid\n";
	    $tid = $ref_map{$tid};
	}
	
	#print STDERR "TARGETID1: $tid\n";
	
	if (defined $tid) {
	    
	    my @array = @{$map{$tid}};
	    
	    foreach my $a(@array) {
		#print STDERR $a->xDB."\n"; 
		
		if (($a->xDB eq "SPTREMBL") || ($a->xDB eq "SPTR") || ($a->xDB eq "REFSEQ")) {
		    #print STDERR "IDT: $queryperc\t$targetperc\n";
		    my $dbentry = Bio::EnsEMBL::IdentityXref->new
			( -adaptor => $adaptor,
			  -primary_id => $a->xAC,
			  -display_id => $a->xID,
			  -version => 1,
			  -release => 1,
			  -dbname => $a->xDB);
		    
		    $dbentry->query_identity($queryperc);
		    $dbentry->target_identity($targetperc);
				    
		    my @synonyms = split (/;/,$a->xSYN);
		    
		
		    foreach my $syn (@synonyms) {
			if ($syn =~ /\S+/) {
			    $dbentry->add_synonym($syn);
			}
		    }
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
		    $adaptor->store($dbentry,$queryid,"Translation");
		    
		}
	    }
	}
	
	else  {
	    print STDERR " not defined\n";
	}  
	
	
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


