use strict;

=head1 dump_textview

=head2 Description

The aim of this script is to produce a flat file containing text information about gathered from different databases (OMIM, Interpro, Swissprot). For OMIM and SP the link from the Ensembl database to the external one is taken from GeneDBLink. For Interpro information is taken from the interpro table.

=head2 Synopsis

get_textview_dump.pl -db ens075 -host ecs1c.sanger.ac.uk -dbuser ensro

=head2 Format of the file dumped

ENSG....|ExternalDB|ExternalAC|Text

=head2 Options

    -db        dbname (dbname= in locator)
    -host      host name (gets put as host= in locator)
    -dbuser    what username to connect to the database (dbuser= in locator)
    -interpro  name of the file which contains the interpro data (XML file)
    -omim      name of the file for OMIM
    -sp        name of the file for SP
    -output    name of the file where the output should be written 

=cut

use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $dbpass = undef;
my $dbuser = 'ensro';
my $ensdbname = 'ens075';
my $host = 'ecs1c.sanger.ac.uk';
my $output;

&GetOptions(
	    'db:s' => \$ensdbname,
	    'host:s'=> \$host,
	    'dbuser:s'=> \$dbuser,
	    'interpro:s'=>\$interpro,
	    'omim:s'=>\$omim,
	    'sp:s'=>\$sp,
	    'output:s' => \$output
	    );


my $enslocator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;dbname=$ensdbname;user=$dbuser;pass=$dbpass;perlonlyfeatures=1";
my $ensdb =  Bio::EnsEMBL::DBLoader->new($enslocator);

my %omim;
my %sp;

open (OUT, ">$output");

####################
#READ INTERPRO FILE#
####################


print STDERR "Getting Interpro\n";

#open (INTERPRO, "interpro_07_02.xml");
open (INTERPRO, "$interpro") || die "Can't open file $interpro\n";

#This file is in XML format

my @interpros;
my $interpro;
my %hash;

#Get the Interpro ID
 if (/\<interpro id\=\"(\S+)\"/) {
     $interpro->{id}=$1;
 }

#For each entry get the abstract,publist and memberlist.
#Each fiels beginns by <fieldname> and ends by <\fieldname>

ENTRY: while (<INTERPRO>) {

    my $abstract;

    if (/\<abstract/) {
	my $string=$_;
	while (<INTERPRO>) {

	    $string .= $_;

#If end of the field, get red out all of the junk and last
	    if (/\<\/abstract/) {
		
		$string =~ s/\<abstract\>//gm;
		$string =~ s/\<\/abstract\>//gm;
		$string =~ s/\<p\>//gm;
		$string =~ s/\<\/p\>//gm;
		$string =~ s/\>//gm;
		$string =~ s/\<//gm;
		$string =~ s/\///gm;
		$string =~ s/\n/ /gm;
		$string =~ s/\s+/ /gm;
		last;
	    }
	}
#Each field is referenced by using the name of the field
	$interpro->{abstract}=$string;
    }

    if (/\<publist/) {
	my $string = $_;
	while (<INTERPRO>) {
	    $string .= $_;
	    if (/\<\/publist/) {
		$string =~ s/\<publist\>//gm;
		$string =~ s/\<\/publist\>//gm;
		$string =~ s/\<authorlist\>//gm;
		$string =~ s/\<\/authorlist\>//gm;
		$string =~ s/\<title\>//gm;
		$string =~ s/\<\/title\>//gm;
		$string =~ s/\<\/publication\>//gm;
		$string =~ s/\<journal\>//gm;
		$string =~ s/\<year\>//gm;
		$string =~ s/\<\/year\>//gm;
		$string =~ s/\>//gm;
		$string =~ s/\<//gm;
		$string =~ s/\///gm;
		$string =~ s/\n/ /gm;
		$string =~ s/\s+/ /gm;
		
		last;
	    }
	}
	$interpro->{pub}=$string;
    }


    if (/\<memberlist/) {
	my $string;
	while (<INTERPRO>) {
	    my ($dbkey,$name) = $_ =~ /dbkey\=\"(\w+)\"\s*name\=\"(\w+)/;
	    
	    if ((defined $dbkey) && (defined $name)) {
		
		$string .= "$dbkey $name ";
		
	    }
	    if (/\<\/memberlist/) {
		
		last;
	    }
	}
	$interpro->{member} = $string;
    }

    
#End of the entry    
    if (/\<interpro id\=\"(\S+)\"/) {
	my $inter_id = $interpro->{id};
	my $total =  $interpro->{abstract}." ". $interpro->{pub}." ". $interpro->{member},"\n";
	$hash{$inter_id} = $total;
	
	$interpro->{id}=$1;
		
	next ENTRY;
    }
}

################
#READ OMIM FILE#
################

print STDERR "Getting OMIM\n";

#open (OMIM, "omim.txt") || die "Can't open\n";
open (OMIM, "$omim") || die "Can't open file $omim\n";

#Change the default of $_, read OMIM file by entry (*RECORD* separator)
$/ = "\*RECORD\*";

while (<OMIM>) { 
#Get OMIM ac
    my ($omim_ac) = $_ =~ /\*FIELD\* NO\n(\d+)/;

#Get all of the text corresponding to this AC    
    my $omim_tot = $_;
    
    $omim_tot =~ s/\n/ /gm;
   
    if (!defined $omim{$omim_ac}) {
	$omim{$omim_ac} = [];
    }
    if ($omim_ac) {
	$omim{$omim_ac} = $omim_tot;
    }
}

##############
#READ SP FILE#
##############


print STDERR "Getting SP\n";

#Put back the default for $_
$/ = "\n";

#Use bioperl to read SP
#my $in  = Bio::SeqIO->new(-file => 'sp_sptrembl_10_00.txt', '-format' =>'swiss');
my $in  = Bio::SeqIO->new(-file => '$sp', '-format' =>'swiss');

while ( my $seq = $in->next_seq() ) {
    my @refs;
    my @comments;

#Get AC, RC, DE, KW and comment lines
    my $ac =  $seq->accession;

    my $ke = $seq->keywords;
    
    my $desc = $seq->desc;

	
    foreach my $ref ( $seq->annotation->each_Reference ) {
	push (@refs,$ref->comment());
    }

    foreach my $comment ( $seq->annotation->each_Comment ) {
	push (@comments,$comment->text)
	}

    my $tot_sp = "$ke|$desc|@refs|@comments";
 
    if (!defined $sp{$ac}) {
	$sp{$ac} = [];
    }

    $sp{$ac} = $tot_sp;

}

#Dump files
print STDERR "Dumping file\n";

#Get all of the distinct gene ac present in genedblink
my $sth2 = $ensdb->prepare ("select distinct (gene_id) from genedblink");

$sth2->execute;

while (my @row2 = $sth2->fetchrow) {
    if ($row2[0]) {

	my $sth3 = $ensdb->prepare ("select external_id from genedblink where gene_id = '$row2[0]' and external_db = 'MIM'");
	$sth3->execute;
	my $seen3=0;
	while (my @row3 = $sth3->fetchrow) {
	    $seen3=1;

	    #If this gene ac has a link to an OMIM ac dump the text corresponding to this OMIM ac
	    if ($omim{$row3[0]}) {
		print OUT "$row2[0]\|OMIM\|$row3[0]\|$omim{$row3[0]}\n";
	    }
	}
		

	my $sth4 = $ensdb->prepare ("select external_id from genedblink where gene_id = '$row2[0]' and external_db = 'SPTREMBL'");
	$sth4->execute;

	my $seen4=0;
	while (my @row4 = $sth4->fetchrow) {
	    $seen4 = 1;
#Same for SP
	    print OUT "$row2[0]\|SP\|$row4[0]\|$sp{$row4[0]}\n";
	}
    }
}	    
	
#Select all of the gene ac having an intepro domain
my $sth5 = $ensdb->prepare ("select t.gene,i.interpro_ac from protein_feature as pf, transcript as t,interpro as i,interpro_description as id where pf.translation = t.translation and pf.hid = i.id and i.interpro_ac = id.interpro_ac");

$sth5->execute;

while (my @row5 = $sth5->fetchrow) {
#If the gene has an interpro domain, dump the text describing this domain
    if ($hash{$row5[1]}) {
	print OUT "$row5[0]\|INTERPRO\|$row5[1]|$hash{$row5[1]}\n";
    }
    
}

	









