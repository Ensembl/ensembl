use strict;

=head1 dump_textview

=head2 Description

The aim of this script is to produce a flat file containing text information about gathered from different databases (OMIM, Interpro, Swissprot). For OMIM and SP the link from the Ensembl database to the external one is taken from GeneDBLink. For Interpro information is taken from the interpro table.

=head2 Synopsis

dump_textview_dump.pl -db ens075 -host ecs1c.sanger.ac.uk -dbuser ensro

=head2 Format of the file dumped

ENSG....|ExternalDB|ExternalAC|Text

=head2 Options

    -db          dbname (dbname= in locator)
    -host        host name (gets put as host= in locator)
    -dbuser      what username to connect to the database (dbuser= in locator)
    -interpro    name of the file which contains the interpro data (XML file)
    -omim        name of the file for OMIM
    -sp          name of the file for SP
    -domainout   name of the file where the output for domains indexing should be written 
    -extout      name of the file where the output for external DB indexing should be written

=cut


use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $dbpass = undef;
my $dbuser = 'ensro';
my $ensdbname = 'mus_musculus_core_5_3';
my $host = 'ecs1d';
my $output;
my $interpro;
my $omim;
my $sp;
my $domainout;

&GetOptions(
	    'db:s' => \$ensdbname,
	    'host:s'=> \$host,
	    'dbuser:s'=> \$dbuser,
	    'interpro:s'=>\$interpro,
	    'omim:s'=>\$omim,
	    'sp:s'=>\$sp,
	    'domainout:s'=>\$domainout,
	    'extout:s' => \$output
	    );


my $enslocator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;dbname=$ensdbname;user=$dbuser;pass=$dbpass;perlonlyfeatures=1";
my $ensdb =  Bio::EnsEMBL::DBLoader->new($enslocator);

my %omim;
my %sp;


#In this file will be printed out the results of the parsing of Interpro flat file (domains)
open (DOMAINOUT, ">$domainout") || die "Can't open $domainout\n";

#In this file will be printed out the results of the parsing of the external DB (eg: Swiss-Prot, OMIM, ...)
open (OUT, ">$output") || die "Can't open $output\n";

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


#For each entry get the abstract,publist and memberlist.
#Each fiels beginns by <fieldname> and ends by <\fieldname>

my $interpro_id;
my $interpro_abstract;
my $interpro_pub;
my $interpro_member;

ENTRY: while (<INTERPRO>) {
    if (/\<interpro id\=\"(\S+)\"/) {
	$interpro->{id}=$1;
	#print "$1\n";
    }
    my $abstract;
    
    if (/\<abstract/) {
	my $string=$_;
	while (<INTERPRO>) {

	    $string .= $_;
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
	$interpro->{abstract}=$string;
    }

    if (/\<publist/) {
	my $string = $_;
	while (<INTERPRO>) {
	    $string .= $_;
	    if (/\<\/publist|\<memberlist/) {
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
		 
    if (/\<\/interpro\>/) {
	
	my $inter_id = $interpro->{id};
	my $total =  $interpro->{abstract}." ". $interpro->{pub}." ". $interpro->{member},"\n";
	$hash{$inter_id} = $total;
	#print "$inter_id\t$total\n\n";
	#$interpro->{id}=$1;
	#$count++;
	my $check = $total;
	my $in;
	#while ($check =~ s/interpro id\=//) {
	 #   $in++;
	#}
	#if ($in > 1 ) {
	 #   print "Found $in in $total\n";
	#}
     
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
my $in  = Bio::SeqIO->new(-file => $sp, '-format' =>'swiss');

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

    #foreach my $comment ( $seq->annotation->each_Comment ) {
	#push (@comments,$comment->text);
	#}

    #my $tot_sp = "$ke|$desc|@refs|@comments";
 
    my $tot_sp = "$ke|$desc|@refs";
 

    if (!defined $sp{$ac}) {
	$sp{$ac} = [];
    }

    $sp{$ac} = $tot_sp;

}

#Dump files
print STDERR "Dumping file\n";

#Get all of the distinct gene ac present in genedblink
#my $sth2 = $ensdb->prepare ("select distinct (gene_id) from genedblink");

my $sth2 = $ensdb->prepare ("select distinct (ensembl_id) from objectXref");

$sth2->execute;

while (my @row2 = $sth2->fetchrow) {
    if ($row2[0]) {
	#print STDERR "$row2[0]\n";
#	print STDERR "Dump for MIM\n";

	#my $sth3 = $ensdb->prepare ("select external_id from genedblink where gene_id = '$row2[0]' and external_db = 'MIM'");

	my $sth3 = $ensdb->prepare ("select x.dbprimary_id, gs.stable_id from gene_stable_id gs, Xref as x, externalDB as e, objectXref as o, transcript as t where gs.gene_id=t.gene_id and o.xrefId = x.xrefId and o.ensembl_id = '$row2[0]' and t.translation_id = '$row2[0]'");

	$sth3->execute;
	my $seen3=0;
	while (my @row3 = $sth3->fetchrow) {
	    $seen3=1;

	    #If this gene ac has a link to an OMIM ac dump the text corresponding to this OMIM ac
	    if ($omim{$row3[0]}) {
		print OUT "$row3[1]\|MIM\|$row3[0]\|$omim{$row3[0]}\n";
	    }

	    if ($sp{$row3[0]}) {

		print OUT "$row3[1]\|SPTR\|$row3[0]\|$sp{$row3[0]}\n";
	    }
	}
		
	#print STDERR "Dump for MIM\n";

	#my $sth4 = $ensdb->prepare ("select external_id from genedblink where gene_id = '$row2[0]' and (external_db = 'SPTREMBL' or external_db = 'SWISS')");

#my $sth4 = $ensdb->prepare ("select x.dbprimary_id, t.gene from Xref as x, externalDB as e, objectXref as o, transcript as t where o.xrefId = x.xrefId and o.ensembl_id = '$row2[0]' and t.translation = '$row2[0]'");

	#$sth4->execute;

	#my $seen4=0;
	#while (my @row4 = $sth4->fetchrow) {
	    #$seen4 = 1;
#Same for SP
	    #if ($sp{$row4[0]}) {

		#print OUT "$row4[1]\|SPTR\|$row4[0]\|$sp{$row4[0]}\n";
	    #}
	#}

    }
}	    
	


print STDERR "Dump for Interpro\n";

#Select all of the gene ac having an intepro domain
my $sth5 = $ensdb->prepare ("select gs.stable_id,i.interpro_ac from protein_feature as pf, gene_stable_id gs,transcript as t,interpro as i where gs.gene_id=t.gene_id and pf.translation = t.translation_id and pf.hid = i.id");

$sth5->execute;

my %saw;
while (my @row5 = $sth5->fetchrow) {
#If the gene has an interpro domain, dump the text describing this domain
  
    #if ($row5[1] eq "IPR002348") {
	#print STDERR "Got IPR002348\t$hash{IPR00248}\n";

#}
    my $ke = "$row5[1]:$hash{$row5[1]}";
    if ($hash{$row5[1]}) {
	if (!defined $saw{$ke}){
	    print DOMAINOUT "$row5[0]\|INTERPRO\|$row5[1]|$hash{$row5[1]}\n";
	    $saw{$ke}=1;
	}
    }
    else {
	print "$row5[1]\n";
    }
}    
	

	









