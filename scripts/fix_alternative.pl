use strict;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Temp;

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'ensembl07';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

my %exons;

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

print STDERR "Getting exon info out of db...\n";
my $sth = $db->prepare("select exon,transcript,rank from exon_transcript");
$sth->execute;
my @out;

while (my $rowhash = $sth->fetchrow_hashref) {
    my $id=$rowhash->{'exon'};
    my $transcript=$rowhash->{'transcript'};
    my $rank=$rowhash->{'rank'};
    my $temp=Bio::EnsEMBL::Temp->new();
    $temp->id($id);
    $temp->transcript($transcript);
    $temp->rank($rank);
    push @out, $temp;
}
@out=sort {$b->transcript <=> $a->transcript} @out;
my $oldtranscript=$out[0]->transcript;
my @ranks;
my $old_clone="";
my @clones;
print STDERR "Finished mysql, checking exons...\n";
foreach my $exon (@out) {
    #print STDERR "Got exon ".$exon->id." on trans ".$exon->transcript."and rank ".$exon->rank."\n";	
    eval {
    if ($exon->transcript eq $oldtranscript) {
	foreach my $rank (@ranks) {
	    if ($exon->rank eq $rank) {
		my $gene_obj=Bio::EnsEMBL::DBSQL::Gene_Obj->new($db);
	        my $real_exon=$gene_obj->get_Exon($exon->id);
		my $clone = $real_exon->clone_id;
		print STDERR "Got clone $clone\n";
		if ($clone ne $old_clone) {
		    push @clones,$clone;
		}
		$old_clone=$clone;
	    }
	}
	push @ranks,$exon->rank;
    }
    else {
	@ranks=[];
    }
    $oldtranscript=$exon->transcript;
};
    if ($@) {
	print STDERR "Could not process ".$exon->id>" because of $@\n";
    }
}
print STDERR "Writing list of buggy clons in bad_clones.list\n";
open (FILE,">bad_clones.list") || die("Could not open bad_clones.list");
foreach my $clone (@clones) {
    print FILE $clone."\n";
}
close FILE;
