use strict;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Temp;

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'ensembl';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

my %exons;

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my $sth = $db->prepare("select exon,transcript,rank from exon_transcript where transcript='ENST00000043895'");
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

foreach my $exon (@out) {
    #print STDERR "Got exon ".$exon->id." on trans ".$exon->transcript." and rank ".$exon->rank."\n";
    eval {
    if ($exon->transcript eq $oldtranscript) {
	foreach my $rank (@ranks) {
	    if ($exon->rank eq $rank) {
		my $real_exon=$db->get_Exon($exon->id);
		my $clone = $real_exon->clone_id;
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

open (FILE,">bad_clones.list") || die("Could not open bad_clones.list");
foreach my $clone (@clones) {
    print FILE $clone."\n";
}
close FILE;
