use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;

my $dbpass = undef;
my $dbuser = 'ensro';
my $ensdbname = 'ensembl080';
my $host = 'hcs2g.sanger.ac.uk';

my $enslocator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;dbname=$ensdbname;user=$dbuser;pass=$dbpass;perlonlyfeatures=1";
my $ensdb =  Bio::EnsEMBL::DBLoader->new($enslocator);


open (INPUT,"textview.dump") || die "Can't open File\n";
open (OUT,">new_dump.txt") || die "Can't open File\n";

my %seen;

while (<INPUT>) {
    my ($ens,$db,$ac,@text) = split(/\|/,$_);

    if ($db eq "INTERPRO") {
	if (!defined $seen{$ac}) {
	    $seen{$ac} = 1;
	    my $sth5 = $ensdb->prepare ("select distinct(t.gene) from protein_feature as pf, transcript as t,interpro as i where pf.translation = t.translation and pf.hid = i.id and i.interpro_ac = '$ac'");
	    
	    
	    
	    $sth5->execute;
	    my @genes;
	    
	    while (my @row5 = $sth5->fetchrow) {
		
		push (@genes,$row5[0]);
	    }
	


	my $total = join(";",@genes);
	print OUT "$total|$db|$ac|@text";
	}

    }

    else {
	print OUT $_;
    }
}
