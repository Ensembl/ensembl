use Test;
use strict;

BEGIN { $| = 1; plan tests => 9 }

use Bio::EnsEMBL::Test::MultiTestDB;

ok(1);

# Database will be dropped when this
# object goes out of scope
my $ens_test = Bio::EnsEMBL::Test::MultiTestDB->new;

ok($ens_test);

my $dba = $ens_test->get_DBAdaptor("core");

ok($dba);

my $sth = $dba->dbc->prepare("select * from gene");
$sth->execute;

ok(scalar($sth->rows) == 20);


# now hide the gene table i.e. make an empty version of it
$ens_test->hide("core","gene");
$sth->execute;
ok($sth->rows == 0);


# restore the gene table
$ens_test->restore();
$sth->execute;
ok(scalar($sth->rows) == 20);


# now save the gene table i.e. make a copy of it
$ens_test->save("core","gene");
$sth->execute;
ok(scalar($sth->rows) == 20);


# delete 9 genes from the db
$sth = $dba->dbc->prepare("delete from gene where gene_id >= 18266");
$sth->execute;

$sth = $dba->dbc->prepare("select * from gene");
$sth->execute;

ok(scalar($sth->rows) == 10);


# check to see whether the restore works again
$ens_test->restore();
$sth->execute;
ok(scalar($sth->rows) == 20);


$sth->finish;


