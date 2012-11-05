use strict;
use warnings;

BEGIN {
  $| = 1;
  use Test;
  plan tests => 17;
}

use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Analysis;

our $verbose = 0;    #turn on or off debug statements

#
# Test new and getters
#
my $start      = 10;
my $end        = 100;
my $hstart     = 1;
my $hend       = 90;
my $hstrand    = 1;
my $hseqname   = 'RF1231';
my $percent_id = 90.8;
my $p_value    = '1.52';
my $score      = 50;
my $species    = 'Homo_sapiens';
my $hspecies   = 'Mus_musculus';

my $idesc       = 'interpro description';
my $interpro_ac = 'interpro accession';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test');

my $f = Bio::EnsEMBL::ProteinFeature->new(-START       => $start,
										  -END         => $end,
										  -ANALYSIS    => $analysis,
										  -HSTART      => $hstart,
										  -HEND        => $hend,
										  -HSEQNAME    => $hseqname,
										  -PERCENT_ID  => $percent_id,
										  -P_VALUE     => $p_value,
										  -SCORE       => $score,
										  -SPECIES     => $species,
										  -HSPECIES    => $hspecies,
										  -IDESC       => $idesc,
										  -INTERPRO_AC => $interpro_ac);

ok($f && $f->isa('Bio::EnsEMBL::ProteinFeature'));

ok($f->start == $start);
ok($f->end == $end);
ok($f->analysis == $analysis);

ok($f->hstart == $hstart);
ok($f->hend == $hend);
ok($f->hseqname eq $hseqname);
ok($f->percent_id == $percent_id);
ok($f->p_value == $p_value);
ok($f->score == $score);
ok($f->species  eq $species);
ok($f->hspecies eq $hspecies);

ok($f->idesc       eq $idesc);
ok($f->interpro_ac eq $interpro_ac);

# check that the strand is 0
ok($f->strand == 0);

#
# Test setters
#
ok(test_getter_setter($f, 'idesc',       'interpro desc1'));
ok(test_getter_setter($f, 'interpro_ac', 'interpro ac1'));

