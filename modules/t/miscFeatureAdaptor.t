use strict;

use lib 't';
use TestUtils qw(test_getter_setter debug);

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}

use MultiTestDB;

our $verbose = 1;

my $multi = MultiTestDB->new;

# get a core DBAdaptor
my $dba = $multi->get_DBAdaptor("core");


#
# Test get_MiscFeatureAdaptor works
#
my $mfa = $dba->get_MiscFeatureAdaptor();

ok($mfa && ref($mfa) && $mfa->isa('Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor'));


#
# Test fetching by slice
#


my $chr_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', '20');
my $features = $mfa->fetch_all_by_Slice($chr_slice);
debug('--- chr 20 misc_features ---');
debug("Got " . scalar(@$features));
ok(@$features == 7);
print_features($features);


$features = $mfa->fetch_all_by_Slice_and_set_code($chr_slice,'ntctgs');
debug('--- chr 20 ntcontigs set---');
debug("Got " . scalar(@$features));
ok(@$features == 7);
print_features($features);

$features = $mfa->fetch_all_by_Slice_and_set_code($chr_slice,'cloneset');
debug('--- chr 20 cloneset set---');
debug("Got " . scalar(@$features));
ok(@$features == 0);
print_features($features);

my $feature = $mfa->fetch_by_dbID(741);
debug('--- fetch_by_dbID ---');
ok($feature->dbID() == 741);
ok($feature->start() == 29307627);
ok($feature->end() == 29541965);
ok($feature->strand() == 1);
print_features([$feature]);



sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
    if(defined($f)) {
      my @attrib_types = $f->get_attribute_types;
      my @attrib_vals  = map {"$_ => " . join(',',$f->get_attribute($_))} @attrib_types;

      my @set_codes = $f->get_set_codes;
      my @set_names = map {"$_ => " . $f->get_set($_)->name()} @set_codes;

      my $attrib_string = join(":", @attrib_vals);
      my $set_string = join(":", @set_names);

      my $seqname = $f->slice->seq_region_name();
      debug($seqname . ' ' . $f->start().'-'.$f->end().'('.$f->strand().
            ') ['. $f->dbID.'] '." ($attrib_string) ($set_string)");
    } else {
      debug('UNDEF');
    }
  }
}
