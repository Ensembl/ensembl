use strict;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

package ExternalFF1;
  use Bio::EnsEMBL::External::ExternalFeatureAdaptor;
  use Bio::EnsEMBL::Test::MultiTestDB;
  use Bio::EnsEMBL::Feature;

  use vars qw(@ISA);
  @ISA = qw(Bio::EnsEMBL::External::ExternalFeatureAdaptor);

  sub coordinate_systems {
    return ('CONTIG');
  }

  sub ensembl_db {
    my $self = shift;
    return $multi->get_DBAdaptor("core");
  }

  sub fetch_all_by_contig_name {
    my $self = shift;
    my $name = shift;
    my $start = shift;
    my $end   = shift;

    my @out;

    push @out, Bio::EnsEMBL::Feature->new(-START => 10_000,
                                         -END   => 11_000,
                                         -STRAND => 1);

    return \@out;
  }




package ExternalFF2;
  use Bio::EnsEMBL::External::ExternalFeatureAdaptor;
  use Bio::EnsEMBL::Test::MultiTestDB;
  use Bio::EnsEMBL::Feature;

  use vars qw(@ISA);
  @ISA = qw(Bio::EnsEMBL::External::ExternalFeatureAdaptor);

  sub coordinate_systems {
    return ('CLONE');
  }

  sub ensembl_db {
    my $self = shift;
    return $multi->get_DBAdaptor("core");
  }

  sub fetch_all_by_clone_accession {
    my $self = shift;
    my $name = shift;
    my $version = shift;
    my $start = shift;
    my $end   = shift;

    my @out;

    push @out, Bio::EnsEMBL::Feature->new(-START => 10_000,
                                         -END   => 11_000,
                                         -STRAND => 1);

    return \@out;
  }




package Test;

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;


my $xff = ExternalFF1->new();

my $db = $multi->get_DBAdaptor('core');

my $slice_adaptor= $db->get_SliceAdaptor();
my $slice = $slice_adaptor->fetch_by_region('chromosome', '20', 30e6,35e6);

my $feats = $xff->fetch_all_by_Slice($slice);

foreach my $f (@$feats) {
  my ($start, $end, $strand) = ($f->start, $f->end, $f->strand);
  debug("F: $start-$end ($strand)\n");
}

ok(@$feats == 12);

my $cln_slice = $slice_adaptor->fetch_by_region('clone', 'AL359765.6');

$feats = $xff->fetch_all_by_Slice($cln_slice);

foreach my $f (@$feats) {
  my ($start, $end, $strand) = ($f->start, $f->end, $f->strand);
  debug("F: $start-$end ($strand)\n");
}

ok(@$feats == 1);

$xff = ExternalFF2->new();

$feats = $xff->fetch_all_by_Slice($slice);


foreach my $f (@$feats) {
  my ($start, $end, $strand) = ($f->start, $f->end, $f->strand);
  debug("F: $start-$end ($strand)\n");
}

ok(@$feats == 12);

$feats = $xff->fetch_all_by_Slice($cln_slice);

foreach my $f (@$feats) {
  my ($start, $end, $strand) = ($f->start, $f->end, $f->strand);
  debug("F: $start-$end ($strand)\n");
}

ok(@$feats == 1);

$feats = $xff->fetch_all_by_supercontig_name('NT_028392');

foreach my $f (@$feats) {
  my ($start, $end, $strand) = ($f->start, $f->end, $f->strand);
  debug("F: $start-$end ($strand)\n");
}

ok(@$feats == 12);

$multi=undef;

done_testing();
