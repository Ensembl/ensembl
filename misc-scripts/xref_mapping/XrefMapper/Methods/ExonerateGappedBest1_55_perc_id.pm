package XrefMapper::Methods::ExonerateGappedBest1_55_perc_id;

use XrefMapper::Methods::ExonerateBasic;

use vars '@ISA';

@ISA = qw{XrefMapper::Methods::ExonerateBasic};



sub options {
  return ('--gappedextension FALSE', '--model', 'affine:local', '--subopt', 'no', '--bestn', '1');
}

sub query_identity_threshold {
  return 55;
}

sub target_identity_threshold {
  return 55;
}


1;

