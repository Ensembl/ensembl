package XrefMapper::Methods::ExonerateGappedBest1;

use XrefMapper::Methods::ExonerateBasic;

use vars '@ISA';

@ISA = qw{XrefMapper::Methods::ExonerateBasic};



sub options {

  return ('--model', 'affine:local', '--subopt', 'no', '--bestn', '1');

}

sub query_identity_threshold {

  return 100;

}

sub target_identity_threshold {

  return 100;

}


1;
