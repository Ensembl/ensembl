package XrefMapper::Methods::ExonerateBest1;

use XrefMapper::Methods::ExonerateBasic;

use vars '@ISA';

@ISA = qw{XrefMapper::Methods::ExonerateBasic};



sub options {

  return ('--bestn', '1');

}


1;
