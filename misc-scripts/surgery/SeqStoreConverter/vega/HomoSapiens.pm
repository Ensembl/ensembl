package SeqStoreConverter::vega::HomoSapiens;

use strict;
use warnings;

use SeqStoreConverter::HomoSapiens;
use SeqStoreConverter::vega::VBasicConverter;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::vega::VBasicConverter SeqStoreConverter::HomoSapiens);

1;
