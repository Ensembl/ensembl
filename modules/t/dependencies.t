
use strict;
use warnings;


use Bio::EnsEMBL::Test::TestUtils;


BEGIN { $| = 1;
	use Test;
	plan tests => 1;
}

#test for dependencies on Variation, Compara and Funcgen APIs

my @result = `egrep -r "use Bio::EnsEMBL::(Variation|Compara|Funcgen){1}" ../Bio/EnsEMBL/`;

my %result = map{$_ => 1} @result;

my @exceptions = ('/Bio/EnsEMBL/Utils/TranscriptAlleles.pm', 
   	       	 '/Bio/EnsEMBL/Utils/ensembl_init.example');

my $exceptions = join("|",@exceptions);
$exceptions =~ s/\//\\\//g;


foreach my $key (keys %result) {
	if ( $key =~ /($exceptions)/ ) {
	   delete($result{$key});
	}
}

ok(!%result);

if (%result) { 
   warn "Dependencies found in the following files:\n";
   warn keys %result;
}
