
use strict;
use warnings;

use Test::More;
use File::Spec;
#use Bio::EnsEMBL::Test::TestUtils;

my ($volume, $directory, $file) = File::Spec->splitpath(__FILE__);
$directory = File::Spec->rel2abs($directory);
my $modules_dir = File::Spec->catdir($directory, File::Spec->updir(), qw/Bio EnsEMBL/);

#test for dependencies on Variation, Compara and Funcgen APIs
my @result = `egrep -r "use Bio::EnsEMBL::(Variation|Compara|Funcgen){1}" $modules_dir`;

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

done_testing();