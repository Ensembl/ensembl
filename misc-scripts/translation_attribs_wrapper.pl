#!/software/bin/perl

=head1 NAME

translation_attribs_wrapper.pl - script to calculate peptide statistics, if the first aminoacid is
                                methionine and there is a stop codon in the aminoacid sequence
                                and store them in translation_attrib table. This is mainly a wrapper
                                around the translation_attribs script to submit several jobs to the farm

=head1 SYNOPSIS

translation_attribs_wrapper.pl [arguments]

Required arguments:

  --user=user                         username for the database

  --pass=pass                         password for database

  --release=release                   release number

Optional arguments:

  --binpath=PATH                      directory where the binary script to calculate 
                                      pepstats is stored (default: /software/pubseq/bin/emboss)

  --tmpfile=file                      file to store tmp results of pepstats (default=/tmp)

  --host=host                         server where the core databases are stored (default: ens-staging)

  --port=port                         port (default=3306)

  --pepstats_only                     when used, will only run the pepstats calculation

  --met_and_stop_only                 when used, will only run the methionine and stop codon calculation
  
  --help                              print help (this message)

=head1 DESCRIPTION

This script will calculate the peptide statistics, if the first aminoacid is methionine and
there is a stop codon  for all core databases in the server and store them as a 
translation_attrib values. This is a wraper around the translation_attrib and will simply
submit jobs to the farm grouping the core databases in patterns

=head1 EXAMPLES

Calculate translation_attributes for all databases in ens-staging 

  $ ./translation_attribs_wrapper.pl --user ensadmin --pass password --release 51

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Daniel Rios <dani@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Utils::Exception qw(throw);


## Command line options

my $binpath = '/software/pubseq/bin/emboss'; 
my $tmpdir = '/tmp';
my $host = 'ens-staging';
my $release = undef;
my $user = undef;
my $pass = undef;
my $port = 3306;
my $help = undef;
my $pepstats_only = undef;
my $met_and_stop_only = undef;

GetOptions('binpath=s' => \$binpath,
	   'tmpdir=s' => \$tmpdir,
	   'host=s'    => \$host,
	   'user=s'    => \$user,
	   'pass=s'    => \$pass,
	   'port=s'    => \$port,
	   'release=i' => \$release,
	   'pepstats_only' => \$pepstats_only,
	   'met_and_stop_only' => \$met_and_stop_only,
	   'help'    => \$help
	   );

pod2usage(1) if($help);
throw("--user argument required") if (!defined($user));
throw("--pass argument required") if (!defined($pass));
throw("--release argument required") if (!defined($release));

my $queue = 'long';
my $memory = "'select[mem>4000] rusage[mem=4000]' -M4000000";
my $options = '';
if (defined $binpath){
    $options .= "--binpath $binpath ";
}
if (defined $tmpdir){
    $options .= "--tmpdir $tmpdir "
}
if (defined $host){
    $options .= "--host $host ";
}
if (defined $port){
    $options .= "--port $port ";
}
if (defined $pepstats_only){
    $options .= "--pepstats_only ";
}
if (defined $met_and_stop_only){
    $options .= "--met_and_stop_only "
}
my @ranges = ('^[a-b]','^[c-d]','^[e-f]','^[g-h]','^[i-m]','^[n-p]','^[q-r]','^[s-t]','^[u-z]');
my $core_db = ".*core_$release\_.*";
my $call;
foreach my $pattern (@ranges){
    $call = "bsub -o /lustre/scratch1/ensembl/dr2/tmp_smallfiles/output_translation_$pattern.txt -q $queue -R$memory ./translation_attribs.pl -user $user -pass $pass $options";
    $call .= " --pattern '" . $pattern . $core_db. "'";

    system($call);
#print $call,"\n";
}

#we now need to run it for the otherfeatures|vega databases, but only the pepstats

my $vega_db = ".*_vega_$release\_.*";
$call = "bsub -o /lustre/scratch1/ensembl/dr2/tmp_smallfiles/output_translation_vega.txt -q $queue -R$memory ./translation_attribs.pl -user $user -pass $pass $options";
$call .= " --pattern '" . $vega_db. "'";

system($call);
#print $call,"\n";

$options .= "--pepstats_only ";

my $other_db = ".*_otherfeatures_$release\_.*";
$call = "bsub -o /lustre/scratch1/ensembl/dr2/tmp_smallfiles/output_translation_other.txt -q $queue -R$memory ./translation_attribs.pl -user $user -pass $pass $options";
$call .= " --pattern '" . $other_db. "'";

system($call);
#print $call,"\n";
