use strict;

=head1 NAME

    grep_sptr_entries.pl

=head1 SYNOPSIS
 
    perl grep_sptr_entries.pl

=head1 DESCRIPTION
    
    The purpose of this script is to grep the right sequences from a sptr flat file (given the ox number), it will also translate these sequences in fasta format.

=head1 CONTACT

    mongin@ebi.ac.uk

=cut

use Bio::SeqIO;


BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
  require "mapping_conf.pl";
}

my $match;

# Get the mapping infos from mapping.conf
my %conf     = %::mapping_conf;

my $input         = $conf{'total_sptr'};
my $swiss_output  = $conf{'sptr_swiss'};
my $fasta_output  = $conf{'sptr_fa'};
my $ox            = $conf{'ox'};

if ((!defined $input) || (!defined $swiss_output) || (!defined $fasta_output) || (!defined $ox)) {
    die "One of the following field hasn't been defined in mapping.conf:\ntotal_sptr: $input\nsptr_swiss: $swiss_output\nsptr_fa: $fasta_output\nOX: $ox\n";
}

print STDERR "INPUT: $input\n";

open (INPUT,"$input");
open (SWISSOUT,">$swiss_output");
open (FASOUT,">$fasta_output"); 

my $pattern = "(\\n\\*\\*OX.* $ox\\;)|(\\nOX.*[\\=\\ 0]$ox)";

# Read an entire trembl record at a time
$/ = "\/\/\n";

while(<INPUT>){
    $match = /$pattern/mo;
    if ($match){
             s/$pattern/\e[7m$&\e[m/mg;
	
            print SWISSOUT;
        }
}    

close (SWISSOUT);
close (INPUT);

#put the default back
$/ = "\n";

my $in  = Bio::SeqIO->new(-file => $swiss_output, '-format' =>'swiss');

while ( my $seq = $in->next_seq() ) {
   print FASOUT ">".$seq->accession."\n".$seq->seq."\n";

}

