#!/usr/local/bin/perl

use strict;

=head1 NAME - clone2contig

Loops over each argv giving you out the contigs contained within it

=head1 DESCRIPTION - 

=head2 Running the script

Make sure you have the latest (0.06 series) bioperl in your @INC. At
sanger this is at /nfs/disk21/birney/prog/bioperl/bioperl-live

This needs to be run with -I ../modules if you don't have the 
modules in place.

=cut

use Bio::EnsEMBL::DB::Obj;

my $db;
eval {
    $db = new Bio::EnsEMBL::DB::Obj( -user => 'root', -db => 'pog' , -host => 'croc.sanger.ac.uk');
};

if( $@ ) {
    die "Could not connect to database\n$@\n";
}

foreach my $clone ( @ARGV ) {

    my @contigs;


    my $clone = $db->get_Clone($clone);
    @contigs = $clone->get_all_Contigs();


    foreach my $contig ( @contigs ) {
	print $contig->id(), "\n";
    }
}


