#!/usr/local/bin/perl

=head1 NAME - contig2gff

Takes one contig id and provides GFF output
from it

=head1 DESCRIPTION - 

Makes GFF output from a contig id by connecting
to a database, retrieving objects a SeqFeature
objects and dumping using the dump gff routines
provided therein.

=head2 Running the script

Make sure you have the latest (0.06 series) bioperl in your @INC. At
sanger this is at /nfs/disk21/birney/prog/bioperl/bioperl-live

This needs to be run with -I ../modules if you don't have the 
modules in place.

=cut

# dinky script... ;

use strict;
use Bio::EnsEMBL::DB::Obj;

my $db = new Bio::EnsEMBL::DB::Obj( -user => 'root', -db => 'pog' , -host => 'croc' );
my $contig = $db->get_Contig(shift);

foreach my $sf ( $contig->get_all_SeqFeatures ) {
    # $sf is Bio::SeqFeature::Generic object.
    print $sf->gff_string, "\n";
}
