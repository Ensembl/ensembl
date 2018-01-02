=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


package XrefParser::Broad_GenericParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

#
# Broad gene annotation file with gene predictions, their gene names and descriptions
#
#
# This is the based Broad parser, with the logic code for parsing, but doesn't know anything about the URL and the source. Needs to add an extra module that inherits from it with the source logic.
#

sub parse {
    
    my ($self, $ref_arg) = @_;
    my $species_id   = $ref_arg->{species_id};
    my $source_id    = $ref_arg->{source_id};
    my $files        = $ref_arg->{files};
    my $verbose      = $ref_arg->{verbose};
    
    if( (!defined $species_id) or (!defined $files) ){
	croak "Need to pass source_id, species_id and files as pairs";
    }
    $verbose |=0;
    
    my $file = @{$files}[0];
    
    # Very hacky code, but can't get better!
    # Problem is that the filename from the URI is 'download', which turns up to be compressed
    # so need to uncompress it after download.
    # Because it doesn't match .gz, not done yet automatically, so we have to do it now

    # Other problem is that should it not try to uncompress it when already uncompressed (ie when not downloaded again)
    # Would be nice to capture the error code when trying to uncompress an already uncompressed file, 
    # but haven't figured out yet how to !
    # To do so, it seems that somehow we need to capture STDERR
    
    system ("mv $file $file.gz");
    print STDERR "file: $file\n";
    eval {
	system ("gunzip", "-f", "$file.gz");
    };
    if ($@) {
	# doesn't catch anything actually, why ?
	print STDERR "catching eval error code!\n";
    }
    system ("mv $file.gz $file");

    my @xrefs = ();
    
    if ($verbose) {
	print STDERR "Parsing broad source, $source_id\n";
    }

    my $file_io = $self->get_filehandle($file);
    
    if ( !defined $file_io ) {
	print STDERR "ERROR: Could not open $file\n";
	return 1;    # 1 is an error
    }
    
    while ( $_ = $file_io->getline() ) {
	
	my $line = $_;
	chomp $line;

	next if ($line =~ /^LOCUS/);   # skip header
	
	my $xref = {};
	
	my ($stable_id,$symbol,$synonym,$length,$start,$end,$strand,$desc,$chromosome) = $line =~ /^([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t*.*/ or warn("Can't parse Broad entry: $_\n"); 

	# build the xref object and store it
	$xref->{ACCESSION}     = $stable_id;
	$xref->{LABEL}         = $stable_id;
	$xref->{DESC}          = $desc;
	$xref->{SPECIES_ID}    = $species_id;
	$xref->{SOURCE_ID}     = $source_id;
	
	push @xrefs, $xref;
	
    }

    $file_io->close();
    
    print STDERR scalar(@xrefs) . " BROAD_ xrefs succesfully parsed\n" if($verbose);
    
    return \@xrefs;
}

sub store {
    my ($self, $xrefs_aref) = @_;

    foreach my $xref_href (@$xrefs_aref) {

	my $gene_xref_id = $self->add_xref({ acc         => $xref_href->{ACCESSION},
					     label       => $xref_href->{LABEL},
					     desc        => $xref_href->{DESC},
					     source_id   => $xref_href->{SOURCE_ID},
					     species_id  => $xref_href->{SPECIES_ID},
					     version     => 1,
					     info_type   => "DIRECT"} );
	$self->add_direct_xref($gene_xref_id, $xref_href->{ACCESSION}, "Gene", "DIRECT");
    }

    return 0; # success

}

1;
