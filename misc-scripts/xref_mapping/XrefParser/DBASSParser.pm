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

package XrefParser::DBASSParser;

use strict;

use base qw( XrefParser::BaseParser);
use DBI;
use Carp;

# This parser will read direct xrefs from a simple comma-delimited file downloaded from the DBASS database.
# The columns of the file should be the following:
#
# 1)    DBASS Gene ID
# 2)    DBASS Gene Name
# 3)    Ensembl Gene ID


my $dbi;
my $xref_id;
my $source_id;


sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $filename = @{$files}[0];


    my $file_io = $self->get_filehandle($filename);
    if ( !defined($file_io) ) {
        return 1;
    }

     my $parsed_count = 0;

     $file_io->getline();
    
    while ( defined( my $line = $file_io->getline() ) ) {
 
	$line =~ s/\s*$//;
        # csv format can come with quoted columns, remove them
        $line =~ s/"//g;

        my ( $dbass_gene_id, $dbass_gene_name, $ensembl_id) = split( /,/, $line );

        if ( !defined($dbass_gene_id) || !defined($ensembl_id) ) {
            printf STDERR ( "Line %d contains  has less than two columns.\n",
                    1 + $parsed_count );
            print STDERR ("The parsing failed\n");
            return 1;
        }
	
	
	my $first_gene_name = $dbass_gene_name;
	my $second_gene_name;
	
	
	if ($dbass_gene_name =~ /.\/./){
		($first_gene_name, $second_gene_name) = split( /\//, $dbass_gene_name );
	}
	
	if ($dbass_gene_name =~ /(.*)\((.*)\)/){
		$first_gene_name = $1;
	 	$second_gene_name = $2;
#		print $first_gene_name, "\n", $second_gene_name, "\n" if($verbose);
	}
       
       
        my $label       = $first_gene_name;
	my $type        = 'gene';
        my $description = '';
        my $version     = '1';
	my $synonym	= $second_gene_name;

        ++$parsed_count;

        my $xref_id = $self->get_xref( $dbass_gene_id, $source_id, $species_id, $dbi );

        if ( !defined($xref_id) || $xref_id eq '' ) {
            $xref_id = $self->add_xref({ acc        => $dbass_gene_id,
					 version    => $version,
					 label      => $label,
					 desc       => $description,
					 source_id  => $source_id,
                                         dbi        => $dbi,
					 species_id => $species_id,
					 info_type => "DIRECT"} );
        }
	
	$self->add_direct_xref( $xref_id, $ensembl_id, $type, '', $dbi);
	
	if (defined ($synonym)) {
		$self->add_synonym($xref_id, $synonym, $dbi);
	}
	elsif ($synonym =~ /^\s/){
		print "There is white space! \n" if($verbose);	
	}
	else {
		next;
	}
	
	
    } ## end while ( defined( my $line...

    printf( "%d direct xrefs succesfully parsed\n", $parsed_count ) if($verbose);

    $file_io->close();


    return 0;
} ## end sub run


1;
