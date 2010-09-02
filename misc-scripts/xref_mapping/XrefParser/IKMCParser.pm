package XrefParser::IKMCParser;

use strict;

use base qw( XrefParser::BaseParser );

# This parser will read Direct Xrefs from a simple tab-delimited file.
# The columns of the file should be the following:
#
# 1)    Accession ID
# 2)    label
# 3)    source type
# 4)    stable_id
#

sub new {
    my $proto = shift;

    my $class = ref $proto || $proto;
    my $self = bless {}, $class;

    return $self;
}

sub run {
    my $self = shift;

    my $source_id = shift;
    my $species_id = shift;
    my $files_ref  = shift;
    my $rel_file   = shift;
    my $verbose = shift;
    
    my $filename = @{$files_ref}[0];

    my $file_io = $self->get_filehandle($filename);
    if ( !defined($file_io) ) {
        return 1;
    }

    my $parsed_count = 0;

    printf( STDERR "source = %d\t species = %d, file is %s\n",
            $source_id, $species_id, $filename );

    my %type2id;
    foreach my $t ("ES cells available", "Vector available", "No products available yet", "Mice available"){
      my $ikmc = "IKMC_".$t;
      $ikmc =~ s/ /_/g;
      $type2id{$t}  = XrefParser::BaseParser->get_source_id_for_source_name($ikmc);
      print $ikmc."\t".$type2id{$t}."\n";
      if(!defined( $type2id{$t})){
	die  "Could not get source id for $ikmc\n";
      }
    }	


    while ( defined( my $line = $file_io->getline() ) ) {
        chomp $line;

        my ( $accession, $label, $source_type, $ensembl_id)
          = split( /\t/, $line );

        if ( !defined($accession)) {
            printf( "Line %d contains  has less than one column.\n",
                    1 + $parsed_count );
            print("The parsing failed\n");
            return 1;
        }

        my $type        = 'gene';
        $label       ||= $accession;

	my $source_id = $type2id{$source_type};
        ++$parsed_count;

        my $xref_id =
          XrefParser::BaseParser->get_xref( $accession, $source_id, $species_id );

        if ( !defined($xref_id) || $xref_id eq '' ) {
            $xref_id =
              XrefParser::BaseParser->add_xref(
                                   $accession,   undef,   $label,
                                   '', $source_id, $species_id, "DIRECT"
              );
        }
        XrefParser::BaseParser->add_direct_xref( $xref_id, $ensembl_id,
                                                 $type, $accession );
    } ## end while ( defined( my $line...

    printf( "%d direct xrefs succesfully parsed\n", $parsed_count );

    $file_io->close();

    print "Done\n";

    return 0;
} ## end sub run

1;
