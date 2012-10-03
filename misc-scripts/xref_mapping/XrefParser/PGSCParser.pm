package XrefParser::PGSCParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

sub run {

 my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $gene_source_id = $self->get_source_id_for_source_name("PGSC_GENE");

  my $pgsc_io = $self->get_filehandle($file);

  if ( !defined $pgsc_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  while ( $_ = $pgsc_io->getline() ) {

    chomp;

    if ($_ =~ /^([^\t]+)\t([^\t]+)$/) {
	my ($pgsc_id, $desc) = ($1,$2);
	
	if ($verbose) {
	    #print STDERR "pgsc_id, $pgsc_id\n";
	    #print STDERR "desc: $desc\n";
	}
	
	my $locus_name = $pgsc_id;
	
	my $gene_xref_id = $self->add_xref({ acc        => $pgsc_id,
					     label      => $locus_name,
					     desc       => $desc,
					     source_id  => $gene_source_id,
					     species_id => $species_id,
					     info_type  => "DIRECT"} );
	$self->add_direct_xref($gene_xref_id, $pgsc_id, "Gene", "DIRECT");
	
	$xref_count++;
	
    }
    else {
	if ($verbose) {
	    print STDERR "failed to parse line, $_\n\n";
	}
    }
  }
 
 $pgsc_io->close();
 
 print $xref_count." PGSC Xrefs added\n" if($verbose);
 return 0; #successful
}

1;
