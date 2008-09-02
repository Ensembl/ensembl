package XrefParser::VbDirectParser;

use strict;
use File::Basename;
use Bio::SeqIO;
use XrefParser::BaseParser;


use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);

# OParser for FASTA-format probe mappings from Agilent
# >A_23_P253586
# CTGTCAGGATTCTAGAACTTCTAAAATTAAAGTTTGGGGAAATCAGTAGCTCTGATGAGA
# >A_23_P217507
# AGAAAGACGTTTTCCAACATGTAGAACTGCTTTTTAACTGGAGGAAAAATACTTCAGGAG

sub run {

	my ($self, $file, $source_id, $species_id) = @_;

	open (INFILE, "<$file");
	my $i=0;
	print STDERR "source: ".$source_id."\tspecies: ".$species_id."\n";
	

	my $type = "transcript";

	while (my $ln = <INFILE>) {
		chomp($ln);
		my ($probe,$id, $version, $description, $ensembl_id) = split("\t",$ln);      
		$i++;

		my $xref_id = XrefParser::BaseParser->get_xref($probe, $source_id, $species_id);
		if (!defined($xref_id) || $xref_id eq "") {
			$xref_id = XrefParser::BaseParser->add_xref($probe, 1, $probe, $description, $source_id, $species_id);
			}
		XrefParser::BaseParser->add_direct_xref($xref_id, $ensembl_id, $type, $probe);
		}

	print $i." VB direct xrefs succesfully parsed\n" if($verbose);

	close(INFILE);

	return 0;

	}


sub new {
  my $self = {};
  bless $self, "XrefParser::VbDirectParser";
  return $self;

}

1;
