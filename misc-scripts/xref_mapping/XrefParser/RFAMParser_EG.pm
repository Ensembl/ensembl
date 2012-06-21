package XrefParser::RFAMParser_EG;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;
use Bio::EnsEMBL::Registry;

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  #get direct RFAM xrefs from core

  my $registry = "Bio::EnsEMBL::Registry";

  $registry->load_registry_from_multiple_dbs( 
      {
        '-host'     => 'mysql-eg-staging-1.ebi.ac.uk',
	'-port'     => 4160,
        '-user'     => 'ensro',
      },
      {
        '-host'     => 'mysql-eg-staging-2.ebi.ac.uk',
	'-port'     => 4275,
        '-user'     => 'ensro',
      },
 
  );

  #get the species name
  my %id2name = $self->species_id2name;
  my $species_name = $id2name{$species_id}[0];

  my $dba = $registry->get_DBAdaptor($species_name, 'core');

  my $rfam_sql = "select distinct t.stable_id, hit_name from analysis a join transcript t on (a.analysis_id = t.analysis_id and a.logic_name = 'ncRNA' and t.biotype != 'miRNA') join exon_transcript et on (t.transcript_id = et.transcript_id) join supporting_feature sf on (et.exon_id = sf.exon_id and sf.feature_type = 'dna_align_feature' ) join dna_align_feature df on (sf.feature_id = df.dna_align_feature_id) order by hit_name";

  my $sth = $dba->dbc->prepare($rfam_sql);
  $sth->execute();

  #hash keyed on RFAM accessions, value is an array of ensembl transcript stable_ids
  my %rfam_transcript_stable_ids;

  while (my ($stable_id, $hit_name) = $sth->fetchrow_array ) {

      my $rfam_id;
      if ($hit_name =~ /^(RF\d+)/) {
	  $rfam_id = $1;
      }
      if ($rfam_id) {
	  push @{$rfam_transcript_stable_ids{$rfam_id}}, $stable_id; 
      }
  }
  $sth->finish;     

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my @xrefs;

  local $/ = "//\n";


  my $xref_count;
  my $direct_count;

  while ($_ = $file_io->getline()) {

    my $xref;

    my $entry = $_;
    chomp $entry;

    next if (!$entry);

    my ($accession) = $entry =~ /\n#=GF\sAC\s+(\w+)/;
    my ($label) = $entry =~ /\n#=GF\sID\s+([^\n]+)/;
    my ($description) = $entry =~ /\n#=GF\sDE\s+([^\n]+)/;

    if (exists($rfam_transcript_stable_ids{$accession})){

	#add xref
	my $xref_id = $self->add_xref({ acc        => $accession,
				      version    => 0,
				      label      => $label || $accession ,
				      desc       => $description,
				      source_id  => $source_id,
				      species_id => $species_id,
				      info_type  => "DIRECT"} );

	my @transcript_stable_ids = @{$rfam_transcript_stable_ids{$accession}};

	foreach my $stable_id (@transcript_stable_ids){
	    $self->add_direct_xref($xref_id, $stable_id, "Transcript", "");
	    $direct_count++;
	}	

      $xref_count++;

    }

  }

  $file_io->close();

  print "Added $xref_count RFAM xrefs and $direct_count direct xrefs\n" if($verbose);
  if ( !$xref_count ) {
      return 1;    # 1 error
  }

  return 0; # successfull
 

}

1;
