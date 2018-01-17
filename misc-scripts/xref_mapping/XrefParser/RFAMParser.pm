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

package XrefParser::RFAMParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use Bio::EnsEMBL::Registry;


sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $species_name = $ref_arg->{species};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $core_db      = $ref_arg->{dba};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $wget;
  my $user = "ensro";
  my $host;
  my $port = 3306;
  my $dbname;
  my $pass;

  if($file =~ /wget[=][>](\S+?)[,]/){
    $wget = $1;
  }
  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }
  if($file =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }


  #get direct RFAM xrefs from core
  my $registry = "Bio::EnsEMBL::Registry";
  my $dba;

  #get the species name
  my %id2name = $self->species_id2name($dbi);
  if (defined $species_name) { push @{$id2name{$species_id}}, $species_name; }
  if (!defined $id2name{$species_id}) { next; }
  $species_name = $id2name{$species_id}[0];

  if ($host) {
      $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
         '-host'     => $host,
         '-user'     => $user,
         '-pass'     => $pass,
         '-dbname'   => $dbname,
         '-port'     => $port,
         '-species'  => $species_name,
         '-group'    => 'core',
       );
  } elsif (defined $core_db) {
    $dba = $core_db;
  } else {
      $registry->load_registry_from_multiple_dbs( 
      {
        '-host'    => 'mysql-ens-sta-1',
	'-port'    => 4519,
        '-user'    => 'ensro',
      },
      );
      $dba = $registry->get_DBAdaptor($species_name, 'core');
  }

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

  my @lines;
  if (defined $wget) {
    my $ua = LWP::UserAgent->new();
    $ua->timeout(10);
    $ua->env_proxy();
    my $request = HTTP::Request->new(GET => $wget);
    my $response = $ua->request($request);

    if ( !$response->is_success() ) {
      warn($response->status_line);
      return 1;
    }
    @lines = split(/\n\n\n/, $response->decoded_content);
  } else {
    my $file_io = $self->get_filehandle($file);
    if ( !defined $file_io ) {
      print "ERROR: Can't open HGNC file $file\n";
      return 1;
    }
    while (my $line = $file_io->getline()) {
      push(@lines, $line); 
    }
  }

  my @xrefs;
  my $xref_count = 0;
  my $direct_count = 0;

  while (my $entry = shift @lines) {

    my $xref;
    chomp $entry;
    next if (!$entry);

    my ($accession) = $entry =~ /#=GF\sAC\s+(\w+)/ ;
    my ($label) = $entry =~ /\n#=GF\sID\s+([^\n]+)/; 
    my ($description) = $entry =~ /\n#=GF\sDE\s+([^\n]+)/;
    if ($accession) {
      if (exists($rfam_transcript_stable_ids{$accession})){
      #add xref
        my $xref_id = $self->add_xref({ acc      => $accession,
  				      version    => 0,
				      label      => $label || $accession ,
				      desc       => $description,
				      source_id  => $source_id,
				      species_id => $species_id,
                                      dbi        => $dbi,
				      info_type  => "DIRECT"} );

        my @transcript_stable_ids = @{$rfam_transcript_stable_ids{$accession}};
        foreach my $stable_id (@transcript_stable_ids){
           $self->add_direct_xref($xref_id, $stable_id, "Transcript", "", $dbi);
           $direct_count++;
         }	
         $xref_count++;
      }
    }
  }

  print "Added $xref_count RFAM xrefs and $direct_count direct xrefs\n" if($verbose);
  if ( !$xref_count ) {
      return 1;    # 1 error
  }

  return 0; # successfull
 

}

1;
