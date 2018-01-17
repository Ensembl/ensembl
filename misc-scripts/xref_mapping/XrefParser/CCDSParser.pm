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

package XrefParser::CCDSParser;

use strict;
use Carp;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;
# Parse file of CCDS records and assign direct xrefs
# All assumed to be linked to transcripts
# The same CCDS may be linked to more than one transcript, but need to only
# add the xref once, so check if it already exists before adding it.

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $db           = $ref_arg->{dba};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file  as pairs";
  }
  $verbose |=0;

  my $user = "ensro";
  my $host;
  my $port = 3306;
  my $dbname;
  my $pass;

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

  my ($ccds_db, $dbi2);
  if (defined $host) {
    $ccds_db =  XrefParser::Database->new({ host   => $host,
					     port   => $port,
					     user   => $user,
					     dbname => $dbname,
					     pass   => $pass});
    $dbi2 = $ccds_db->dbi();
  } elsif (defined $db) {
    $dbi2 = $db->dbc();
  }
  if(!defined($dbi2)){
    return 1;
  }
  
  my $sql =(<<'SCD');
SELECT t.stable_id, x.dbprimary_acc 
  FROM xref x, object_xref ox, transcript t, external_db e
    WHERE x.xref_id=ox.xref_id AND
          ox.ensembl_object_type = "Transcript" AND
          ox.ensembl_id = t.transcript_id AND
          e.external_db_id = x.external_db_id AND
          e.db_name like "Ens\_%\_transcript"
SCD

  my %seen;

  my $sth = $dbi2->prepare($sql) or die "Could not prepare sql $sql\n";
  $sth->execute() or die "Could not execute $sql\n";
  my $xref_count = 0;
  my $direct_count=0;
  my ($stable_id, $display_label);
  $sth->bind_columns( \$display_label,\$stable_id);
  while ( $sth->fetch ) {

    my ($acc, $version) = split (/\./,$display_label);

    my $xref_id;
    if (!defined($seen{$display_label})) {
      $xref_id = $self->add_xref({ acc        => $acc,
				   version    => $version,
				   label      => $display_label,
				   source_id  => $source_id,
				   species_id => $species_id,
                                   dbi        => $dbi,
				   info_type  => "DIRECT"} );
      $xref_count++;
      $seen{$display_label} = $xref_id;
    }
    else{
      $xref_id = $seen{$display_label};
    }

    $self->add_direct_xref($xref_id, $stable_id, "Transcript", "", $dbi);
    $direct_count++;
  }

  print "Parsed CCDS identifiers from $file, added $xref_count xrefs and $direct_count direct_xrefs\n" if($verbose);

  return 0;
}

1;
