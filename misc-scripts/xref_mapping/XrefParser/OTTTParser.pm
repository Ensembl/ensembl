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

package XrefParser::OTTTParser;

use strict;
use warnings;
use Carp;
use File::Basename;
use XrefParser::Database;

use base qw( XrefParser::BaseParser );

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my ($type, $my_args) = split(/:/,$file);

  my $user = "ensro";
  my $host ="ens-staging";
  my $port = "3306";
  my $dbname = "homo_sapiens_vega_51_36m";
  my $pass;

  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($my_args =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($my_args =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($my_args =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }


  my $sql =(<<'SQL');
SELECT t.stable_id, x.display_label 
  FROM xref x, object_xref ox , transcript t, external_db e 
    WHERE e.external_db_id = x.external_db_id AND
          x.xref_id = ox.xref_id AND
          t.transcript_id = ox.ensembl_id AND
          e.db_name like ?
SQL

  my %ott_to_enst;

  my $db =  XrefParser::Database->new({ host   => $host,
					port   => $port,
					user   => $user,
					dbname => $dbname,
					pass   => $pass});

  my $dbi2 = $db->dbi();

  if(!defined($dbi2)){
    return 1;
  }

  my $sth = $dbi2->prepare($sql);   # funny number instead of stable id ?????
  $sth->execute("ENST_CDS") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[0]} = $row[1];
  }
  $sth->finish;

  $sth = $dbi2->prepare($sql);
  $sth->execute("ENST_ident") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[0]} = $row[1];
  }
  $sth->finish;

  my $xref_count = 0;
  foreach my $ott (keys %ott_to_enst){
    my $xref_id = $self->add_xref({ acc        => $ott,
				    label      => $ott,
				    source_id  => $source_id,
				    species_id => $species_id,
				    info_type  => "DIRECT"} );
    $xref_count++;
    $self->add_direct_xref($xref_id, $ott_to_enst{$ott}, "transcript", "");
  }
  return 0;
}

1;

