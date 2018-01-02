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

package XrefParser::ncRNA_DBParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use XrefParser::Database;

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my ($type, $my_args) = split(/:/,$file);

  my $user  ="ensro";
  my $host;
  my $port;
  my $dbname;
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
  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }

  my %source;
  $source{"RFAM"}     = $self->get_source_id_for_source_name('RFAM')
    || die "Could not get source_id for RFAM";
  $source{"miRBase"}  = $self->get_source_id_for_source_name('miRBase') 
    || die "Could not get source_id for miRBase";



  my $ccds_db =  XrefParser::Database->new({ host   => $host,
					     port   => $port,
					     user   => $user,
					     dbname => $dbname,
					     pass   => $pass});
  my $dbi2 = $ccds_db->dbi();

  if(!defined($dbi2)){
    return 1;
  }

  my $sql=(<<SQL);
  SELECT transcript_stable_id, transcript_dbid, source, xref_name, xref_primary_id, xref_description
    FROM ncRNA_Xref
      WHERE taxonomy_id = $species_id
SQL

  my ($stable_id, $dbid, $source_name, $label, $acc, $desc);
  my $sth = $dbi2->prepare($sql);
  $sth->execute() or croak( $dbi2->errstr() );
  $sth->bind_columns(\$stable_id, \$dbid, \$source_name, \$label, \$acc, \$desc);
  my $added = 0;
  while ( $sth->fetch() ) {
    my ($description,$junk) = split("[[]Source:",$desc);

    if(!defined($source{$source_name})){
      print STDERR "Could not find source_id for source $source_name for xref $acc\n";
      next;
    }
    my $xref_id = $self->get_xref($acc,$source{$source_name}, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref({ acc        => $acc,
				   label      => $label,
				   desc       => $description,
				   source_id  => $source{$source_name},
				   species_id => $species_id,
				   info_type  => "DIRECT"} );
      $added++;
    }
    my $transcript_id = $dbid;
    if(defined($stable_id) and $stable_id ne ""){
      $transcript_id = $stable_id;
    }
    $self->add_direct_xref($xref_id, $transcript_id, "Transcript", "") if (defined($transcript_id));    

  }
  $sth->finish;

  print "Added $added Xrefs for ncRNAs\n" if($verbose);

  if ($added > 0) {
      return 0;
  } else {
      return 1;
  }
}

1;
