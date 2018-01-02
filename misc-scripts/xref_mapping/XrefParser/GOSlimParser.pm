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

package XrefParser::GOSlimParser;

use strict;
use warnings;
use DBI;
use Carp;
use base qw( XrefParser::BaseParser );
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";
use XrefParser::Database;

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

  my $user = "ensro";
  my $host;
  my $port = "3306";
  my $dbname;
  my $pass;

  if($file =~ /host[=][>](\S+?)[,]/x){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/x){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/x){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/x){
    $pass = $1;
  }
  if($file =~ /user[=][>](\S+?)[,]/x){
    $user = $1;
  }


  my $dbi2;
  if(!defined($dbname)){
    $reg->load_registry_from_db(
                                -host => $host,
                                -user => $user,
			        -group => "ontology");
    my $dbc = $reg->get_adaptor("multi","ontology","GOTerm");
    $dbi2 = $dbc->dbc;
  }
  else{
    my $db =  XrefParser::Database->new({ host   => $host,
					  port   => $port,
					  user   => $user,
					  dbname => $dbname,
					  pass   => $pass});
    $dbi2 = $db->dbi();
  }

  my $add_dependent_xref_sth = $self->dbi->prepare("INSERT INTO dependent_xref  (master_xref_id,dependent_xref_id, linkage_annotation, linkage_source_id) VALUES (?, ?, 'IEA', $source_id)");


  if(!defined($dbi2)){
    print STDERR "Could not connect to ontology database\n";
    return 1;
  }

  my (%go)   = %{$self->get_valid_codes("GO",$species_id)};


  my $count = 0;

  my $xref_sth = $dbi2->prepare("SELECT t.accession, s.name, s.accession FROM term t, term s, aux_GO_goslim_generic_map ts WHERE ts.term_id = t.term_id and ts.subset_term_id = s.term_id");

  $xref_sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $xref_sth->fetchrow_array() ) {
    my $term_acc     = $row[0];
    my $desc         = $row[1];
    my $subterm_acc  = $row[2];

    if(defined($go{$term_acc})){
      foreach my $go_xref_id (@{$go{$term_acc}}) {
	my $xref_id = $self->add_xref({ acc        => $subterm_acc,
					label      => $subterm_acc,
					desc       => $desc,
					source_id  => $source_id,
					species_id => $species_id,
					info_type  => "DEPENDENT"} );
	$add_dependent_xref_sth->execute($go_xref_id, $xref_id);
	$count++;
      }
    }

  }
  $xref_sth->finish;
  print "Parsed GOSlim Generic identifiers from $file, added $count dependent_xrefs\n" if($verbose);

  return 0;
}

1;
