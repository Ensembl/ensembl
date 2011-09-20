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
  my $host = "ens-staging1";
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

  my $add_dependent_xref_sth = $self->dbi->prepare("INSERT INTO dependent_xref  (master_xref_id,dependent_xref_id, linkage_source_id) VALUES (?,?, $source_id)");


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
	my $xref_id = $self->add_xref($subterm_acc, undef, $subterm_acc, $desc, $source_id, $species_id, "DEPENDENT");
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
