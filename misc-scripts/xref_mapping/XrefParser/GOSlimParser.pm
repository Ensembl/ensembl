package XrefParser::GOSlimParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";


sub run_script {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id  = shift;
  my $species_id = shift;
  my $verbose    = shift;

  my $user = "ensro";
  my $host;
  my $port = 3306;
  my $dbname;
  my $pass;
  my $tran_name;


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
  if($file =~ /tran_name[=][>](\S+?)[,]/){
    $tran_name = $1;
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
    $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  }

  my $add_dependent_xref_sth = $self->dbi->prepare("INSERT INTO dependent_xref  (master_xref_id,dependent_xref_id, linkage_source_id) VALUES (?,?, $source_id)");


  if(!defined($dbi2)){
    return 1;
  }

  my (%go)   = %{XrefParser::BaseParser->get_valid_codes("GO",$species_id)};


  my $count = 0;

  my $xref_sth = $dbi2->prepare("SELECT t.accession, s.name, s.accession FROM term t, term s, aux_GO_goslim_goa_map ts WHERE ts.term_id = t.term_id and ts.subset_term_id = s.term_id");

  $xref_sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $xref_sth->fetchrow_array() ) {
    my $term_acc     = $row[0];
    my $desc         = $row[1];
    my $subterm_acc  = $row[2];

    if(defined($go{$term_acc})){
      my $xref_id = $self->add_xref($subterm_acc, undef, $subterm_acc, $desc, $source_id, $species_id, "DEPENDENT");
      $add_dependent_xref_sth->execute($go{$term_acc}, $xref_id);
      $count++;
    }

  }
  $xref_sth->finish;
  print "Parsed GOSlim GOA identifiers from $file, added $count dependent_xrefs\n" if($verbose);

  return 0;
}

1;
