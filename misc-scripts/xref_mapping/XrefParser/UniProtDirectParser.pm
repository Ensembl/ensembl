package XrefParser::UniProtDirectParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Uniprot records and assign direct xrefs
# All assumed to be linked to translation

my $verbose;

# --------------------------------------------------------------------------------

sub run {

  my $self = shift;
  my $source_id  = shift;
  my $species_id = shift;
  my $files  = shift;
  my $rel_file   = shift;
  $verbose       = shift;

  my %prefix = (9606 => "ENSP0", 10090 => "ENSMUSP0", 10116 => "ENSRNOP0");

  if(!defined($prefix{$species_id})){
    print "No prefix known for this species $species_id???\n";
    return 1;
  }

  my $filename = @{$files}[0];

  my $file_io = $self->get_filehandle($filename);
  if ( !defined($file_io) ) {
    return 1;
  }

  my $parsed_count = 0;


  my %prot2ensembl;

  my $count = 0;
  while ( defined( my $line = $file_io->getline() ) ) {
    my ($prot, $ens) = split /\s+/,$line;
    if($ens =~ /$prefix{$species_id}/){
      push @{$prot2ensembl{$prot}}, $ens;
   }
  }
  my $dbi = XrefParser::BaseParser->dbi();

  my $sw_source_id =  XrefParser::BaseParser->get_source_id_for_source_name("uniprot/swissprot","sequence_mapped");
  if($sw_source_id < 1){
    die "Could not find source id for uniprot/swissprot ???\n";
  }
  else{
    print "Source_id = $sw_source_id\n";
  }
  my $get_desc_sth = $dbi->prepare("select xref_id, version, label, description from xref where source_id = $sw_source_id and accession = ?");


  my $get_dependents_sth = $dbi->prepare("select dependent_xref_id, linkage_annotation, linkage_source_id  from dependent_xref where master_xref_id = ?");

  my $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref (master_xref_id,dependent_xref_id,linkage_annotation, linkage_source_id) VALUES (?,?,?,?)");


  my $get_aliases_sth =  $dbi->prepare("select synonym from synonym where xref_id = ?");
  my $add_alias_sth   =  $dbi->prepare("INSERT INTO synonym (xref_id, synonym) VALUES (?, ?)");



  my $err_count;
  foreach my $key (keys %prot2ensembl){

    #
    # get the descrptions etc for the uniprot entry
    #
    $get_desc_sth->execute($key);
    my ($old_xref_id, $version, $label, $description);
    $get_desc_sth->bind_columns(\$old_xref_id, \$version, \$label, \$description);
    $get_desc_sth->fetch;
    if(!defined($old_xref_id)){
      print STDERR "Could not find $key in the database\n" if ($err_count <10);
      $err_count++;
      next;
    }
    $count++;

    #
    # get the dependents
    #
    my %linkage_anotation=();
    my %linkage_source_id=();
    my ($dependent_xref_id, $linkage_annotation, $linkage_source_id);
    $get_dependents_sth->execute($old_xref_id);
    $get_dependents_sth->bind_columns(\$dependent_xref_id, \$linkage_annotation, \$linkage_source_id);
    while($get_dependents_sth->fetch){
      $linkage_anotation{$dependent_xref_id} =  $linkage_annotation;
      $linkage_source_id{$dependent_xref_id} =  $linkage_source_id;
    }

#    print $key."\t";
    #
    # Add the new xref
    #

    my $xref_id = XrefParser::BaseParser->add_xref($key, $version, $label, $description, $source_id, $species_id, "DIRECT");


    #
    # Add the synonyms
    #
    my $synonym;
    $get_aliases_sth->execute($old_xref_id);
    $get_aliases_sth->bind_columns(\$synonym);
    while($get_aliases_sth->fetch()){
      $add_alias_sth->execute($xref_id, $synonym) || die "Could not add synonym for $xref_id, $synonym";
    }


    foreach my $trans (@{$prot2ensembl{$key}}){
      #
      #add the direct xref entry
      #

      XrefParser::BaseParser->add_direct_xref( $xref_id, $trans, "Translation", '');
#      print ":".$trans;

      #
      #add the dependents
      #
      foreach my $dep (keys %linkage_anotation){
	$add_dependent_xref_sth->execute($xref_id, $dep, $linkage_anotation{$dep}, $linkage_source_id{$dep});	
      }
    }
  }


  print $count." entrys added\n".$err_count." not found\n";
  return 0;
}


1;
