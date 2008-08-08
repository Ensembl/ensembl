package XrefParser::RefSeq_CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Refseq records and assign direct xrefs

sub run {

  my ($self, $source_id, $species_id, $file) = @_;

  my $refseq_io = $self->get_filehandle($file);

  my $dna_pred = XrefParser::BaseParser->get_source_id_for_source_name("RefSeq_dna_predicted");

  if ( !defined $refseq_io ) {
    print "Could not open $file\n";
    return 1;
  }

  # becouse the direct mapping have no descriptions etc
  # we have to steal these from the previous Refseq parser.

  my %label;
  my %version;
  my %description;

  my $dbi = $self->dbi();  
  my $sql = "select xref.accession, xref.label, xref.version,  xref.description from xref, source where xref.source_id = source.source_id and source.name = 'RefSeq_dna'";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver, $desc);
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc;
  }
  $sth->finish;
 
  $sql = 'select x.accession, x.xref_id, d.ensembl_stable_id, d.type 
            from xref x, direct_xref d, source s 
             where s.source_id = x.source_id and 
                   x.xref_id = d.general_xref_id and s.name like "CCDS"'; 
 
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($access, $old_xref_id, $stable_id, $type);
  $sth->bind_columns(\$access, \$old_xref_id, \$stable_id, \$type);
  my %ensembl_stable_id;
  my %ensembl_type;
  my %old_xref;
  while (my @row = $sth->fetchrow_array()) {
      $ensembl_stable_id{$access} = $stable_id;
      $ensembl_type{$access} = $type;
      $old_xref{$access} = $old_xref_id; 
  }
  $sth->finish;
  
 

  my $line_count = 0;
  my $xref_count = 0;
  my %seen;
  my %old_to_new;

  $refseq_io->getline();    # header

  while ( $_ = $refseq_io->getline() ) {
      chomp;
      my ($ccds,$refseq) = split;
    
      $line_count++;
      if(!defined($seen{$refseq})){
	  $seen{$refseq} = 1;
	  my $key = "CCDS".$ccds;
	  if(defined($ensembl_stable_id{$key})){
	    my $new_source_id = $source_id;
	    if($refseq =~ /^XM/){
	      $new_source_id = $dna_pred;
	    }
	    my $xref_id = $self->add_xref($refseq, $version{$refseq} , $label{$refseq}||$refseq , 
					    $description{$refseq}, $new_source_id, $species_id);
	    $self->add_direct_xref($xref_id, $ensembl_stable_id{$key}, $ensembl_type{$key}, "");
	    $old_to_new{$old_xref{$refseq}} = $xref_id;
	    $xref_count++;
	  }
      }
  }
  
#for each one seen get all its dependent xrefs and load them fro the new one too;

  my $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
  my $get_dependent_xref_sth = $dbi->prepare("SELECT dependent_xref_id, linkage_annotation "
					    .  "FROM  dependent_xref where master_xref_id = ?");

  foreach my $old_xref (keys %old_to_new){
      my $linkage;
      my $dependent_id;
      $get_dependent_xref_sth->execute($old_xref);
      $get_dependent_xref_sth->bind_columns(\$dependent_id, \$linkage);
      while(my @row = $get_dependent_xref_sth->fetchrow_array()){
	  $add_dependent_xref_sth->execute($old_to_new{$old_xref}, $dependent_id, $linkage, $source_id); 
      }   
  }

  $refseq_io->close();

  print "Parsed $line_count RefSeq_dna identifiers from $file, added $xref_count xrefs and $xref_count direct_xrefs  from $line_count lines.\n";


  return 0;

}

1;
