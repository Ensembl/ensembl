# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#
# ID   CYC_PIG                 Reviewed;         104 AA.  Swissprot
# ID   Q3ASY8_CHLCH            Unreviewed;     36805 AA.  SPTrEMBL



package XrefParser::UniProtParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );


my $verbose;

sub run {

  my $self = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  $verbose       = shift;

  my $file = @{$files}[0];

  my $species_name;

  my ( $sp_source_id, $sptr_source_id, $sp_release, $sptr_release );

  if(!defined($species_id)){
    ($species_id, $species_name) = $self->get_species($file);
  }

  $sp_source_id =
    $self->get_source_id_for_source_name('Uniprot/SWISSPROT',"sequence_mapped");
  $sptr_source_id =
    $self->get_source_id_for_source_name('Uniprot/SPTREMBL');

  print "SwissProt source id for $file: $sp_source_id\n" if ($verbose);
  print "SpTREMBL source id for $file: $sptr_source_id\n" if ($verbose);
 

  my @xrefs =
    $self->create_xrefs( $sp_source_id, $sptr_source_id, $species_id,
      $file );

  if ( !@xrefs ) {
      return 1;    # 1 error
  }

  # delete previous if running directly rather than via BaseParser
  if (!defined(caller(1))) {
    print "Deleting previous xrefs for these sources\n" if($verbose);
    $self->delete_by_source(\@xrefs);
  }

  # upload
  if(!defined($self->upload_xref_object_graphs(@xrefs))){
    return 1; 
  }


    if ( defined $release_file ) {
        # These two lines are duplicated from the create_xrefs() method
        # below...
        my $sp_pred_source_id =
          $self->get_source_id_for_source_name(
            'Uniprot/SWISSPROT_predicted');
        my $sptr_pred_source_id =
          $self->get_source_id_for_source_name(
            'Uniprot/SPTREMBL_predicted');

        # Parse Swiss-Prot and SpTrEMBL release info from
        # $release_file.
        my $release_io = $self->get_filehandle($release_file);
        while ( defined( my $line = $release_io->getline() ) ) {
            if ( $line =~ m#(UniProtKB/Swiss-Prot Release .*)# ) {
                $sp_release = $1;
                print "Swiss-Prot release is '$sp_release'\n" if($verbose);
            } elsif ( $line =~ m#(UniProtKB/TrEMBL Release .*)# ) {
                $sptr_release = $1;
                print "SpTrEMBL release is '$sptr_release'\n" if($verbose);
            }
        }
        $release_io->close();

        # Set releases
        $self->set_release( $sp_source_id,        $sp_release );
        $self->set_release( $sptr_source_id,      $sptr_release );
        $self->set_release( $sp_pred_source_id,   $sp_release );
        $self->set_release( $sptr_pred_source_id, $sptr_release );
    }


  return 0; # successfull
}

# --------------------------------------------------------------------------------
# Get species (id and name) from file
# For UniProt files the filename is the taxonomy ID

sub get_species {
  my $self = shift;
  my ($file) = @_;

  my ($taxonomy_id, $extension) = split(/\./, basename($file));

  my $sth = $self->dbi()->prepare("SELECT species_id,name FROM species WHERE taxonomy_id=?");
  $sth->execute($taxonomy_id);
  my ($species_id, $species_name);
  while(my @row = $sth->fetchrow_array()) {
    $species_id = $row[0];
    $species_name = $row[1];
  }
  $sth->finish;

  if (defined $species_name) {

    print "Taxonomy ID " . $taxonomy_id . " corresponds to species ID " . $species_id . " name " . $species_name . "\n" if($verbose);

  } else {

    print STDERR "Cannot find species corresponding to taxonomy ID " . $species_id . " - check species table\n";
    exit(1);

  }

  return ($species_id, $species_name);

}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects

sub create_xrefs {
  my $self = shift;

  my ( $sp_source_id, $sptr_source_id, $species_id, $file ) = @_;

  my $num_sp = 0;
  my $num_sptr = 0;
  my $num_sp_pred = 0;
  my $num_sptr_pred = 0;

  my %dependent_sources = $self->get_dependent_xref_sources(); # name-id hash
  my %GeneNameSynonym;

  if(defined($dependent_sources{'MGI'})){
    $dependent_sources{'MGI'} = XrefParser::BaseParser->get_source_id_for_source_name("MGI","uniprot");
  }	


  # Get predicted equivalents of various sources used here
    my $sp_pred_source_id =
      $self->get_source_id_for_source_name(
        'Uniprot/SWISSPROT_predicted');

    my $sptr_pred_source_id =
      $self->get_source_id_for_source_name(
        'Uniprot/SPTREMBL_predicted');

#  my $go_source_id = $self->get_source_id_for_source_name('GO');
  my $embl_pred_source_id = $dependent_sources{'EMBL_predicted'};
  my $protein_id_pred_source_id = $dependent_sources{'protein_id_predicted'};
  print "Predicted SwissProt source id for $file: $sp_pred_source_id\n" if($verbose);
  print "Predicted SpTREMBL source id for $file: $sptr_pred_source_id\n" if($verbose);
  print "Predicted EMBL source id for $file: $embl_pred_source_id\n" if($verbose);
  print "Predicted protein_id source id for $file: $protein_id_pred_source_id\n" if($verbose);
#  print "GO source id for $file: $go_source_id\n";

    my (%genemap) =
      %{ $self->get_valid_codes( "mim_gene", $species_id ) };
    my (%morbidmap) =
      %{ $self->get_valid_codes( "mim_morbid", $species_id ) };

    my $uniprot_io = $self->get_filehandle($file);
    if ( !defined $uniprot_io ) { return }

  my @xrefs;

  local $/ = "//\n";

  # Create a hash of all valid taxon_ids for this species
  my %species2tax = $self->species_id2taxonomy();
  my @tax_ids = @{$species2tax{$species_id}};
  my %taxonomy2species_id = map{ $_=>$species_id } @tax_ids;


#
# MGI data needed---------------------------------------------------------
#
  my %mgi_acc_to_desc;
  my %mgi_acc_to_label;
  my %mgi_label_to_desc;
  my %mgi_label_to_acc;

  my $sth = $self->dbi()->prepare("SELECT x.accession, x.label, x.description from xref x, source s where x.source_id = s.source_id and s.name like 'MGI' and s.priority_description like 'descriptions'");
  
  $sth->execute() or croak( $self->dbi()->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $mgi_acc_to_desc{$row[0]}   = $row[2];
    $mgi_acc_to_label{$row[0]}  = $row[1];
    $mgi_label_to_desc{$row[1]} = $row[2];
    $mgi_label_to_acc{$row[1]}  = $row[0];
  }
  $sth->finish;


  #
  # Get the MGI synonyms
  #

  $sth = $self->dbi()->prepare("SELECT sy.synonym, x.accession, x.description from xref x, source s, synonym sy where sy.xref_id = x.xref_id and x.source_id = s.source_id and s.name like 'MGI' and s.priority_description like 'descriptions'");
  
  $sth->execute() or croak( $self->dbi()->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $mgi_label_to_desc{$row[0]} = $row[2];
    $mgi_label_to_acc{$row[0]}  = $row[1];
  }
  $sth->finish;

  #
  # end MGI data needed -------------------------------------------------
  #


  my %dependent_xrefs;

  while ( $_ = $uniprot_io->getline() ) {

    # if an OX line exists, only store the xref if the taxonomy ID that the OX
    # line refers to is in the species table
    # due to some records having more than one tax_id, we need to check them 
    # all and only proceed if one of them matches.
    #OX   NCBI_TaxID=158878, 158879;
    #OX   NCBI_TaxID=103690;

    my ($ox) = $_ =~ /OX\s+[a-zA-Z_]+=([0-9 ,]+);/;
    my @ox = ();
    my $found = 0;

    if ( defined $ox ) {
        @ox = split /\, /, $ox;

        # my %taxonomy2species_id = $self->taxonomy2species_id();

        foreach my $taxon_id_from_file (@ox) {
          if ( exists $taxonomy2species_id{$taxon_id_from_file} ){
            $found = 1;
          }
        }
    }

    next if (!$found); # no taxon_id's match, so skip to next record
    my $xref;

    # set accession (and synonyms if more than one)
    # AC line may have primary accession and possibly several ; separated synonyms
    # May also be more than one AC line
    my ($acc) = $_ =~ /(AC\s+.+)/s; # will match first AC line and everything else

    my @all_lines = split /\n/, $acc;

    # extract ^AC lines only & build list of accessions
    my @accessions;
    foreach my $line (@all_lines) {
      my ($accessions_only) = $line =~ /^AC\s+(.+)/;
      push(@accessions, (split /;\s*/, $accessions_only)) if ($accessions_only);

    }


    if(lc($accessions[0]) eq "unreviewed"){
      print "WARNING: entries with accession of $acc not allowed will be skipped\n";
      next;
    }
    $xref->{INFO_TYPE} = "SEQUENCE_MATCH";
    $xref->{ACCESSION} = $accessions[0];
    for (my $a=1; $a <= $#accessions; $a++) {
      push(@{$xref->{"SYNONYMS"} }, $accessions[$a]);
    }

    # Check for CC (caution) lines containing certain text
    # if this appears then set the source of this and and dependent xrefs to the predicted equivalents
    my $is_predicted = /CC.*EMBL\/GenBank\/DDBJ whole genome shotgun \(WGS\) entry/;

    my ($label, $sp_type) = $_ =~ /ID\s+(\w+)\s+(\w+)/;

    # SwissProt/SPTrEMBL are differentiated by having STANDARD/PRELIMINARY here
    if ($sp_type =~ /^Reviewed/i) {

      $xref->{SOURCE_ID} = $sp_source_id;
      if ($is_predicted) {
	$xref->{SOURCE_ID} = $sp_pred_source_id;
	$num_sp_pred++;
      } else {
	$xref->{SOURCE_ID} = $sp_source_id;
	$num_sp++;
      }
    } elsif ($sp_type =~ /Unreviewed/i) {

      if ($is_predicted) {
	$xref->{SOURCE_ID} = $sptr_pred_source_id;
	$num_sptr_pred++;
      } else {
	$xref->{SOURCE_ID} = $sptr_source_id;
	$num_sptr++;
      }

    } else {

      next; # ignore if it's neither one nor t'other

    }



    # some straightforward fields
    $xref->{LABEL} = $label;
    $xref->{SPECIES_ID} = $species_id;
    $xref->{SEQUENCE_TYPE} = 'peptide';
    $xref->{STATUS} = 'experimental';

    # May have multi-line descriptions
    my ($description_and_rest) = $_ =~ /(DE\s+.*)/s;
    @all_lines = split /\n/, $description_and_rest;

    # extract ^DE lines only & build cumulative description string
    my $description = " ";
    my $name        = "";
    my $sub_description = "";

    foreach my $line (@all_lines) {

      next if(!($line =~ /^DE/));

      # get the data
      if($line =~ /^DE   RecName: Full=(.*);/){
        $name .= $1;
      }
      elsif($line =~ /RecName: Full=(.*);/){
        $description .= $1;
      }
      elsif($line =~ /SubName: Full=(.*);/){
        $name .= $1;
      }


      $description =~ s/^\s*//g;
      $description =~ s/\s*$//g;

      
      my $desc = $name.$description;
      if(!length($desc)){
	$desc = $sub_description;
      }
      $desc =~ s/\(\s*EC\s*\S*\)//g;
      $xref->{DESCRIPTION} = $desc;
    }

    # extract sequence
    my ($seq) = $_ =~ /SQ\s+(.+)/s; # /s allows . to match newline
      my @seq_lines = split /\n/, $seq;
    my $parsed_seq = "";
    foreach my $x (@seq_lines) {
      $parsed_seq .= $x;
    }
    $parsed_seq =~ s/\/\///g;   # remove trailing end-of-record character
    $parsed_seq =~ s/\s//g;     # remove whitespace
    $parsed_seq =~ s/^.*;//g;   # remove everything before last ;

    $xref->{SEQUENCE} = $parsed_seq;
    #print "Adding " . $xref->{ACCESSION} . " " . $xref->{LABEL} ."\n";

    
    my ($gns) = $_ =~ /(GN\s+Name.+)/; # /s allows . to match newline
    my @gn_lines = ();
    if ( defined $gns ) { @gn_lines = split /\n/, $gns }
  
    foreach my $gn (@gn_lines){
      my $gene_name = undef;
      my %depe;
      
      if($gn =~ /Name=(\S+);/){
	$depe{LABEL} = uc($1);
	$depe{ACCESSION} = $self->get_name($xref->{ACCESSION},$depe{LABEL});
	$gene_name = $depe{ACCESSION};

	$depe{SOURCE_NAME} = "Uniprot_genename";
	$depe{SOURCE_ID} = $dependent_sources{"Uniprot_genename"};
	$depe{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
	push @{$xref->{DEPENDENT_XREFS}}, \%depe;
	$dependent_xrefs{"Uniprot_genename"}++;
	my @syn;
	if($gn =~ /Synonyms=([^;]+);/){
	  my $syn = $1;
 	  $syn =~ s/\s+//g;
	  @syn= split(/,/,$syn);
          push (@{$depe{"SYNONYMS"}}, @syn);
	}
      }
    }	

    my ($deps) = $_ =~ /(DR\s+.+)/s; # /s allows . to match newline
    my @dep_lines = ();
    if ( defined $deps ) { @dep_lines = split /\n/, $deps }

    # dependent xrefs - only store those that are from sources listed in the source table
    my ($deps) = $_ =~ /(DR\s+.+)/s; # /s allows . to match newline

    my @dep_lines = ();
    if ( defined $deps ) { @dep_lines = split /\n/, $deps }

    my %seen=();  # per record basis

    foreach my $dep (@dep_lines) {
      #both GO and UniGene have the own sources so ignore those in the uniprot files
      #as the uniprot data should be older
      if($dep =~ /GO/ || $dep =~ /UniGene/){
	next;
      }
      if ($dep =~ /^DR\s+(.+)/) {
	my ($source, $acc, @extra) = split /;\s*/, $1;
	if($source =~ "RGD"){  #using RGD file now instead.
	  next;
	}
	if($source =~ "IPI"){  
	  next;
	}
	if($source =~ "UCSC"){  
	  next;
	}
	if($source =~ "SGD"){  
	  next;
	}
	if($source =~ "HGNC"){
	  next;
	}
	if (exists $dependent_sources{$source} ) {
	  # create dependent xref structure & store it
	  my %dep;
          $dep{SOURCE_NAME} = $source;
          $dep{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
          $dep{SOURCE_ID} = $dependent_sources{$source};

	  if($source =~ /HGNC/){
	    $acc =~ s/HGNC://;
	    $extra[0] =~ s/[.]//;
	    $dep{LABEL} = $extra[0];
	  }
	  $dep{ACCESSION} = $acc;

	  if($source =~ /MGI/){
	    $extra[0] =~ s/[.]$//;
            if($extra[0] =~ /ENSMUSG/ or $extra[0] =~ /OTTMUSG/ ){
               next;  # no extra info gained and it could now link to different MGI
            }
	    $dep{LABEL} = $extra[0];
	    if(defined($mgi_acc_to_label{$acc})){
	      $dep{LABEL} = $mgi_acc_to_label{$acc};
	    }
	    if(defined($mgi_acc_to_desc{$acc})){
	      $dep{DESCRIPTION} = $mgi_acc_to_desc{$acc};
	    }
	    elsif(defined($mgi_label_to_desc{$dep{LABEL}})){ # old mgi number ?? use label 
              $dep{DESCRIPTION} = $mgi_label_to_desc{$dep{LABEL}};
              $dep{ACCESSION}   = $mgi_label_to_acc{$dep{LABEL}};
	    }
            else{
               print "Not found $acc, ".$extra[0]."\n" if($verbose);
            }
	  }

	  if($dep =~ /MIM/){
	    $dep{ACCESSION} = $acc;
	    if(defined($morbidmap{$acc}) and $extra[0] eq "phenotype."){
	      $dep{SOURCE_NAME} = "MIM_MORBID";
	      $dep{SOURCE_ID} = $dependent_sources{"MIM_MORBID"};
	    }
	    elsif(defined($genemap{$acc}) and $extra[0] eq "gene."){
	      $dep{SOURCE_NAME} = "MIM_GENE";
	      $dep{SOURCE_ID} = $dependent_sources{"MIM_GENE"};
	    }
	    elsif($extra[0] eq "gene+phenotype."){
	      $dep{SOURCE_NAME} = "MIM_MORBID";
	      $dep{SOURCE_ID} = $dependent_sources{"MIM_MORBID"};
	      if(defined($morbidmap{$acc})){
		$dependent_xrefs{ $dep{SOURCE_NAME} }++; # get count of depenent xrefs.
		push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
	      }
	      my %dep2;
	      $dep2{ACCESSION} = $acc;
	      $dep2{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
	      $dep2{SOURCE_NAME} = "MIM_GENE";
	      $dep2{SOURCE_ID} = $dependent_sources{"MIM_GENE"};	      
	      if(defined($genemap{$acc})){
		$dependent_xrefs{ $dep2{SOURCE_NAME} }++; # get count of depenent xrefs.
		push @{$xref->{DEPENDENT_XREFS}}, \%dep2; # array of hashrefs
	      }
	      next;
	    }
	    else{
#	      print "missed $dep\n";
	      next;
	    }
	  }
	  if ($source eq "EMBL" && $is_predicted) {
	    $dep{SOURCE_ID} = $embl_pred_source_id
	  };

#	  $dep{ACCESSION} = $acc;
	  $dependent_xrefs{ $dep{SOURCE_NAME} }++; # get count of depenent xrefs.
	  if(!defined($seen{$dep{SOURCE_NAME}.":".$dep{ACCESSION}})){
	    push @{$xref->{DEPENDENT_XREFS}}, \%dep; # array of hashrefs
	    $seen{$dep{SOURCE_NAME}.":".$dep{ACCESSION}} =1;
	  }
	  if($dep =~ /EMBL/){
	    my ($protein_id) = $extra[0];
	    if($protein_id ne "-" and !defined($seen{$source.":".$protein_id})){
	      my %dep2;
	      $dep2{SOURCE_NAME} = $source;
	      $dep2{SOURCE_ID} = $dependent_sources{"protein_id"};
	      if ($is_predicted) {
		$dep2{SOURCE_ID} = $protein_id_pred_source_id
	      };
	      $dep2{LINKAGE_SOURCE_ID} = $xref->{SOURCE_ID};
	      # store accession unversioned
	      $dep2{LABEL} = $protein_id;
	      my ($prot_acc, $prot_version) = $protein_id =~ /([^.]+)\.([^.]+)/;
	      $dep2{ACCESSION} = $prot_acc;
	      $dep2{VERSION} = $prot_acc;
	      $dependent_xrefs{ $dep2{SOURCE_NAME} }++; # get count of dependent xrefs.
	      $seen{$source.":".$protein_id} = 1;
	      push @{$xref->{DEPENDENT_XREFS}}, \%dep2; # array of hashrefs
	    }
	  }
	}
      }
    }

    push @xrefs, $xref;

  }

  $uniprot_io->close();

  print "Read $num_sp SwissProt xrefs and $num_sptr SPTrEMBL xrefs from $file\n" if($verbose);
  print "Found $num_sp_pred predicted SwissProt xrefs and $num_sptr_pred predicted SPTrEMBL xrefs\n" if (($num_sp_pred > 0 || $num_sptr_pred > 0) and $verbose);




#  print "$kount gene anmes added\n";

  print "Added the following dependent xrefs:-\n" if($verbose);
  foreach my $key (keys %dependent_xrefs){
    print $key."\t".$dependent_xrefs{$key}."\n" if($verbose);
  }


  return \@xrefs;

  #TODO - currently include records from other species - filter on OX line??
}

sub get_name {
  my $self = shift;
  my $acc  = shift;
  my $label = shift;

  return $acc;
}
1;
