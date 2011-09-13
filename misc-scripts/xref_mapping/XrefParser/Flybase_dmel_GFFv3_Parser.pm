# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#
# This module will read in the fly gff text file and make xrefs from the information in the file.
# First of all, it read knows what all the gene, transcript and translation types are, found in column 3 of the gff file:
# Gene = gene
# Transcript = mRNA ncRNA snRNA tRNA rRNA pseudogene snoRNA miRNA
# Translation = protein (could include CDS here but haven't ?) 
# 
# ID=FBgn => flybase_gene_id
# ID=FBtr => flybase_transcript_id
# ID=FBpp => flybase_polypeptide_id
# Name=CG0123 => FlyBaseName_gene
# Name=CG0123-RA => FlyBaseName_transcript
# Name=CG0123-PA => FlyBaseName_translations
# Dbxref=FlyBase:FBan => flybase_annotation_id
# Dbxref=FlyBase_Annotation_IDs:CG0123 => gadfly_gene_cgid
# Dbxref=FlyBase_Annotation_IDs:CG0123-RA => gadfly_transcript_cgid
# Dbxref=FlyBase_Annotation_IDs:CG0123-PA => gadfly_translation_id
# Alias= => flybase_synonym
#
# For each line in the gff file for the above list of genes, transcript and translations, the following happens:
# The unique_id is read in from ID= (FBgn, FBtr, FBpp). This is the direct xref for all xrefs of this entry.
# An xref is made for the entry, using the ID as the xref's accession. Synonyms from Alias= are added to this xref.
# The Name (Name=) is read in and added as an xref. Synonyms from Alias= are added to this xref.
# All entries from Dbxref= are added in as xrefs for the entry; they have no synonyms.
  
#2L gene [...]  ID=CG11023;Dbxref=FlyBase:FBan0011023,FlyBase:FBgn0031208;gbunit=AE003590;synonym=CG11023
#2L mRNA [...]  ID=CG11023-RA;Dbxref=FlyBase:FBtr008,FlyBase:FBgn003;dbxref_2nd=Gadfly:CG11023-RA;synonym=CG23-RA
#3R     FlyBase gene    8084471 8128509 .       +       .
#ID=FBgn0003651;Name=svp;Alias=FBgn0011337,FBgn0011492,FBgn0011510,FBgn0038010,FBgn0063263;Ontology_term=SO:0000010,SO:0000087,GO:0004872,GO:0007270,GO:0042331,GO:0005515,GO:0007419,GO:0007503,GO:0045449,GO:0004879,GO:0003700,GO:0005634,GO:0007465,GO:0007462,GO:0007464,GO:0007510,GO:0005737,GO:0007507,GO:0007417,GO:0001700,GO:0006357,GO:0007165,GO:0043565,GO:0003707,GO:0008270,GO:0048749,GO:0001752;Dbxref=FlyBase:FBan0011502,FlyBase_Annotation_IDs:CG11502,INTERPRO:IPR013088,GB:AC007724,GB:AE003695,GB_protein:AAF54773,GB_protein:AAN13541,GB_protein:AAF54774,GB:AI108883,GB:AI402121,GB:AY075272,GB_protein:AAL68139,GB:AY119490,GB_protein:AAM50144,GB:AY129452,GB_protein:AAM76194,GB:BG633933,GB:BI167911,GB:CZ468719,GB:CZ472606,GB:CZ475640,GB:CZ475641,GB:CZ477001,GB:CZ482253,GB:CZ485541,GB:CZ485543,GB:G00472,GB:M28863,GB_protein:AAA62770,GB:M28864,GB_protein:AAA03014,UniProt/Swiss-Prot:P16375,UniProt/Swiss-Prot:P16376,UniProt/TrEMBL:Q8MRP3,INTERPRO:IPR000536,INTERPRO:IPR001628,INTERPRO:IPR001723,INTERPRO:IPR003068,INTERPRO:IPR008946,INTERPRO:IPR013629,dedb:9161,flygrid:66603,hybrigenics:521960,if:/newgene/sevenup.htm,orthologs:ensAG:ENSANGG00000002454,orthologs:ensAM:ENSAPMG00000000116,orthologs:ensCF:ENSCAFG00000008076,orthologs:ensDM:CG12744,orthologs:ensDR:ENSDARG00000017168,orthologs:ensFR:SINFRUG00000127451,orthologs:ensGG:ENSGALG00000007000,orthologs:ensHS:ENSG00000185551,orthologs:ensMM:ENSMUSG00000030551,orthologs:ensPT:ENSPTRG00000007484,orthologs:ensRN:ENSRNOG00000010308,orthologs:ensTN:GSTENG00006911001,orthologs:modCB:WBGene00030075;cyto_range=87B4-87B5;gbunit=AE014297;


package XrefParser::Flybase_dmel_GFFv3_Parser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;
use Bio::EnsEMBL::Utils::Exception;

use base qw( XrefParser::BaseParser );

my %cache_source =();

my $verbose;

sub new {
  my $proto = shift;
  my $self = $proto->SUPER::new(@_);

  $self->external_source_db_name('flybase_gff');

  #  my @gff_obj =qw( CDS exon gene mRNA ncRNA pseudogene rRNA snRNA snoRNA tRNA );
  # this array may need to change between releases so check that it's updated
  my @gff_obj =qw( gene mRNA ncRNA snRNA tRNA rRNA pseudogene snoRNA miRNA);
  $self->gff_object_types(\@gff_obj);

  #
  # hard-coded field separators out of gff
  #

  $self->gff_name("Name=");
  $self->gff_ontology("Ontology_term=");
  $self->gff_synonym("Alias=");
  $self->gff_dbxref("Dbxref=");

  #
  # hard-coded source-names for different objects out of ./sql/populate_metadata.sql
  #
  # For Alias
  $self->source_name_synonym('flybase_synonym'); # source for any Alias
  # For Name
  $self->source_name_name_prefix('FlyBaseName_'); # source for any Name
  # For Dbxref
  $self->source_name_fbgn('flybase_gene_id');         # source-name for ID=FBgn
  $self->source_name_fbtr('flybase_transcript_id');   # source-name for ID=FBtr
  $self->source_name_fbpp('flybase_polypeptide_id');  # source-name for ID=FBpp
  $self->source_name_fban('flybase_annotation_id');   # source-name for ID=FBan
  $self->source_name_gadfly_gene('gadfly_gene_cgid');                # For Dbxref=FlyBase_Annotation_IDs
  $self->source_name_gadfly_transcript('gadfly_transcript_cgid');    # For Dbxref=FlyBase_Annotation_IDs
  $self->source_name_gadfly_translation('gadfly_translation_cgid');  # For Dbxref=FlyBase_Annotation_IDs
  $self->source_name_affymetrix('AFFY_DrosGenome1');    # For Dbxref=Affymetrix
  $self->source_name_dgrc1('DGRC-1');   # For Dbxref=DGRC-1 
  $self->source_name_dgrc2('DGRC-2');   # For Dbxref=DGRC-2
  $self->source_name_drsc('DRSC');      # For Dbxref=DRSC
  $self->source_name_epd('EPD');        # For Dbxref=EPD
  $self->source_name_flyreg('FlyReg');  # For Dbxref=FlyReg
  $self->source_name_gb('EMBL');        # For Dbxref=GB
  $self->source_name_gbprotein('protein_id');   # For Dbxref=GB_protein
  $self->source_name_gcr('GPCR');       # For Dbxref=GCR
  $self->source_name_gi('GI');          # For Dbxref=GI
  $self->source_name_go('GO');          # For Dbxref=GO
  $self->source_name_genomeRNAi('GenomeRNAi');  # For Dbxref=GenomeRNAi
  $self->source_name_interpro('Interpro');      # For Dbxref=INTERPRO
  $self->source_name_merops('MEROPS');  # For Dbxref=MEROPS
  $self->source_name_miRBase('miRBase');        # For Dbxref=miRBase
  $self->source_name_mitodrome('MitoDrome');    # For Dbxref=MitoDrome
  $self->source_name_nrl3d('PDB');      # For Dbxref=NRL_3D
  $self->source_name_pdb('PDB');        # For Dbxref=PDB
  $self->source_name_rfam('RFAM');      # For Dbxref=Rfam
  $self->source_name_tf('TransFac');    # For Dbxref=TF
  $self->source_name_uniprotsp('Uniprot/SWISSPROT');    # For Dbxref=UniProt/Swiss-Prot
  $self->source_name_uniprottr('Uniprot/SPTREMBL');     # For Dbxref=UniProt/TrEMBL
  $self->source_name_bdgpinsituexpr('BDGP_insitu_expr');        # For Dbxref=bdgpinsituexpr
  $self->source_name_dedb('DEDb');      # For Dbxref=dedb
  $self->source_name_drosdel('DrosDel');        # For Dbxref=drosdel
  $self->source_name_flygrid('FlyGrid');        # For Dbxref=flygrid
  $self->source_name_hybrigenics('hybrigenics');        # For Dbxref=hybrigenics
  $self->source_name_if('InteractiveFly');      # For Dbxref=if
  $self->source_name_prefix_ensAGgene('Ens_Ag_gene');    # For Dbxref=ensAG
  $self->source_name_prefix_ensAMgene('Ens_Am_gene');    # For Dbxref=ensAM
  $self->source_name_prefix_ensCEgene('Ens_Ce_gene');    # For Dbxref=ensCE
  $self->source_name_prefix_ensCFgene('Ens_Cf_gene');    # For Dbxref=ensCF
  $self->source_name_prefix_ensDMgene('Ens_Dm_gene');    # For Dbxref=ensDM
  $self->source_name_prefix_ensDRgene('Ens_Dr_gene');    # For Dbxref=ensDR
  $self->source_name_prefix_ensFRgene('Ens_Fr_gene');    # For Dbxref=ensFR
  $self->source_name_prefix_ensGGgene('Ens_Gg_gene');    # For Dbxref=ensGG
  $self->source_name_prefix_ensHSgene('Ens_Hs_gene');    # For Dbxref=ensHS
  $self->source_name_prefix_ensMMgene('Ens_Mm_gene');    # For Dbxref=ensMM
  $self->source_name_prefix_ensPTgene('Ens_Pt_gene');    # For Dbxref=ensPT
  $self->source_name_prefix_ensRNgene('Ens_Rn_gene');    # For Dbxref=ensRN
  $self->source_name_prefix_ensTNgene('Ens_Tn_gene');    # For Dbxref=ensTN
  $self->source_name_prefix_modCBgene('modCB_gene');     # For Dbxref=modCB
  $self->source_name_prefix_modCEgene('modCE_gene');     # For Dbxref=modCE
  $self->source_name_prefix_modDDgene('modDD_gene');     # For Dbxref=modDD

  my @gene_types = qw (gene) ;
  my @translation_types = qw (protein);
  # The transcript_types may change from release to release so check that this list is up-to-date
  my @transcript_types = qw (mRNA ncRNA snRNA tRNA rRNA pseudogene snoRNA miRNA);

  $self->gene_types(\@gene_types) ;
  $self->translation_types(\@translation_types) ;
  $self->transcript_types(\@transcript_types) ;

  $self->{'_xrefs'}=[];
  $self->{'_direct_xrefs'}=[];
  $self->{'_synonyms'}={};
  
  return $self;
}


# --------------------------------------------------------------------------------



# large number of calls to SQL should now be speeded up as cached.
sub get_source{
  my ($self, $name) =@_;

  if(!defined($cache_source{$name})){
    $cache_source{$name} = XrefParser::BaseParser->get_source_id_for_source_name($name)
  }

  return $cache_source{$name};
}

sub run {


  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $species_name;

  if(!defined($species_id)){
    ($species_id, $species_name) = get_species($file);
  }
  $self->species_id($species_id) ;

  my $external_source_db_name = $self->external_source_db_name() ;
  my $flybase_source_id = $self->get_source($external_source_db_name);

  if(!$self->create_xrefs($flybase_source_id, $file)){
    return 1;
  }

  my @xrefs = @{$self->xrefs};

  $self->relink_synonyms_to_xrefs();

  my @direct_xrefs = @{ $self->direct_xrefs } ;

  # delete previous if running directly rather than via BaseParser
  if (!defined(caller(1))) {
    print "Deleting previous xrefs for these sources\n" if($verbose);
    XrefParser::BaseParser->delete_by_source(\@xrefs);
  }
  print "... parsed.\n" if($verbose);
  print STDERR "uploading ".scalar(@xrefs)." xrefs's\n" if($verbose);
  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print STDERR "uploading ".scalar(@direct_xrefs)." direct-xrefs's\n" if($verbose);
  XrefParser::BaseParser->upload_direct_xrefs(\@direct_xrefs);

  return 0;
}

sub relink_synonyms_to_xrefs{
  my $self = shift;
  foreach my $x (@{$self->xrefs} ){
    my $src_name = XrefParser::BaseParser->get_source_name_for_source_id($x->{SOURCE_ID});
    if ($src_name =~ m/^FlyBaseName_/ || $src_name =~ m/^flybase_.*_id$/) {
      $x->{SYNONYMS} = $self->get_synonyms($x->{ENSEMBL_STABLE_ID});
    }
  }
}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects
# parse

sub create_xrefs {
  my ($self, $flybase_source_id, $file) = @_;

  print STDERR "starting to parse $file...." if($verbose);

  my $gff_io = $self->get_filehandle($file);

  if ( !defined $gff_io ) {
    print STDERR "ERROR: Can't open the GFF file $file\n";
    return 0;
  }

  while ( $_ = $gff_io->getline() ) {
    chomp;
    my @col = split /\t/;
    if($col[3]){

      # test if line contains information for object wanted (CDS,mRNA,gene,..)
      if ( $self->line_contains_object_to_process( $col[2] ) ){
        # work out if we have a gene, transcript or translation
 	my $type = $self->set_ensembl_object_type($col[2]);
        # the 9th column contains all the attributes
	my @desc = split /\;/,$col[8];
        # the ID= is always the first element of this array
	my $unique_id = shift @desc;
	if(!$unique_id=~m/ID=/){
	  throw("parse-error: There seems to be no Identifier: $unique_id. Suspicous!");
	  # print "parse-error: There seems to be no Identifier: $unique_id. Suspicous!";
	  # return 0;
	}
        # for a gene, this will be FBgn, for a transcript this will be FBtr, etc
	$unique_id =~s/ID=//g;
        $self->make_id_xref($unique_id,$type);
	# set up xref-entry for EVERY single item
	foreach my $item (@desc) {
          $self->set_flybase_synonyms($item,$unique_id);
	  # make all xrefs for type "Name=" in desc-field
          # these are FlyBaseName_gene for genes, FlyBaseName_transcript for transcripts, etc
	  $self->make_name_xref($item,$unique_id,$type);

	  # make all xrefs for type "Name=" in desc-field
	  $self->make_dbxref_xref($item,$unique_id,$type);
	}
      }
	  
    } # we don't want to read the line otherwise
	
  } # while ( $_ = $gff_io->getline() ) {

  $gff_io->close();

  return 1;
}

sub set_ensembl_object_type{
  my ($self,$t) = @_ ; # $t is identifier in gff for object : CDS,mRNA,gene,pseudogene,snRNA,....

  for my $hc (@{ $self->gene_types } ){
    if ($t=~m/^$hc$/){
      return 'gene';
    }
  }
  for my $hc (@{ $self->translation_types } ){
    if ($t=~m/^$hc$/){
      return 'translation';
    }
  }
  for my $hc (@{ $self->transcript_types} ){
    if ($t=~m/^$hc$/){
      return 'transcript';
    }
  }
}


sub make_dbxref_xref{
  my ($self,$item,$unique_id,$type) = @_;
  # item = attribute 
  # unique_id = ID
  # type = gene, transcript, translation
  my ($xref);
  my $tg1 = $self->gff_dbxref ;
  my $tg2 = $self->gff_ontology;

  if ($item=~/$tg1/ || $item=~/$tg2/){ # Dbxref=
    # split the xrefs up into a list
    my $dbx1 = get_fields($item,$tg1);
    my @dbx;
    push @dbx, @{$dbx1} if $dbx1;

    foreach my $dbx (@dbx) {
      my $src_id = undef;
      my $source_type = undef;

      if ($dbx =~m/FlyBase:/){
	$dbx =~s/FlyBase://g;

	if($dbx=~m/FBgn/ and $type eq "gene"){
	  $src_id = $self->get_source($self->source_name_fbgn);
	}elsif ($dbx =~m/FBtr/ and $type eq "transcript"){
	  $src_id = $self->get_source($self->source_name_fbtr);
	}elsif ($dbx =~m/FBpp/ and $type eq "translation"){
	  $src_id = $self->get_source($self->source_name_fbpp);
	}elsif ($dbx =~m/FBan/){
	  $src_id = $self->get_source($self->source_name_fban);
	}
      }elsif($dbx =~m/FlyBase_Annotation_IDs:/){
	$dbx =~s/FlyBase_Annotation_IDs://g;
	if($type eq "gene"){
	  $src_id = $self->get_source($self->source_name_gadfly_gene) ;
	}
	elsif($type eq "translation"){
	  $src_id = $self->get_source($self->source_name_gadfly_translation);
	}
	elsif($type eq "transcript"){
	  $src_id = $self->get_source($self->source_name_gadfly_transcript); 
	}
      } elsif ($dbx =~m/Affymetrix:/) {
        $dbx =~s/Affymetrix://g;
        $src_id = $self->get_source($self->source_name_affymetrix) ;
      } elsif ($dbx =~m/DGRC-1:/) {
        $dbx =~s/DGRC-1://g;
        $src_id = $self->get_source($self->source_name_dgrc1) ;
      } elsif ($dbx =~m/DGRC-2:/) {
        $dbx =~s/DGRC-2://g;
        $src_id = $self->get_source($self->source_name_dgrc2);
      } elsif ($dbx =~m/DRSC:/) {
        $dbx =~s/DRSC://g;
        $src_id = $self->get_source($self->source_name_drsc);
      } elsif ($dbx =~m/EPD:/) {
        $dbx =~s/EPD://g;
        $src_id = $self->get_source($self->source_name_epd);
      } elsif ($dbx =~m/FlyReg:/) {
        $dbx =~s/FlyReg://g;
        $src_id = $self->get_source($self->source_name_flyreg);
      } elsif ($dbx =~m/GB:/) {
        $dbx =~s/GB://g;
        $src_id = $self->get_source($self->source_name_gb);
      } elsif ($dbx =~m/GB_protein:/) {
        $dbx =~s/GB_protein://g;
        $src_id = $self->get_source($self->source_name_gbprotein);
      } elsif ($dbx =~m/GCR:/) {
        $dbx =~s/GCR://g;
        $src_id = $self->get_source($self->source_name_gcr);
      } elsif ($dbx =~m/GI:/) {
        $dbx =~s/GI://g;
        $src_id = $self->get_source($self->source_name_gi);
      } elsif ($dbx =~m/GO:/) {
        # this is an ontology_term
        $dbx =~s/GO://g;
        $src_id = $self->get_source($self->source_name_go);
      } elsif ($dbx =~m/GenomeRNAi:/) {
        $dbx =~s/GenomeRNAi://g;
        $src_id = $self->get_source($self->source_name_genomeRNAi);
      } elsif ($dbx =~m/INTERPRO:/) {
        $dbx =~s/INTERPRO://g;
        $src_id = $self->get_source($self->source_name_interpro);
      } elsif ($dbx =~m/MEROPS:/) {
        $dbx =~s/MEROPS://g;
        $src_id = $self->get_source($self->source_name_merops);
      } elsif ($dbx =~m/MIR:/) {
        $dbx =~s/MIR://g;
        $src_id = $self->get_source($self->source_name_miRBase);
      } elsif ($dbx =~m/MITODROME:/) {
        $dbx =~s/MITODROME://g;
        $src_id = $self->get_source($self->source_name_mitodrome);
      } elsif ($dbx =~m/NRL_3D:/) {
        $dbx =~s/NRL_3D://g;
        $src_id = $self->get_source($self->source_name_nrl3d);
      } elsif ($dbx =~m/PDB:/) {
        $dbx =~s/PDB://g;
        $src_id = $self->get_source($self->source_name_pdb);
      } elsif ($dbx =~m/Rfam:/) {
        $dbx =~s/Rfam://g;
        $src_id = $self->get_source($self->source_name_rfam);
      } elsif ($dbx =~m/SO:/) {
        # do nothing, we don't collect these
      } elsif ($dbx =~m/TF:/) {
        $dbx =~s/TF://g;
        $src_id = $self->get_source($self->source_name_tf);
      } elsif ($dbx =~m/UniProt\/Swiss-Prot:/) {
        $dbx =~s/UniProt\/Swiss-Prot://g;
        $src_id = $self->get_source($self->source_name_uniprotsp);
      } elsif ($dbx =~m/UniProt\/TrEMBL:/) {
        $dbx =~s/UniProt\/TrEMBL://g;
        $src_id = $self->get_source($self->source_name_uniprottr);
      } elsif ($dbx =~m/bdgpinsituexpr:/) {
        $dbx =~s/bdgpinsituexpr://g;
        $src_id = $self->get_source($self->source_name_bdgpinsituexpr);
      } elsif ($dbx =~m/dedb:/) {
        $dbx =~s/dedb://g;
        $src_id = $self->get_source($self->source_name_dedb);
      } elsif ($dbx =~m/drosdel:/) {
        $dbx =~s/drosdel://g;
        $src_id = $self->get_source($self->source_name_drosdel);
      } elsif ($dbx =~m/flygrid:/) {
        $dbx =~s/flygrid://g;
        $src_id = $self->get_source($self->source_name_flygrid);
      } elsif ($dbx =~m/hybrigenics:/) {
        $dbx =~s/hybrigenics://g;
        $src_id = $self->get_source($self->source_name_hybrigenics);
      } elsif ($dbx =~m/if:/) {
        $dbx =~s/if://g;
        $src_id = $self->get_source($self->source_name_if);
      } elsif ($dbx =~m/orthologs:ensAG:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensAGgene);
      } elsif ($dbx =~m/orthologs:ensAM:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensAMgene);
      } elsif ($dbx =~m/orthologs:ensCE:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensCEgene);
      } elsif ($dbx =~m/orthologs:ensCF:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensCFgene);
      } elsif ($dbx =~m/orthologs:ensDM:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensDMgene);
      } elsif ($dbx =~m/orthologs:ensDR:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensDRgene);
      } elsif ($dbx =~m/orthologs:ensFR:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensFRgene);
      } elsif ($dbx =~m/orthologs:ensGG:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensGGgene);
      } elsif ($dbx =~m/orthologs:ensHS:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensHSgene);
      } elsif ($dbx =~m/orthologs:ensMM:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensMMgene);
      } elsif ($dbx =~m/orthologs:ensPT:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensPTgene);
      } elsif ($dbx =~m/orthologs:ensRN:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensRNgene);
      } elsif ($dbx =~m/orthologs:ensTN:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_ensTNgene);
      } elsif ($dbx =~m/orthologs:modCB:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_modCBgene);
      } elsif ($dbx =~m/orthologs:modCE:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_modCEgene);
      } elsif ($dbx =~m/orthologs:modDD:/) {
        $dbx =~s/orthologs://g;
        $src_id = $self->get_source($self->source_name_prefix_modDDgene);
      } else {
        warning("Dbxref type not recognised : $dbx");
      }

      if ($src_id){ # only add xref entry for FBgn FBtr...
	my $xref ;
	$xref->{ACCESSION} = $dbx ;
	$xref->{LABEL} = $dbx;
	$xref->{SOURCE_ID} = $src_id;
	$xref->{SPECIES_ID} = $self->species_id();
	#$xref->{SYNONYMS} = $self->get_synonyms($unique_id);
	$self->add_xref($xref);

	if ($type){
	  my $direct_xref;
	  $direct_xref = $xref ;
	  $direct_xref->{ENSEMBL_STABLE_ID} = $unique_id;
	  $direct_xref->{ENSEMBL_TYPE} = $type;
	  #$direct_xref->{LINKAGE_XREF}=undef;
	  $self->add_direct_xref($direct_xref) if $type ;
	}
      }
    }
    return;
  }
}

sub set_flybase_synonyms {
  my ($self,$item,$unique_id) = @_; 
  my $syn1 = $self->gff_synonym;

  if ($item=~/$syn1/){
    my $s1 = get_fields($item,$syn1);
    my @syns;
    push @syns, @{$s1} if $s1;
    $self->add_synonym($unique_id,\@syns);
    return \@syns;
  }
  return;
}

sub make_id_xref{
  my ($self,$unique_id,$type) = @_;
  my $xref=undef;

  # make an xref
  $xref->{ACCESSION} = $unique_id;
  $xref->{LABEL} = $unique_id;
  $xref->{SPECIES_ID} = $self->species_id();
  $xref->{SYNONYMS} = $self->get_synonyms($unique_id);
  my $type_s = $type;
  if ($type eq "gene") {
    $type_s = $self->source_name_fbgn();
  } elsif ($type eq "transcript") {
    $type_s = $self->source_name_fbtr();
  } elsif ($type eq "translation") {
    $type_s =  $self->source_name_fbpp();
  } else {
    throw ("Type $type not recognised");
  }

  $xref->{SOURCE_ID} =  $self->get_source($type_s);
	$self->add_xref($xref);

  # only allow Name on genes. This is a fix for Biomart really.
  if (defined($xref) and $type){
    my $direct_xref;
    $direct_xref = $xref ;
    $direct_xref->{ENSEMBL_STABLE_ID} = $unique_id;
    $direct_xref->{ENSEMBL_TYPE} = $type;
    $direct_xref->{LINKAGE_TYPE}='bla';
    $direct_xref->{SYNONYMS} = $self->get_synonyms($unique_id);
    $self->add_direct_xref($direct_xref);
  }
  return;
}

sub make_name_xref{
  my ($self,$item,$unique_id,$type) = @_;
  my $xref=undef;
  my $target = $self->gff_name ;
  if($item=~m/$target/){  ##Name=
    #print "having $$gff_gene_name[0]\n" ;
    # remove the Name= bit and split the names on a ','
    my $gff_gene_name = get_fields ( $item, $target ) ;
    throw("there is more than one id for item $item\n") if $$gff_gene_name[1];
    $xref->{ACCESSION} = $$gff_gene_name[0];
    $xref->{LABEL} = $$gff_gene_name[0];
    $xref->{SPECIES_ID} = $self->species_id();
    $xref->{SYNONYMS} = $self->get_synonyms($unique_id);
    my $type_s = $type;
    if($type eq "translation"){
      $type_s = $type."s";
    }
    $xref->{SOURCE_ID} =  $self->get_source($self->source_name_name_prefix().$type_s);
	$self->add_xref($xref);
  }
  # only allow Name on genes. This is a fix for Biomart really.
  if (defined($xref) and $type){
    my $direct_xref;
    $direct_xref = $xref ;
    $direct_xref->{ENSEMBL_STABLE_ID} = $unique_id;
    $direct_xref->{ENSEMBL_TYPE} = $type;
    $direct_xref->{LINKAGE_TYPE}='bla';
    $direct_xref->{SYNONYMS} = $self->get_synonyms($unique_id);
    $self->add_direct_xref($direct_xref);
  }
  return;
}

sub get_fields {
  my ($item,$target) =@_;

  my @entrys;
  if ($item =~m/$target/){
    $item =~s/$target//g;

    # check if there is more than one synonym / dbxref ...
    if ($item =~/,/){
      @entrys = split (/\,/,$item);
    } else{
      push @entrys, $item;
    }
    return \@entrys;

    # if the item does not hold information of specific field
  }else{
    return;
  }
}

sub source_name_name{
  my $self = shift;

  $self->{_source_name_name} = shift if @_ ;
  return $self->{_source_name_name};
}

sub source_name_name_prefix{
  my $self = shift;

  $self->{_source_name_name_prefix} = shift if @_ ;
  return $self->{_source_name_name_prefix};
}


sub source_name_synonym{
  my $self = shift;

  $self->{_source_name_synonym} = shift if @_ ;
  return $self->{_source_name_synonym};
}


sub source_name_fbgn{
  my $self = shift;

  $self->{_source_name_gene} = shift if @_ ;
  return $self->{_source_name_gene};
}


sub source_name_gadfly_gene{
  my $self = shift;

  $self->{_source_name_gadfly_gene} = shift if @_ ;
  return $self->{_source_name_gadfly_gene};
}

sub source_name_gadfly_transcript{
  my $self = shift;

  $self->{_source_name_gadfly_transcript} = shift if @_ ;
  return $self->{_source_name_gadfly_transcript};
}
sub source_name_gadfly_translation{
  my $self = shift;

  $self->{_source_name_gadfly_translation} = shift if @_ ;
  return $self->{_source_name_gadfly_translation};
}


sub source_name_fbtr{
  my $self = shift;

  $self->{_source_name_transcript} = shift if @_ ;
  return   $self->{_source_name_transcript}  ;
}

sub source_name_fbpp{
  my $self = shift;

  $self->{_source_name_fbpp} = shift if @_ ;
  return $self->{_source_name_fbpp};
}

sub source_name_fban{
  my $self = shift;

  $self->{_sn_fban} = shift if @_ ;
  return $self->{_sn_fban};
}

sub source_name_affymetrix {
  my $self = shift;
  $self->{_sn_affymetrix} = shift if @_ ;
  return $self->{_sn_affymetrix};
}

sub source_name_dgrc1 {
  my $self = shift;
  $self->{_sn_dgrc1} = shift if @_ ;
  return $self->{_sn_dgrc1};
}

sub source_name_dgrc2 {
  my $self = shift;
  $self->{_sn_dgrc2} = shift if @_ ;
  return $self->{_sn_dgrc2};
}

sub source_name_drsc {
  my $self = shift;
  $self->{_sn_drsc} = shift if @_ ;
  return $self->{_sn_drsc};
}

sub source_name_epd {
  my $self = shift;
  $self->{_sn_epd} = shift if @_ ;
  return $self->{_sn_epd};
}

sub source_name_flyreg {
  my $self = shift;
  $self->{_sn_flyreg} = shift if @_ ;
  return $self->{_sn_flyreg};
}

sub source_name_gb {
  my $self = shift;
  $self->{_sn_gb} = shift if @_ ;
  return $self->{_sn_gb};
}

sub source_name_gbprotein {
  my $self = shift;
  $self->{_sn_gbprotein} = shift if @_ ;
  return $self->{_sn_gbprotein};
}

sub source_name_gcr {
  my $self = shift;
  $self->{_sn_gcr} = shift if @_ ;
  return $self->{_sn_gcr};
}
sub source_name_gi {
  my $self = shift;
  $self->{_sn_gi} = shift if @_ ;
  return $self->{_sn_gi};
}
sub source_name_go {
  my $self = shift;
  $self->{_sn_go} = shift if @_ ;
  return $self->{_sn_go};
}

sub source_name_genomeRNAi {
  my $self = shift;
  $self->{_sn_genomeRNAi} = shift if @_ ;
  return $self->{_sn_genomeRNAi};
}

sub source_name_interpro {
  my $self = shift;
  $self->{_sn_interpro} = shift if @_ ;
  return $self->{_sn_interpro};
}

sub source_name_merops {
  my $self = shift;
  $self->{_sn_merops} = shift if @_ ;
  return $self->{_sn_merops};
}

sub source_name_miRBase {
  my $self = shift;
  $self->{_sn_miRBase} = shift if @_ ;
  return $self->{_sn_miRBase};
}

sub source_name_mitodrome {
  my $self = shift;
  $self->{_sn_mitodrome} = shift if @_ ;
  return $self->{_sn_mitodrome};
}

sub source_name_nrl3d {
  my $self = shift;
  $self->{_sn_nrl3d} = shift if @_ ;
  return $self->{_sn_nrl3d};
}

sub source_name_pdb {
  my $self = shift;
  $self->{_sn_pdb} = shift if @_ ;
  return $self->{_sn_pdb};
}

sub source_name_rfam {
  my $self = shift;
  $self->{_sn_rfam} = shift if @_ ;
  return $self->{_sn_rfam};
}

sub source_name_tf {
  my $self = shift;
  $self->{_sn_tf} = shift if @_ ;
  return $self->{_sn_tf};
}

sub source_name_uniprotsp {
  my $self = shift;
  $self->{_sn_uniprotsp} = shift if @_ ;
  return $self->{_sn_uniprotsp};
}

sub source_name_uniprottr {
  my $self = shift;
  $self->{_sn_uniprottr} = shift if @_ ;
  return $self->{_sn_uniprottr};
}

sub source_name_bdgpinsituexpr {
  my $self = shift;
  $self->{_sn_bdgpinsituexpr} = shift if @_ ;
  return $self->{_sn_bdgpinsituexpr};
}

sub source_name_dedb {
  my $self = shift;
  $self->{_sn_dedb} = shift if @_ ;
  return $self->{_sn_dedb};
}

sub source_name_drosdel {
  my $self = shift;
  $self->{_sn_drosdel} = shift if @_ ;
  return $self->{_sn_drosdel};
}

sub source_name_flygrid {
  my $self = shift;
  $self->{_sn_flygrid} = shift if @_ ;
  return $self->{_sn_flygrid};
}

sub source_name_hybrigenics {
  my $self = shift;
  $self->{_sn_hybrigenics} = shift if @_ ;
  return $self->{_sn_hybrigenics};
}

sub source_name_if {
  my $self = shift;
  $self->{_sn_if} = shift if @_ ;
  return $self->{_sn_if};
}

sub source_name_prefix_ensAGgene {
  my $self = shift;
  $self->{_sn_prefix_ensAG} = shift if @_ ;
  return $self->{_sn_prefix_ensAG};
}

sub source_name_prefix_ensAMgene {
  my $self = shift;
  $self->{_sn_prefix_ensAM} = shift if @_ ;
  return $self->{_sn_prefix_ensAM};
}

sub source_name_prefix_ensCEgene {
  my $self = shift;
  $self->{_sn_prefix_ensCE} = shift if @_ ;
  return $self->{_sn_prefix_ensCE};
}

sub source_name_prefix_ensCFgene {
  my $self = shift;
  $self->{_sn_prefix_ensCF} = shift if @_ ;
  return $self->{_sn_prefix_ensCF};
}

sub source_name_prefix_ensDMgene {
  my $self = shift;
  $self->{_sn_prefix_ensDM} = shift if @_ ;
  return $self->{_sn_prefix_ensDM};
}

sub source_name_prefix_ensDRgene {
  my $self = shift;
  $self->{_sn_prefix_ensDR} = shift if @_ ;
  return $self->{_sn_prefix_ensDR};
}

sub source_name_prefix_ensFRgene {
  my $self = shift;
  $self->{_sn_prefix_ensFR} = shift if @_ ;
  return $self->{_sn_prefix_ensFR};
}

sub source_name_prefix_ensGGgene {
  my $self = shift;
  $self->{_sn_prefix_ensGG} = shift if @_ ;
  return $self->{_sn_prefix_ensGG};
}

sub source_name_prefix_ensHSgene {
  my $self = shift;
  $self->{_sn_prefix_ensHS} = shift if @_ ;
  return $self->{_sn_prefix_ensHS};
}

sub source_name_prefix_ensMMgene {
  my $self = shift;
  $self->{_sn_prefix_ensMM} = shift if @_ ;
  return $self->{_sn_prefix_ensMM};
}

sub source_name_prefix_ensPTgene {
  my $self = shift;
  $self->{_sn_prefix_ensPT} = shift if @_ ;
  return $self->{_sn_prefix_ensPT};
}

sub source_name_prefix_ensRNgene {
  my $self = shift;
  $self->{_sn_ensRN} = shift if @_ ;
  return $self->{_sn_ensRN};
}

sub source_name_prefix_ensTNgene {
  my $self = shift;
  $self->{_sn_ensTN} = shift if @_ ;
  return $self->{_sn_ensTN};
}

sub source_name_prefix_modCBgene {
  my $self = shift;
  $self->{_sn_modCB} = shift if @_ ;
  return $self->{_sn_modCB};
}

sub source_name_prefix_modCEgene {
  my $self = shift;
  $self->{_sn_modCE} = shift if @_ ;
  return $self->{_sn_modCE};
}

sub source_name_prefix_modDDgene {
  my $self = shift;
  $self->{_sn_modDD} = shift if @_ ;
  return $self->{_sn_modDD};
}

sub gff_name{
  my $self = shift;
  $self->{_gff_name} = shift if @_ ;
  return $self->{_gff_name};
}

sub gff_dbxref{
  my $self = shift;
  $self->{_gff_dbxref} = shift if @_ ;
  return $self->{_gff_dbxref};
}

sub gff_synonym{
  my $self = shift;
  $self->{_gff_synonym} = shift if @_ ;
  return $self->{_gff_synonym};
}

sub gff_ontology{
  my $self = shift;
  $self->{_gff_ontology} = shift if @_ ;
  return $self->{_gff_ontology};
}

sub species_id {
  my $self = shift;
  $self->{_species_id} = shift if @_ ;
  return $self->{_species_id};
}

sub xrefs{
  my $self = shift;

  $self->{_xrefs} = shift if @_ ;
  return $self->{_xrefs};
}

sub add_xref{
    my ($self,$add_xref) = @_;
    push @{$self->xrefs() }, $add_xref;
    return;
}


sub direct_xrefs{
  my $self = shift;

  $self->{_direct_xrefs} = shift if @_ ;
  return $self->{_direct_xrefs};
}

sub add_direct_xref{
    my ($self,$dr) = @_;

    push @{$self->direct_xrefs() }, $dr;
    return;
}





sub line_contains_object_to_process{
  my ($self,$type_of_line) = @_;  # shoud be mRNA, gene, pseudogene, CDS,...

  for my $check_types ( @{$self->gff_object_types}) {
    if ($check_types =~/^$type_of_line$/){
      return 1;
    }
  }
  return 0;
}


=pod

=head2 gff_object_types

  Title       : gff_object_types
  Usage       : $obj->gff_object_types(array-ref)
  Function    : contains gff-type-identifiers of gff-objects which have to be processed
  Arguments   : array-ref
  Return-Val  : array-ref

=cut

sub gff_object_types{
  my $self = shift;

  $self->{_gff_object_types} = shift if @_ ;
  return $self->{_gff_object_types};
}


=pod

=head2 external_source_db_name

  Title       : external_source_db_name
  Usage       : $obj->external_source_db_name(external db name)
  Function    : returns name of hardcoded external source db name 
  Arguments   : external db name
  Return-Val  : string

=cut

sub external_source_db_name{
  my $self = shift;

  $self->{_external_source_db_name} = shift if @_ ;
  return $self->{_external_source_db_name};
}



# --------------------------------------------------------------------------------
# Get species (id and name) from file
# For UniProt files the filename is the taxonomy ID

sub get_species {
  my ($file) = @_;
  my ($taxonomy_id, $extension) = split(/\./, basename($file));
  my $sth = XrefParser::BaseParser->dbi()->prepare("SELECT species_id,name FROM species WHERE taxonomy_id=?");
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
    throw("Cannot find species corresponding to taxonomy ID " . $species_id . " - check species table\n");
  }

  return ($species_id, $species_name);

}

sub add_synonym{
  my ($self,$unique_id,$synref) = @_;
  #print "adding synonym for -$unique_id-:".join(" " , @$synref)."\n" ; ;
  ${$self->synonyms}{$unique_id}=$synref   if($synref);
  return;
}


sub get_synonyms{
  my ($self,$unique_id) = @_;

  return ${$self->synonyms}{$unique_id};
}


sub synonyms{
  my $self = shift;
  $self->{_synonyms} = shift if @_ ;
  return $self->{_synonyms};
}




sub gene_types{
  my $self = shift;

  $self->{_gene_types} = shift if @_ ;
  return $self->{_gene_types};
}

sub transcript_types{
  my $self = shift;

  $self->{_trans_types} = shift if @_ ;
  return $self->{_trans_types};
}

sub translation_types{
  my $self = shift;

  $self->{_tl_types} = shift if @_ ;
  return $self->{_tl_types};
}

 1;

 #  Drosophila v5.3 : xrefs 
 #  Gff_file        external_db_id  db_name
 #  ==
 #  Affymetrix      3120    AFFY_DrosGenome1
 #  DGRC-1  830     DGRC-1 
 #  DGRC-2  831     DGRC-2 
 #  DRSC    840     DRSC
 #  EPD     10100   EPD 
 #  FlyBase 800     flybase_gene_id
 #  FlyBase_Annotation_IDs  804     flybase_annotation_id
 #  FlyReg  850     FlyReg
 #  GB      700     EMBL
 #  GB_protein      1700    protein_id
 #  GCR     10200   GPCR
 #  GI      10900   GI
 #  GO      1000    GO
 #  GenomeRNAi      860     GenomeRNAi
 #  INTERPRO        1200    Interpro
 #  MEROPS  10300   MEROPS 
 #  MIR     10400   miRBase 
 #  MITODROME       870     MitoDrome 
 #  NRL_3D  1600    PDB
 #  PDB     1600    PDB
 #  Rfam    4200    RFAM
 #  TF      10500   TransFac 
 #  UniProt/Swiss-Prot      2200    Uniprot/SWISSPROT
 #  UniProt/TrEMBL  2000    Uniprot/SPTREMBL
 #  bdgpinsituexpr  880     BDGP_insitu_expr 
 #  dedb    890     DEDb 
 #  drosdel 881     DrosDel 
 #  flygrid 882     FlyGrid 
 #  hybrigenics     883     hybrigenics
 #  if      884     InteractiveFly
 #  ensAG   6600    Ens_Ag_gene # Anopheles gambiae
 #  ensAM   6630    Ens_Am_gene # apis mellifera?
 #  ensCE   6660    Ens_Ce_gene # C Elegans
 #  ensCF   5700    Ens_Cf_gene # Canis familiaris
 #  ensDM   6690    Ens_Dm_gene #
 #  ensDR   5800    Ens_Dr_gene # Danio rerio
 #  ensFR   6720    Ens_Fr_gene # Takifugu rubripes
 #  ensGG   6400    Ens_Gg_gene # Gallus gallus
 #  ensHS   2700    Ens_Hs_gene # Homo sapiens
 #  ensMM   5000    Ens_Mm_gene # mus musculus
 #  ensPT   6750    Ens_Pt_gene # Pan troglodytes
 #  ensRN   6200    Ens_Rn_gene # Rattus norvegicus
 #  ensTN   6810    Ens_Tn_gene # Tetraodon nigroviridis
 #  modCB   10600   modCB # InParanoid Model organism database, Caenorhabditis briggsae
 #  modCE   10700   modCE # Caenorhabditis elegans
 #  modDD   10800   modDD # Dictyostelium discoideum
