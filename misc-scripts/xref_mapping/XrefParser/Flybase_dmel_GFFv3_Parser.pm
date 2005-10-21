# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#

#2L gene [...]  ID=CG11023;Dbxref=FlyBase:FBan0011023,FlyBase:FBgn0031208;gbunit=AE003590;synonym=CG11023
#2L mRNA [...]  ID=CG11023-RA;Dbxref=FlyBase:FBtr008,FlyBase:FBgn003;dbxref_2nd=Gadfly:CG11023-RA;synonym=CG23-RA


package XrefParser::Flybase_dmel_GFFv3_Parser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use XrefParser::BaseParser;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Exception;


@ISA = qw(XrefParser::BaseParser);

my %cache_source =();


# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: Flybase_dmel_GFFv3_Parser.pm file.gff\n";
    print scalar(@ARGV);
    exit(1);
  }

  run($ARGV[0], -1);

}

# --------------------------------------------------------------------------------

sub new {
  my ($class,@args) = @_;
  my $self={};
  bless $self,$class;

  $self->external_source_db_name('flybase_gff');

  #  my @gff_obj =qw( CDS exon gene mRNA ncRNA pseudogene rRNA snRNA snoRNA tRNA );
  my @gff_obj =qw( CDS gene mRNA);
  $self->gff_object_types(\@gff_obj);

  #
  # hard-coded field separators out of gff
  #

  $self->gff_name("Name=");
  $self->gff_dbxref("Dbxref=");
  $self->gff_2nd_dbxref("dbxref_2nd=");
  $self->gff_synonym("synonym=");
  $self->gff_2nd_synonym("synonym_2nd=");

  #
  # hard-coded source-names for different objects out of ./sql/populate_metadata.sql
  #

  $self->source_name_fbgn('flybase_gene_id');         # source-name for FBgn
  $self->source_name_fbtr('flybase_transcript_id');   # source-name for FBtr
  $self->source_name_fbpp('flybase_polypeptide_id');  # source-name for FBpp
  $self->source_name_fban('flybase_annotation_id');   # source-name for FBan
  $self->source_name_symbol('flybase_synonym');
  $self->source_name_name('flybase_name');
  $self->source_name_name_prefix('FlyBaseName_');

  # gadfly-CG-ids
  $self->source_name_gadfly_gene('gadfly_gene_cgid');                # cg-id from genome annotation drosphila CG0123
  $self->source_name_gadfly_transcript('gadfly_transcript_cgid');    # cg-id from genome annotation drosphila CG0123-RA
  $self->source_name_gadfly_translation('gadfly_translation_cgid');  # cg-id from genome annotation drosphila CG0123-PA
  

  my @gene_types = qw (gene) ;
  my @translation_types = qw (CDS);
  my @transcript_types = qw (mRNA ncRNA snRNA tRNA rRNA pseudogene);

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
  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;

  my $species_name;

  if(!defined($species_id)){
    ($species_id, $species_name) = get_species($file);
  }
  $self->species_id($species_id) ;

  my $external_source_db_name = $self->external_source_db_name() ;
  my $flybase_source_id = $self->get_source($external_source_db_name);

  $self->create_xrefs($flybase_source_id, $file);

  my @xrefs = @{$self->xrefs};

  $self->relink_synonyms_to_xrefs();

  my @direct_xrefs = @{ $self->direct_xrefs } ;

  # delete previous if running directly rather than via BaseParser
  if (!defined(caller(1))) {
    print "Deleting previous xrefs for these sources\n";
    XrefParser::BaseParser->delete_by_source(\@xrefs);
  }
  print "... parsed.\n";
  print STDERR "uploading ".scalar(@xrefs)." xrefs's\n";
  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print STDERR "uploading ".scalar(@direct_xrefs)." direct-xrefs's\n";
  XrefParser::BaseParser->upload_direct_xrefs(\@direct_xrefs);

}

sub relink_synonyms_to_xrefs{
  my $self = shift;
  foreach my $x (@{$self->xrefs} ){
   $x->{SYNONYMS} = $self->get_synonyms($x->{ENSEMBL_STABLE_ID});
  }
}

# --------------------------------------------------------------------------------
# Parse file into array of xref objects
# parse

sub create_xrefs {
  my ($self, $flybase_source_id, $file) = @_;

  print STDERR "starting to parse $file...." ;
  open(GFF, $file) || die "Can't open the GFF file $file\n";

  while (<GFF>) {
    chomp;
    my @col = split /\s+/;
    if($col[3]){
     

      # test if line contains information for object wanted (CDS,mRNA,gene,..)
      if ( $self->line_contains_object_to_process( $col[2] ) ){

 	my $type = $self->set_ensembl_object_type($col[2]);

	my @desc = split /\;/,$col[8];
	my $cgid = shift @desc;
	throw("parse-error: There seems to be no Identifier: $cgid. Suspicous!") unless ($cgid=~m/ID=/);
	$cgid =~s/ID=//g;

	# set up xref-entry for EVERY single item
	foreach my $item (@desc) {


	  # make all xrefs for type "Name=" in desc-field
	  $self->make_name_xref($item,$cgid,$type);

	  $self->set_flybase_synonyms($item,$cgid);

	  # make all xrefs for type "Name=" in desc-field
	  $self->make_dbxref_xref($item,$cgid,$type);
	}
      }
    }
  }
  close (GFF);
  return;
}

sub set_ensembl_object_type{
  my ($self,$t) = @_ ; # $t is identifier in gff for object : CDS,mRNA,gene,pseudogene,snRNA,....

  for my $hc (@{ $self->gene_types } ){
    if ($t=~m/$hc/){
      return 'gene';
    }
  }
  for my $hc (@{ $self->translation_types } ){
    if ($t=~m/$hc/){
      return 'translation';
    }
  }
  for my $hc (@{ $self->transcript_types} ){
    if ($t=~m/$hc/){
      return 'transcript';
    }
  }
}


sub make_dbxref_xref{
  my ($self,$item,$cgid,$type) = @_;
  my ($xref);
  my $tg1 = $self->gff_dbxref ;
  my $tg2 = $self->gff_2nd_dbxref ;

  if ($item=~/$tg1/ || $item=~/$tg2/){
    my $dbx1 = get_fields($item,$tg1);
    my $dbx2 = get_fields($item,$tg2);
    my @dbx;
    push @dbx, @{$dbx1} if $dbx1;
    push @dbx, @{$dbx2} if $dbx2;

    foreach my $dbx (@dbx) {
      my $src_id = undef;
      my $source_type = undef;

      if ($dbx =~m/FlyBase:/){
	$dbx =~s/FlyBase://g;

	if($dbx=~m/FBgn/){
	  $src_id = $self->get_source($self->source_name_fbgn);
	}elsif ($dbx =~m/FBtr/){
	  $src_id = $self->get_source($self->source_name_fbtr);
	}elsif ($dbx =~m/FBpp/){
	  $src_id = $self->get_source($self->source_name_fbpp);
	}elsif ($dbx =~m/FBan/){
	  $src_id = $self->get_source($self->source_name_fban);
	}
      }elsif($dbx =~m/Gadfly:/){
	$dbx =~s/Gadfly://g;
	if($type eq "gene"){
	  $src_id = $self->get_source($self->source_name_gadfly_gene) ;
	}
	elsif($type eq "translation"){
	  $src_id = $self->get_source($self->source_name_gadfly_translation);
	}
	elsif($type eq "transcript"){
	  $src_id = $self->get_source($self->source_name_gadfly_transcript); 
	}
      }

      if ($src_id){ # only add xref entry for FBgn FBtr...
	my $xref ;
	$xref->{ACCESSION} = $dbx ;
	$xref->{LABEL} = $dbx;
	$xref->{SOURCE_ID} = $src_id;
	$xref->{SPECIES_ID} = $self->species_id();
	$xref->{SYNONYMS} = $self->get_synonyms($cgid);
	$self->add_xref($xref);

	if ($type){
	  my $direct_xref;
	  $direct_xref = $xref ;
	  $direct_xref->{ENSEMBL_STABLE_ID} = $cgid;
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
  my ($self,$item,$cgid) = @_; 
  my @syns;

  my $syn1 = $self->gff_synonym;
  my $syn2 = $self->gff_2nd_synonym;

  if ($item=~/$syn1/ || $item=~/$syn2/){
    my $s1 = get_fields($item,$syn1);
    my $s2 = get_fields($item,$syn2);
    my @syns;
    push @syns, @{$s1} if $s1;
    push @syns, @{$s2} if $s2;
    
    $self->add_synonym($cgid,\@syns);
    return \@syns;
  }
  return undef;
}

sub make_name_xref{
  my ($self,$item,$cgid,$type) = @_;
  my $xref=undef;
  my $target = $self->gff_name ;

  if($item=~m/$target/){  ##Name=

    my $gff_gene_name = get_fields ( $item, $target ) ;
    #print "having $$gff_gene_name[0]\n" ; 
    throw("there is more than one id for item $item\n") if $$gff_gene_name[1];
    $xref->{ACCESSION} = $$gff_gene_name[0];
    $xref->{LABEL} = $$gff_gene_name[0];
    $xref->{SPECIES_ID} = $self->species_id();
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
    $direct_xref->{ENSEMBL_STABLE_ID} = $cgid;
    $direct_xref->{ENSEMBL_TYPE} = $type;
    $direct_xref->{LINKAGE_TYPE}='bla';
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
    return undef;
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


sub source_name_symbol{
  my $self = shift;

  $self->{_source_name_symbol} = shift if @_ ;
  return $self->{_source_name_symbol};
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

sub gff_2nd_dbxref{
  my $self = shift;
  $self->{_gff_2nd_dbxref} = shift if @_ ;
  return $self->{_gff_2nd_dbxref};
}

sub gff_synonym{
  my $self = shift;
  $self->{_gff_synonym} = shift if @_ ;
  return $self->{_gff_synonym};
}

sub gff_2nd_synonym{
  my $self = shift;
  $self->{_gff_2nd_synonym} = shift if @_ ;
  return $self->{_gff_2nd_synonym};
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
    print "Taxonomy ID " . $taxonomy_id . " corresponds to species ID " . $species_id . " name " . $species_name . "\n";
  } else {
    throw("Cannot find species corresponding to taxonomy ID " . $species_id . " - check species table\n");
  }

  return ($species_id, $species_name);

}

sub add_synonym{
  my ($self,$cgid,$synref) = @_;
  #print "adding synonym for -$cgid-:".join(" " , @$synref)."\n" ; ;
  ${$self->synonyms}{$cgid}=$synref   if($synref);
  return;
}


sub get_synonyms{
  my ($self,$cgid) = @_;

  return ${$self->synonyms}{$cgid};
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
