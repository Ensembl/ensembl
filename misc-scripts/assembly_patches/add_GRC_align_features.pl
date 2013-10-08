#!/usr/bin/env perl
#
# Adds GRC genomic mapping to the dna align feature table
# in the HAP pipeline with their cigar strings. Note that
# in some cases there is more than one alignment per patch.
#
# Example:
#
# perl add_GRC_align_features.pl -dbhost genebuildn \
#      -dbname homo_sapiens_core_nn_nn -dbuser user -dbpass pass \
#      -patch_release GRCh37.p8 -verbose

use strict;
use warnings;

use Getopt::Long;
use Net::FTP;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis;

$| = 1;

my $dbname         = '';
my $dbhost         = '';
my $dbuser         = 'ensro';
my $dbpass         = '';
my $dbport         = '3306';
my $patch_release  = '';
my $store          = 0;
my $external_db_id = '50692'; # GRC_alignment_import
my $syn_external_db_id = '50710'; # seq_region_synonym slice type - i.e. INSDC

my $verbose        = 0; 

my @patch_types = ('PATCH_FIX','PATCH_NOVEL');
my @dna_align_features = ();

&GetOptions( 
  'dbhost:s'                 => \$dbhost,
  'dbuser:s'                 => \$dbuser,
  'dbpass:s'                 => \$dbpass,
  'dbname:s'                 => \$dbname,
  'dbport:n'                 => \$dbport,
  'patch_release:s'          => \$patch_release,
  'external_db_id:n'         => \$external_db_id,
  'write!'                   => \$store,
  'verbose!'                 => \$verbose,
);

if(!$patch_release){
  throw ("Need to specify assembly version with -patch_release.\n");
}

# get alt_scaffold_placement.txt to generate filename that we will need
# to retrieve files from the NCBI ftp site. Populates...

my %ftp_filename=();                 

my ($content, $remote_file_handle) = "";
open($remote_file_handle, '>', \$content);

my $ftp = Net::FTP->new('ftp.ncbi.nlm.nih.gov', Debug => 0)
  or die "Can't connect to NCBI FTP: $@";

$ftp->login('anonymous', '-anonymous@')
  or die 'Cannot login ', $ftp->message;

chomp $patch_release;

my $ncbi_wd =
  "genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/".
    $patch_release.
      "/PATCHES/alt_scaffolds/";

$ftp->cwd($ncbi_wd)
  or die 'Cannot change working directory ', $ftp->message;
  
$ftp->get('alt_scaffold_placement.txt', $remote_file_handle)
  or die "get failed ", $ftp->message;

my @asp_lines = split /\n/, $content;

foreach my $asp_line (@asp_lines) {
  next if $asp_line =~ /^#/;
  my @elem = split /\t/, $asp_line;
  my $patch_name = $elem[2];
  my $file_name = $elem[3]."_".$elem[6].".gff";
  print "Filename:  $file_name\t\tPatchname:  $patch_name\n" if $verbose;
  $ftp_filename{$patch_name} = $file_name;
}

# change directory to where the GRC alignments are kept:

$ncbi_wd = "alignments"; 
$ftp->cwd($ncbi_wd) or die 'Cannot change working directory ', $ftp->message;

# hash of arrays - there way be more than one alignment per file if they
# have been manually annotated, However the GRC may change all to one
# line in the near future, in the meantime, we need to deal with them.

my %align_str =(); 

foreach my $patch (keys %ftp_filename) {
  close $remote_file_handle;
  open($remote_file_handle, '>', \$content);
  $ftp->get($ftp_filename{$patch}, $remote_file_handle)
    or die "get failed ", $ftp->message;

  my @lines = split "\n", $content;
  foreach my $line (@lines) {
    next if $line =~ /^\#/;
    # We'll parse the data later because we need most of it.
    push @{$align_str{$patch}},$line; 
  }
  #sleep 1;
}



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dbhost,
                                             -user   => $dbuser,
                                             -pass   => $dbpass,
                                             -port   => $dbport,
                                             -dbname => $dbname );

my $sa = $db->get_SliceAdaptor();

my $analysis = new Bio::EnsEMBL::Analysis( -logic_name => "grc_alignment_import",
                                           -db_version => $patch_release);



# TODO - leave $asm_exc_adaptor in - that way we can compare
# information from file with what we already know as a sanity check.


# Now get the patches, they come in pairs, the assembly exception and the reference
print "Getting patches...\n" if $verbose;
my $asm_exc_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();
my @exceptions = @{$asm_exc_adaptor->fetch_all()};
my @patches;
EXC: foreach my $exc (@exceptions){
  foreach my $type (@patch_types){
    if($exc->type() =~ m/$type/){
      push(@patches, $exc);
      next EXC;
    }
  }
}
# Assuming that AssemblyExceptionFeatureAdaptor's fetch_all will always
# return 2 entries for each patch and that those two entries are adjacent
my $num_patches = scalar(@patches)/2;

print "Have ".$num_patches." patches.\n";

# for each patch
for (my $i = 0; $i < $num_patches; $i++) {

  # get the two slices
  my $ref_slice;
  my $patch_slice;
 
  for(my $j = 0; $j < 2; $j++) {
    my $exc = pop(@patches);
    # if this is the ref version
    if($exc->type =~ m/REF/){
      # alt is only the patch slice
      $patch_slice = $exc->alternate_slice();
    }
    else{
      # alt is replaced region of ref
      $ref_slice = $exc->alternate_slice();
    }
  }
  if(!($patch_slice and $ref_slice)){
    throw("Something is wrong, the patch and ref slices were not set correctly.\n");
  }

  my @patch_vals = split /:/, $patch_slice->display_id; 
  my $patch_name = $patch_vals[2]; 

  foreach my $string ( @{ $align_str{$patch_name}}) {
    my @el = split /\t/, $string;
    
    my $num = $#el;
    throw ("Incorrect number of elements in gtf file: $num") unless $num == 8;
 
    my ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attr) = split /\t/, $string;
    
    $strand = fix_strand($strand);

    my %attribute = ();
    foreach my $kvp (split ';', $attr) {
      my ($key, $value) = split '=', $kvp;
      $attribute{$key} = $value;
    }

    my $target = $attribute{"Target"};
    my ($hseqname, $hstart, $hend, $hstrand ) = split " ", $target;
    
    $hstrand = fix_strand($hstrand); 
    
    my $length = ($hend - $hstart) + 1;    

    my $cigar_line;    
    $cigar_line = $attribute{"Gap"};
    if (defined $cigar_line) {
      sanity_check_cigar_line($cigar_line, $length);
      $cigar_line = reformat_cigar_line($cigar_line);
    } else {
      $cigar_line = $length."M";
    }

                      
    # need the seq_region_id from seq_region_synonym
    my @synonyms = @{$ref_slice->get_all_synonyms()};
    my $seq_region_id = '';
    foreach my $syn (@synonyms) {
      if ($syn->external_db_id() == $syn_external_db_id) {
        $seq_region_id = $syn->seq_region_id();
        last();
      }
    }   
    # ...to obtain the slice:

    my $slice = $sa->fetch_by_seq_region_id($seq_region_id);


    my $daf = new Bio::EnsEMBL::DnaDnaAlignFeature(
      -slice          => $slice,
      -start          => $start,
      -end            => $end,
      -strand         => $strand,
      -analysis       => $analysis,                        
      -score          => $score,
      -hstart         => $hstart,
      -hend           => $hend,
      -hstrand        => $hstrand,
      -hseqname       => $hseqname,
      -hcoverage      => $attribute{"pct_coverage"},
      -percent_id     => $attribute{"pct_identity_ungap"}, 
      -external_db_id => $external_db_id,
      -cigar_string   => $cigar_line,
    );

    push @dna_align_features, $daf;
  }
}


# now store all the dna_align features
if (scalar(@dna_align_features) > 0) {
  write_dna_align_features_to_db($db,\@dna_align_features,$store)
}
print "There are ".scalar (@dna_align_features)." dna_align_features.\n";


sub write_dna_align_features_to_db {
  my ($db,$dna_align_features,$store) = @_;

  DAF: foreach my $dna_align_feat (@$dna_align_features) {
    if ($store) {
      $db->get_DnaAlignFeatureAdaptor->store($dna_align_feat);
      if ($@) {
        throw("ERROR: Can't write dna_align_feat ".$dna_align_feat->hseqname." [$@]");
      }  else {
        print "Written ".$dna_align_feat->hseqname." on chr ".$dna_align_feat->slice->name
              ." strand ".$dna_align_feat->hstrand." with start ".$dna_align_feat->start
              ." end ".$dna_align_feat->end."\n" if $verbose;
      }
    } else {
      print "Not storing ".$dna_align_feat->hseqname."\n" if $verbose;
    }
  } # DAF 
  return 1;  
}

sub fix_strand { 
  my $strand = shift;
  if ($strand eq '+') {
    $strand = 1;
  } elsif ($strand eq '-') {
    $strand = -1;
  } else {
    throw("Strand problem :".$strand);
  }
  return $strand;
}
 
sub sanity_check_cigar_line {
  # ok, it sanity checks the GRCs idea of a cigar line which is close to GFF3 format
  my ($line, $len) = @_;
  my $cl_length = 0;
  throw("Can only sanity check cigar lines with whitespace") unless $line =~ /\s/;
  my @elements = split /\s/, $line;
  foreach my $el (@elements) {
    my ($operator, $num) = ($1, $2) if $el =~ /^(\w{1})(\d+)$/;
    if ($operator =~ /[MI]/) {
      $cl_length += $num;
    } elsif ($operator eq 'D') {
      # nothing to do
    } else {
      throw("Unknown alignment operator: $operator acting on $num");
    }
  }
  if ($cl_length != $len) {
    warn("Cigar_line length: $cl_length does not match length: $len for this line:\n$line\n\n");
  }
}

sub reformat_cigar_line {
  my $line = shift;
  # hack
  # the GRC cigar line format turns out to be back to front - fix it
  # this is a retrospective hack, with hindsight the logic of the script
  # would be different and probably incorporated into the sub above.
  my @elements = split /\s/, $line;  
  $line = ''; 
  foreach my $el (@elements) { 
    my ($operator, $num) = ($1, $2) if $el =~ /^(\w{1})(\d+)$/;   
    $line .= $num.$operator;
  }
  return $line;
}


exit;
