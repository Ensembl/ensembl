use Bio::EnsEMBL::Registry;
use strict;

my $pass         = shift;
my $fasta_file   = shift || "alt.scaf.fa";
my $mapping_file = shift || "alt.scaf.agp";
my $txt_file     = shift || "alt_scaffold_placement.txt";

my $dbname = shift || "ianl_homo_sapiens_core_57_37b";
my $host   = "ens-research";
my $user   = "ensadmin";

#connect to the database

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    '-host'    => $host,
    '-user'    => $user,
    '-pass'    => $pass,
    '-dbname'  => $dbname,
    '-species' => "load"
  );

#Load fasta sequences first. Keep a record of the first so that if things
#go wrong we can delete what has already been done.

my $sth = $dba->dbc->prepare("select max(seq_region_id) from seq_region")
  || die "Could not get max seq_region_id";


$sth->execute || die "priblkem executing";
my $max_seq_region_id;
$sth->bind_columns(\$max_seq_region_id) || die "problem binding";
$sth->fetch() || die "problem fetching";
$sth->finish;


$max_seq_region_id++;
my $start_seq_region_id = $max_seq_region_id;

print "starting new seq_region at seq_region_id of $max_seq_region_id\n";
print "\nTo reset\ndelete from dna where seq_region_id >= $start_seq_region_id\ndelete from seq_region where seq_region_id >= $start_seq_region_id\n\n";

my ($chr_coord_id, $contig_coord_id);

$sth = $dba->dbc->prepare('select coord_system_id from coord_system where name like ? and attrib like "%default_version%"')
  || die "Could not prepare coord_system sql";

$sth->execute("chromosome") || die "Could not execute coord_system sql";
$sth->bind_columns(\$chr_coord_id) || die "What a bind NOT";
$sth->fetch() || die "Could not fetch coord_system_id";

$sth->execute("contig") || die "Could not execute coord_system sql";
$sth->bind_columns(\$contig_coord_id) || die "What a bind NOT";
$sth->fetch() || die "Could not fetch coord_system_id";
$sth->finish;


my $get_seq_region_sth = $dba->dbc->prepare("select seq_region_id, length from seq_region where name like ?")
  || die "Could not prepare get_seq_region";


my $insert_seq_region_sth = $dba->dbc->prepare("insert into seq_region (seq_region_id, name, coord_system_id, length) values(?, ?, ?, ?)") || die "Could not prepare seq_region insert sql";

my $insert_dna_sth =  $dba->dbc->prepare("insert into dna (seq_region_id, sequence) values (?, ?)")
  || die "Could not prepare insert dna sql";


if(!defined($chr_coord_id) or !defined($contig_coord_id)){
  die "No coord_system_id for chromosome ($chr_coord_id) or contig ($contig_coord_id)\n";
}

#      Store contigs as seq_region + dna.

open(FASTA, "<".$fasta_file)    || die "Could not open file $fasta_file";
open(MAPPER,"<".$mapping_file)  || die "Could not open file $mapping_file";

my $seq = "";
my $name = undef;
my $version = 0;
my %name2id;

while(<FASTA>){
  chomp;
  if($_ =~ /^>(\S*)/){
#    print $1."\t".$2."\n";
    if(defined($name)){
      load_seq_region($name, \$seq);
    }
    $name = $1;
    $version = $2;
    $seq = "";
  }	
  else{
    $seq .= $_;
  }
}

if(defined($name)){
  load_seq_region($name, \$seq);
}

close FASTA;


#      Create haplotype seq_region (calulate length from mapping data)
while(<MAPPER>){
  next if($_ ~= /^#/); #ignore comemnts
}	

# add mappings

#      First the ref -> haplotype mappings
#      Then haplo -> contig mappings

# create haplotype in assembly_exception


#Do these things just have names or versions too?


sub load_seq_region{
  my $name = shift;
  my $seq =  shift;

  my ($tmp_id,$tmp_len);
  $get_seq_region_sth->execute($name)
    || die "Could not execute get seq_region for $name";
  $get_seq_region_sth->bind_columns(\$tmp_id,\$tmp_len)
    || die "Could not bind get_seq_region";
  $get_seq_region_sth->fetch();
#    || die "Could not fetch seq_region id for name $name";;

  if(defined($tmp_id)){
    if($tmp_len == length($$seq)){
      print "seq region $name already found and in database\n";
      $name2id{$name} = $tmp_id;
      return;
    }
    die "Seq region found for $name but not the same length. In datatabse length = $tmp_len and in the new fasta file ".length($$seq);
  }
  print "adding new seq region $name length of sequence = ".length($$seq)."\n";

  $insert_seq_region_sth->execute($max_seq_region_id, $name, $contig_coord_id, length($$seq))
    || die "Could not insert seq region for $name";
  $insert_dna_sth->execute($max_seq_region_id, $seq)
    || die "Could not add seq for $name";

  $name2id{$name} = $max_seq_rgion_id;

  $max_seq_region_id++;

  return;
}
