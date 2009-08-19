use strict;
use warnings;

use DBI;
use Getopt::Long;

my $file; # if not specified use stdin;
my $host;
my $user;
my $pass;
my $dbname;
my $port = 3306;



GetOptions ('file=s'      => \$file,
            'host=s'      => \$host,
            'user=s'      => \$user,
            'pass=s'      => \$pass,
            'port=s'      => \$port,
            'dbname=s'    => \$dbname);



#
# NOTE: if no file is given then STDIN is used :-)
#

if(defined($file)){
  open(GFF,"<$file") || die "Could not open $file for reading";
}
else{
  open(GFF,"-") || die "Could not open STDIN for reading";;
}


my $dbi = dbi();



my $name_sql = 'select s.seq_region_id, s.name from seq_region s, seq_region_attrib sra, attrib_type at where at.attrib_type_id = sra.attrib_type_id and at.name like "top_level" and s.seq_region_id  = sra.seq_region_id';

my $name_sth = $dbi->prepare($name_sql);
$name_sth->execute;
my ($seq_region_id, $name);
my %name_to_seq_region_id;
$name_sth->bind_columns(\$seq_region_id, \$name);
while($name_sth->fetch()){
  $name_to_seq_region_id{$name}= $seq_region_id;
}
$name_sth->finish;


my %stable_id_to_dbid; # use for gene,transcripts and exons


foreach my $type (qw(gene transcript exon)){
  my $sth = $dbi->prepare("SELECT ".$type."_id, stable_id  FROM ".$type."_stable_id");
  $sth->execute;
  my ($stable_id, $dbid);
  $sth->bind_columns(\$dbid, \$stable_id);
  while($sth->fetch){
    $stable_id_to_dbid{$stable_id} = $dbid;
  }
  $sth->finish;
}
my %exon_to_transcripts;

my $et_sth = $dbi->prepare("select exon_id, transcript_id from exon_transcript");
$et_sth->execute();
my ($e_id, $t_id);
$et_sth->bind_columns(\$e_id, \$t_id);
while($et_sth->fetch){
  if(!defined($e_id)){
    print "no exon_id?? $t_id\n";
    next;
  }
  if(defined($exon_to_transcripts{$e_id})){
    push @{$exon_to_transcripts{$e_id}}, $t_id;
  }
  else{
    $exon_to_transcripts{$e_id} = [$t_id];
  }
}
$et_sth->finish;



# Remove the previous data

my $sth = $dbi->prepare("delete from splicing_event");
$sth->execute || die "Could not delete entries from table splicing_event";

$sth = $dbi->prepare("delete from splicing_event_feature");
$sth->execute || die "Could not delete entries from table splicing_event_feature";

$sth = $dbi->prepare("delete from splicing_transcript_pair");
$sth->execute || die "Could not delete entries from table splicing_transcript_pair";



my $ins_SE_sql = "INSERT INTO splicing_event (splicing_event_id, name, gene_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, type) VALUES (?, ?, ?, ?, ?, ?, ?, ?)";
my $ins_SE_sth = $dbi->prepare($ins_SE_sql);

my $ins_SEF_sql = "INSERT INTO splicing_event_feature (splicing_event_feature_id, splicing_event_id, exon_id, transcript_id, feature_order, transcript_association, type, start, end) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)";
my $ins_SEF_sth = $dbi->prepare($ins_SEF_sql);

my $ins_STP_sql  = "INSERT INTO splicing_transcript_pair (splicing_transcript_pair_id, splicing_event_id, transcript_id_1, transcript_id_2) VALUES (?, ?, ?, ?)";
my $ins_STP_sth = $dbi->prepare($ins_STP_sql);



my %type_hash = ('CNE' => 0, 'CE' => 0, 'AFE'=> 0 ,"A5SS" => 0, "A3SS" => 0, "MXE" => 0, "EI" => 0, "II"=> 0, "IR" => 0, 'ALE' => 0, );
my %missing_type;


my $sef_id = 1;
my $stp_id = 1;





my $i = 0;
my $at_count = 0;
while( my $line = <GFF>){
  print $line;
  chomp $line;
  $i++;
  my ($chrom, $junk, $type, $start, $end, $junk2, $strand, $junk3, $feat) = split (/\t/, $line);

  if(!defined($type_hash{$type})){
    $missing_type{$type}++;
    next;
  }
  else{
    $type_hash{$type}++;
  }

  if(!defined($feat)){
    next;
  }
  my (@array) = split (/;/,$feat);
  my %feat_hash;
  foreach my $f (@array){
    my ($key, $value) = split(/=/,$f);
    $key =~ s/^ *//g;
    if(defined($feat_hash{$key})){
      $feat_hash{$key} .= "-".$value;
    }
    else{
      $feat_hash{$key} = $value;
    }
  }



  my $gene_id = undef;
  if(defined($feat_hash{"Derives_from"})){
    $gene_id = $stable_id_to_dbid{$feat_hash{"Derives_from"}};
  }
  else{
    die "No Derives from for line\n$line\n";
  }
  my $name = undef;
  if(defined($feat_hash{"ID"})){
    $name = $feat_hash{"ID"};
  }
  
  my $seq_region_id = $name_to_seq_region_id{$chrom};

  #  transcript_event_id     INT(10)  UNSIGNED NOT NULL AUTO_INCREMENT,
  #  name                    VARCHAR(40),
  #  gene_id                 INT(10) UNSIGNED NOT NULL,
  #  seq_region_id             INT(10) UNSIGNED NOT NULL,
  #  seq_region_start          INT(10) UNSIGNED NOT NULL,
  #  seq_region_end            INT(10) UNSIGNED NOT NULL,
  #  seq_region_strand         INT(10) UNSIGNED NOT NULL,
  #  type	                  ENUM('CNE','CE','AFE','A5SS','A3SS','MXE','IR','II','EI', 'AT'),

  my $strand_id = 0;
  if($strand eq "-"){
    $strand_id = -1;
  }
  elsif($strand eq "+"){
    $strand_id =1;
  }
  else{
    print STDERR "No strand info?? $_\n";
  }
  my $okay = 1;
  if(!defined($gene_id)){
    $okay = 0;
    print "ERROR: Could not get gene_id for ".$feat_hash{"Derives_from"}." from line $_\n";
  }
  if(!defined($seq_region_id)){
    $okay = 0;
    print "ERROR: Could not get seq_region_id for $chrom from line $_\n";    
  }
  $ins_SE_sth->execute($i,$name,$gene_id,$seq_region_id,$start,$end,$strand_id,$type) if($okay);
  # print SE "$i\t$name\t$gene_id\t$seq_region_id\t$start\t$end\t$strand_id\t$type\n" if($okay);

  if(defined($feat_hash{"Pair"})){
    my (@pairs) = split(/-/, $feat_hash{"Pair"});
    foreach my $pair (@pairs){
      my ($t1, $t2) = split(/,/,$pair);
      $ins_STP_sth->execute($stp_id, $i, $stable_id_to_dbid{$t1}, $stable_id_to_dbid{$t2});
#      print STP "$stp_id\t$i\t".$stable_id_to_dbid{$t1}."\t".$stable_id_to_dbid{$t2}."\n";
      $stp_id++;
    }
  }

  if($type eq "CNE"){
    if(defined($feat_hash{"ConstitutiveExons"})){
      my (@exons) = split(/,/,$feat_hash{"ConstitutiveExons"});
      foreach my $exon (@exons){
	#  splicing_event_id     INT(10)  UNSIGNED NOT NULL,
        #  splicing_event_feature_id 
	#  exon_id                 INT(10)  UNSIGNED NOT NULL,
        #  transcript_id
	#  feature_order           INT(10)  UNSIGNED NOT NULL,
	#  transcript_association  INT(10)  UNSIGNED NOT NULL,
	#  type                    ENUM('constitutive exon','exon','flanking_exon'),
	#  start                   INT(10)  UNSIGNED NOT NULL,
	#  end                     INT(10)  UNSIGNED NOT NULL,	print TEF "$i\t
	
	my $exon_id = $stable_id_to_dbid{$exon};
	my $ta = 1;
	foreach my $transcript_id (@{$exon_to_transcripts{$exon_id}}){
	  $ins_SEF_sth->execute($sef_id, $i, $exon_id, $transcript_id, 1, $ta, "constitutive_exon", $start, $end);
	  # print SEF "$sef_id\t$i\t$exon_id\t$transcript_id\t1\t$ta\tconstitutive_exon\t$start\t$end\n";
	  $ta++;
	}
      }
      $sef_id++;
    }	
  }
  elsif($type eq "CE"){
    ##FeaturesA=ENSE00001541085[ENST00000374404:ENST00000399974],ENSE00001463431[ENST00000374404:ENST00000399974];     
    my @features  = split(/,/, $feat_hash{"FeaturesA"});
    my @start_end = split(/,/, $feat_hash{"SitesA"});
    my $index = 1;
    while(defined($features[$index-1])){
      my $index_2 = index($features[$index-1], "\[");
      my $exon_id = $stable_id_to_dbid{substr($features[$index-1], 0, $index_2)};
      my (@trans) = split( ":", substr($features[$index-1],$index_2+1, -1));
      my $ta = 1;
      foreach my $t (@trans){
	my $transcript_id = $stable_id_to_dbid{$t};
	$start_end[$index-1] =~ /e\((\d+)-(\d+)/;
	$ins_SEF_sth->execute($sef_id, $i, $exon_id, $transcript_id, $index, $ta, "exon", $1, $2);
#	print SEF "$sef_id\t$i\t$exon_id\t$transcript_id\t$index\t$ta\texon\t$1\t$2\n";
	$ta++;
      }
      $sef_id++;
      $index++;
    }
    

  }
  elsif($type eq "IR"){
    my $index=1;
    ### FeaturesA=ENSE00000673410[ENST00000373020];
    ### FeaturesB=ENSE00001731724[ENST00000431386],ENSE00001805604[ENST00000431386]; 
    ### SitesA=e(99890555-99890743); 
    ### SitesB=e(99890720-99890743),e(99890555-99890587);



    ### FeaturesA=ENSE00001299015[ENST00000341376:ENST00000353205]; 
    ### FeaturesB=ENSE00001550702[ENST00000407855],ENSE00001554966[ENST00000407855]; 
    ### SitesA=e(41065096-41070145); 
    ### SitesB=e(41068761-41068917),e(41069749-41069853)
    foreach my $letter (qw(A B)){
      print "Feature $letter\n";
      my $start_end = $feat_hash{"Sites".$letter};
      my @start;
      my @end;
      if($letter eq "B"){
	$start_end  =~ /[e]\((\d+)-(\d+)\),[e]\((\d+)-(\d+)/;
        $start[0] = $1;
	$end[0] = $2;
	$start[1] = $3;
	$end[1] = $4;
      }
      else{			      
        $start_end = $feat_hash{"Sites".$letter};
        $start_end  =~ /[ie]\((\d+)-(\d+)/;
        $start[0] = $1;
        $end[0] = $2;
        $start[1] = $1;
        $end[1] = $2;
      }

      my @exons = split(/,/,$feat_hash{"Features".$letter});
      my $ec = 0;
      foreach my $exon (@exons){
	my $index_2 = index($exon, "\[");
	my $exon_id = $stable_id_to_dbid{substr($exon, 0, $index_2)};
	my (@trans) = split( ":", substr($exon,$index_2+1, -1));
	my $ta = 1;
	foreach my $t (@trans){
	  my $transcript_id = $stable_id_to_dbid{$t};
	  print  "($ec)  $sef_id\t$i\t$exon_id\t$transcript_id\t$index\t$ta\texon\t".$start[$ec]."\t".$end[$ec]."\n";
	  $ins_SEF_sth->execute($sef_id, $i, $exon_id, $transcript_id, $index, $ta, "exon", $start[$ec], $end[$ec]);
	  $ta++;
	}
      }	
      $sef_id++;
      $index++;
    }
  }
  elsif($type eq "AFE" 
     or $type eq "ALE"
     or $type eq "A5SS" 
     or $type eq "A3SS" 
     or $type eq "MXE"
     or $type eq "EI"
     or $type eq "IR"
     or $type eq "II"){
    my $index=1;
    ##FeaturesA=ENSE00001411023[ENST00000337248:ENST00000374404:ENST00000003583:ENST00000399974]; 
    ##FeaturesB=ENSE00001407060[ENST00000374409];
    ##alt FeaturesB=ENSE00001558073[ENST00000407906],ENSE00001410909[ENST00000404556];
    foreach my $letter (qw(A B)){
      my $start_end = $feat_hash{"Sites".$letter};
      $start_end  =~ /[ie]\((\d+)-(\d+)/;
      my $start = $1;
      my $end = $2;
      my @exons = split(/,/,$feat_hash{"Features".$letter});
      foreach my $exon (@exons){
	my $index_2 = index($exon, "\[");
	my $exon_id = $stable_id_to_dbid{substr($exon, 0, $index_2)};
	my (@trans) = split( ":", substr($exon,$index_2+1, -1));
	my $ta = 1;
	foreach my $t (@trans){
	  my $transcript_id = $stable_id_to_dbid{$t};
	  $ins_SEF_sth->execute($sef_id, $i, $exon_id, $transcript_id, $index, $ta, "exon", $start, $end);
#	  print SEF "$sef_id\t$i\t$exon_id\t$transcript_id\t$index\t$ta\texon\t$start\t$end\n";
	  $ta++;
	}
      }	
      $sef_id++;
      $index++;
    }
  }

  else{
    $missing_type{$type}++;
#    print STDERR "AHH do not know type $type\n";
  }
  
}

close GFF;

foreach my $type (keys %missing_type){
  print STDERR "Skipping ".$missing_type{$type}." entries of type $type as do not know what to do with them\n";
}
print STDERR "\n\n";
foreach my $type (keys %type_hash){
  print STDERR "Added ".$type_hash{$type}." entries of type $type\n";
}



sub dbi
{
    my $self = shift;

    if ( !defined $dbi || !$dbi->ping() ) {
        my $connect_string =
          sprintf( "dbi:mysql:host=%s;port=%s;database=%s",
            $host, $port, $dbname );

        $dbi =
          DBI->connect( $connect_string, $user, $pass,
            { 'RaiseError' => 1 } )
          or croak( "Can't connect to database: " . $DBI::errstr );
        $dbi->{'mysql_auto_reconnect'} = 1; # Reconnect on timeout
    }

    return $dbi;
}

##mysql -hens-research -uensadmin -pensembl ianl_homo_sapiens_core_54_36p -e"delete form splicing_event"
#mysql -hens-research -uensadmin -pensembl ianl_homo_sapiens_core_54_36p -e"load data local infile 'SE.txt' into table splicing_event"
