use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $dbCore;
my $tmp_dir; #tmp directory to store the dumps of the data, by default /ecs2/scratch3/dani
my $affy_array; #hash containing the relation between misc_set_id => affy_array_id
my $probe_feature; #hash containing the information relevant to the affy_feature table
my $affy_probe; #hash containing the information relevant to the affy_probe table: misc_set_id -> [probeset,probename]
my $affy_probe_id = 1; #first affy_probe_id in the affy_probe table

# unbeknown to me there are actually probesets that contain the same probe twice.
# for the conversion we have to take them out, it will work better when the features
# and probes are generate from the fasta files directly

my %kill_probeset = ( '892_at' => 1 );


{
    my ($chost, $cuser, $cpass, $cport, $cdbname); #ensembl core db
    GetOptions('host=s'   => \$chost,
	       'user=s'   => \$cuser,
	       'pass=s'   => \$cpass,
	       'port=i'   => \$cport,
	       'dbname=s' => \$cdbname,
	       'tmpdir=s' => \$tmp_dir
	       );
    #by default, connect to the stagging server at ecs2:3364
    $chost ||= 'ecs2';
    $cuser ||= 'ensadmin';
    $cport ||= 3364;

    $tmp_dir ||= '/ecs2/scratch3/dani';

    usage('-cdbname argument is required.') if(!$cdbname);    
    #connection to the Core database
    $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	(-host   => $chost,
	 -user   => $cuser,
	 -pass   => $cpass,
	 -port   => $cport,
	 -dbname => $cdbname);

}

&populate_affy_array($dbCore);
print STDERR "Populated affy_array table\n";
&populate_probe_info($dbCore);

#get all the data from the misc_set table and copy it to the affy array table 
sub populate_affy_array{
    my $dbCore = shift;
    my $affy_name;
    my $misc_set_id;
    my $sth = $dbCore->dbc()->prepare(qq{SELECT misc_set_id,name from misc_set where code like 'AFFY%'
    });
    my $sth_insert = $dbCore->dbc()->prepare(qq{INSERT INTO affy_array (name) VALUES (?)
					 }
				      );
    $sth->execute();
    $sth->bind_columns(\$misc_set_id,\$affy_name);
    #copy each entry in the affy_array table and keep the affy_array_id assigned for a later use
    while ($sth->fetch){
	$sth_insert->execute($affy_name);
	$affy_array->{$misc_set_id} = $dbCore->dbc()->db_handle()->{'mysql_insertid'};	
    }
    $sth->finish;    
}

sub populate_probe_info{
  my $dbCore = shift;

  my ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, 
      $misc_set_id, $complete_probename, $mismatch);

  my ($affy_name, $affy_probeset, $probename);

  my $previous_seq_region_id = -1;
  my $previous_seq_region_start = -1;
  my $old_affy_probe;
  my $probe_set;
  print STDERR "Going to get affy information....\n";


  my $sql = 
    qq{
	SELECT STRAIGHT_JOIN seq_region_id, seq_region_start, seq_region_end, 
	                     seq_region_strand, mff.misc_set_id, ma1.value, 
                             (IF (ma2.value = 'Mismatch',1,0)) as mismatch
			  FROM misc_feature mf, misc_attrib ma1, misc_attrib ma2, 
                               misc_feature_misc_set mff, attrib_type at1, 
                               attrib_type at2, misc_set ms
			 WHERE ma2.misc_feature_id = mf.misc_feature_id 
			   AND ma2.attrib_type_id = at2.attrib_type_id 
			   AND at2.code = 'matchStatus' 
			   AND ma1.attrib_type_id = at1.attrib_type_id 
			   AND at1.code = 'probeName'
			   AND ma1.misc_feature_id = mf.misc_feature_id 
			   AND mf.misc_feature_id = mff.misc_feature_id 
	                   AND ms.misc_set_id = mff.misc_set_id 
	                   AND ms.code <> 'All_Affy'
                      ORDER BY seq_region_id,seq_region_start };


  print STDERR "Ready to create affy files\n";
  open FEATURE, ">$tmp_dir/affy_feature_$$\.txt";
  open PROBE, ">$tmp_dir/affy_probe_$$\.txt";

  my %stored_probes;
  my $current_probe_id = 1;
  my %merge_cache;
  my $merge_key;

  # merging of probes is only allowes for mismatch = 0
  for my $mismatch_process ( 0..1 ) {
    my $sth = $dbCore->dbc()->prepare( $sql );
    $sth->{mysql_use_result} = 1;
    $sth->execute();
    $sth->bind_columns(\$seq_region_id, \$seq_region_start, \$seq_region_end, 
		       \$seq_region_strand, 
		       \$misc_set_id, \$complete_probename, \$mismatch);
  
  
    my $prev_seq_region_id = -1;
    my $prev_start = -1;

    while($sth->fetch()){

      next unless( $mismatch == $mismatch_process );

      # flush the merge cache regularly
      if( $prev_start != $seq_region_start ||
	  $prev_seq_region_id != $seq_region_id ) {
	%merge_cache = ();
	$prev_start = $seq_region_start;
	$prev_seq_region_id = $seq_region_id;
      }

      my ($affy_name,$affy_probeset,$probename) = split /:/,$complete_probename,3;
      if( $kill_probeset{ $affy_probeset } ) { next; }

      # first check wether we have to store probe information
      my $probe_id = $stored_probes{ $complete_probename };
      if( ! $mismatch ) {
	$merge_key = join( "-", $seq_region_id, $seq_region_start, $seq_region_end,
			   $seq_region_strand, $mismatch, $affy_probeset ); 
      }

      if( ! defined $probe_id ) {
	# probe information needs to be stored, but new probe_id or existing one?
	$probe_id = $merge_cache{ $merge_key };
	if(( ! defined $probe_id ) || $mismatch ) {
	  $probe_id = $current_probe_id++;
	}

	print PROBE join( "\t", $probe_id, 
			  $affy_array->{$misc_set_id},
			  $affy_probeset,
			  $probename),"\n";
	$stored_probes{ $complete_probename } = $probe_id;
      }
      # at this point the probe_id is correct, it might already be clear that the 
      # feature doesnt need storing (there is already a merge cache entry for 
      # this position.


      # do we want to store the feature ?
      # if its already stored with that probe_id its in the 
      # merge_cache no addition feature is needed
      if( $mismatch ) {
	$merge_key = join( "-", $probe_id, $seq_region_id, $seq_region_start, $seq_region_end,
			   $seq_region_strand ); 
      }

      if( exists $merge_cache{ $merge_key } ) {
	# this one is already stored
      } else {
	$merge_cache{ $merge_key } = $probe_id;
	print FEATURE join ("\t",$seq_region_id, $seq_region_start,
			    $seq_region_end, $seq_region_strand,
			    $mismatch, $probe_id ),"\n";
      }
    } 

    $sth->finish();
  }

  close FEATURE;
  close PROBE;

  #and finally import the information
  print STDERR "Loading new affy information\n";
  load($dbCore,"$tmp_dir/affy_feature_$$\.txt",qw(affy_feature seq_region_id seq_region_start seq_region_end seq_region_strand mismatches affy_probe_id));
  load($dbCore,"$tmp_dir/affy_probe_$$\.txt",qw(affy_probe affy_probe_id affy_array_id probeset name));
}


sub load{
    my $dbCore = shift;
    my $file = shift;
    my $tablename = shift;
    my @colnames = @_;
    
    my $cols = join( ",", @colnames );
    my $sql = qq{
	LOAD DATA INFILE '$file'
	    INTO TABLE $tablename ($cols)
	};
    $dbCore->dbc()->do($sql);
    unlink ("$file");
}

sub usage {
    my $msg = shift;
    
    print STDERR <<EOF;
    
  usage: perl affy_data.pl <options>
      
    options:
      -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
      -cuser <user>        username of core Ensembl MySQL database (default = ensadmin)
      -cpass <pass>        password of core Ensembl MySQL database
      -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
      -cdbname <dbname>    dbname of core Ensembl MySQL database
EOF

      die("\n$msg\n\n");
}
