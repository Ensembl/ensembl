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
{
    my ($chost, $cuser, $cpass, $cport, $cdbname); #ensembl core db
    GetOptions('chost=s'   => \$chost,
	       'cuser=s'   => \$cuser,
	       'cpass=s'   => \$cpass,
	       'cport=i'   => \$cport,
	       'cdbname=s' => \$cdbname,
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

    my ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $misc_set_id, $value, $mismatch);
    my ($affy_name, $affy_probeset, $probename); #attributes in the value field for an attrib_type_id = 5
    my $previous_seq_region_id = 0;
    my $previous_seq_region_start = 0;
    my $old_affy_probe;
    my $probe_set;
    print STDERR "Going to get affy information....\n";
    my $sth = $dbCore->dbc()->prepare(qq{SELECT STRAIGHT_JOIN seq_region_id, seq_region_start, seq_region_end, seq_region_strand, mff.misc_set_id, ma1.value, (IF (ma2.value = 'Mismatch',1,0)) as mismatch
				      FROM misc_feature mf, misc_attrib ma1, misc_attrib ma2, misc_feature_misc_set mff, attrib_type at1, attrib_type at2, misc_set ms
				      WHERE ma2.misc_feature_id = mf.misc_feature_id 
				      AND ma2.attrib_type_id = at2.attrib_type_id 
				      AND at2.code = 'matchStatus' 
				      AND ma1.attrib_type_id = at1.attrib_type_id 
				      AND at1.code = 'probeName'
				      AND ma1.misc_feature_id = mf.misc_feature_id 
				      AND mf.misc_feature_id = mff.misc_feature_id 
				      AND ms.misc_set_id = mff.misc_set_id 
				      AND ms.code <> 'All_Affy'
				      ORDER BY seq_region_id,seq_region_start
				  });
    $sth->{mysql_use_result} = 1;
    $sth->execute();
    $sth->bind_columns(\$seq_region_id, \$seq_region_start, \$seq_region_end, \$seq_region_strand, \$misc_set_id, \$value, \$mismatch);
    print STDERR "Ready to create affy files\n";
    open FEATURE, ">$tmp_dir/affy_feature_$$\.txt";
    open PROBE, ">$tmp_dir/affy_probe_$$\.txt";
    while($sth->fetch()){
	#we have a new probe, add the previous one to the database, and flush the structures
	unless ((($previous_seq_region_id == $seq_region_id) or ($previous_seq_region_id == 0)) and (($previous_seq_region_start == $seq_region_start) or ($previous_seq_region_start == 0))){
	    foreach my $key (keys %{$affy_probe}){
		if (!exists $probe_set->{$affy_probeset . ":" . $affy_probe->{$key}}){
		    $probe_set->{$affy_probeset . ":" . $affy_probe->{$key}} = $affy_probe_id;
		    print PROBE join ("\t",$affy_probe_id,$affy_array->{$key},$affy_probeset,$affy_probe->{$key}),"\n";
		}
		$old_affy_probe = $probe_set->{$affy_probeset . ":" . $affy_probe->{$key}};
	    }
	    #insert all the affy_probe values in the file	    
	    print FEATURE join ("\t",$probe_feature->{'seq_region_id'},$probe_feature->{'seq_region_start'},$probe_feature->{'seq_region_end'}, $probe_feature->{'seq_region_strand'},$probe_feature->{'mismatches'},$old_affy_probe),"\n";
	    $affy_probeset = '';
	    $affy_probe_id++;
	    $affy_probe = ();
	    $probe_feature = ();
	}
	$previous_seq_region_id = $seq_region_id;
	$previous_seq_region_start = $seq_region_start;
	($affy_name,$affy_probeset,$probename) = split /:/,$value,3;
	$affy_probe->{$misc_set_id} = $probename;
	$probe_feature->{'seq_region_id'} = $seq_region_id;
	$probe_feature->{'seq_region_start'} = $seq_region_start;
	$probe_feature->{'seq_region_end'} = $seq_region_end;
	$probe_feature->{'seq_region_strand'} = $seq_region_strand;
	$probe_feature->{'mismatches'} = $mismatch;
    }
    $sth->finish();
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
	LOAD DATA LOCAL INFILE '$file'
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
