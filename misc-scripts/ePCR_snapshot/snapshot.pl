# access timdb directly and produce a tab delimited file
# for ePCR hits 

# uses timdb freeze dataset
# queries the clone directories directly
# outputs tab delimited on STDOUT
# contigname \t seq_start \t seq_end \t strand\t 1 \t
# length of the hit \t  HID
# and some progress info on STDERR

# to read this into ensembl you generate a helper table 
# create table map_feature (   contig_name varchar(40) DEFAULT '' NOT NULL, seq_start int(11) DEFAULT '0' NOT NULL,seq_end int(11) DEFAULT '0' NOT NULL,  strand int(11), hstart int(11),  hend int(11), hid varchar(40) DEFAULT '' NOT NULL, PRIMARY KEY (contig_name,hid,seq_start,seq_end));
# and do a 
# load data infile '/tmp/snapshot.txt' ignore into table map_feature( contig_name, seq_start, seq_end, strand, hstart, hend, hid );
# then you have to generate a row in analysis table
# and finally insert into feature select .... from map_feature, contig ...
# to find the internal_ids to the contig names ....

use Bio::EnsEMBL::TimDB::Obj;
use Bio::EnsEMBL::Analysis::ensConf;
use NDBM_File;

# Bio::EnsEMBL::Analysis::ensConf->import;
print join( " ",%$UNFIN_DATA_ROOT_CGP),"\n";

$timdb = Bio::EnsEMBL::TimDB::Obj->new
  ( -freeze => 1, -nogene => 1 );

print STDERR ( "Entering get_all_Clone_id\n" ); 
@cloneIds = $timdb->get_all_Clone_id;
print STDERR ( "Leaving get_all_Clone_id\n" );
# print join( " ", @cloneIds ),"\n";

print STDERR ( "Found ", scalar( @cloneIds ), " clones.\n" );
*LOG = *STDOUT;

print STDERR ( "Reading NDBM files.\n" );
for $cgp ( 'EU', 'EF', 'SU', 'SF' ) {
	my $clone;
	$cgp_dir = $UNFIN_DATA_ROOT_CGP->{$cgp};
        $contig_dbm_file = "$cgp_dir/unfinished_ana.dbm";
        my %unfin_contig;
        unless(tie(%unfin_contig,'NDBM_File',$contig_dbm_file,O_RDONLY,0644)){
		die("Error opening contig dbm file $contig_dbm_file");
	}
    	while (($key,$val) = each %unfin_contig) {
		( $clone ) = ( $key =~ /([^.]+)./ );
		push( @{$cloneContigs->{$clone}}, $key );
	}
}
print STDERR ("Done reading NDBM files\n" );


for my $clone (@cloneIds) {

	$clonecount++;
	my @moreInfo = $timdb->get_id_acc( $clone );
	$cgp = $moreInfo[2];
	$clonedir = $UNFIN_DATA_ROOT_CGP->{$cgp}."/data/".$moreInfo[1];
	@clonefiles  = map( "$clonedir/$_.ePCR", @{$cloneContigs->{$clone}} );
	# glob( "$clonedir*ePCR" );
	# print join( "\n",@clonefiles),"\n";
	# next;
	
	for my $clonefile (@clonefiles) {	
		if ( ! open( EPCR, $clonefile )) {
			print STDERR ( "Couldnt open $clonefile\n" );
			next;
		}

		while( <EPCR> ) {
  			( $contigname, $range, $accnum ) = split;
  			if( ! ( $range =~ /^(\d+)\.\.(\d+)$/ )) {
    				next;
  			} else {
    				( $start, $end ) = ( $1, $2 );
    				$strand = 1;
    				if( $end < $start ) {
      					$strand = -1;
      					( $start, $end ) = ( $end, $start );
    				}
				$length = $end-$start+1;
			}
			$contigname =~ s/([^.]*)\./$clone./g;
			$hitcount++;
			print LOG ( "$contigname\t$start\t$end\t$strand\t1\t$length\t$accnum\n" );
			if(( $clonecount % 100 ) == 0 ) {
				# exit;
				print STDERR ( "$hitcount hits on $clonecount clones.\n" );
			}
		}
		close( EPCR );
	}
}
exit;
