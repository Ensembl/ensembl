use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;

use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $host      = 'ecs1b';
my $dbuser    = 'ensadmin';
my $dbname    = 'drosophila_melanogaster_9_3';
my $dbpass    = 'ensembl';
my $path      = 'FLYBASE';

print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );

my ($dna) = @ARGV;

my $in  = Bio::SeqIO->new(-file => $dna, '-format' =>'Fasta');
my $count = 1;
my $total = 0;
while( (my $seq = $in->next_seq ) ) {
    my $internal_count = 1;
    $total++;
    my $sequence = $seq->seq;
    my $length = length($sequence);
    my $ac = $seq->id;

    my $clone = $total;
    my $chr_id = $clone;

#Load chromosome table
    my $chrsth = $db->prepare('insert into chromosome (chromosome_id,name,length) values (?,?,?)');
    $chrsth->execute(
		     $chr_id,
		     $ac,
		     $length
		     );
    

#Load clone table    
        my $sth = $db->prepare('insert into clone (clone_id,name, embl_acc, version, embl_version, htg_phase, created, modified) values(?, ?, ?, ?, ?, ?,NOW(), NOW())'); 
       
    
    $sth->execute(
		  $clone,
		  $ac,
		  'NULL',
		  1,
		  0,
		  3
		  );
		    
    
    my $div = int($length/25000);
    my $l = int ($length/$div);
    
    print STDERR "AC: $ac\tDIV: $div\tL: $l\n";
    
    my $prev_end;

    while ($internal_count <= $div) {

	my $total_length;
	if ($internal_count == 1) {
	    
	    my $actmp = $ac."_1";
	    my $t = $l;
#	    print STDERR  "AC: $actmp\nAC_CONTIG: $count\nDIV: $count\n";
#	    print STDERR  "$actmp\t$prev_end\t$length\n";
	    
	    my $subseq = $seq->subseq(1,$l);
	    my $subseql = length($subseq);
	    $total_length = $total_length + $subseql;
	    
	    print STDERR "SUB: $subseql\tL: $l\n";
	    
#Load DNA table
		my $statement = $db->prepare("
        insert into dna(sequence,created) 
        values(?, NOW())
        "); 
	    
		my $rv = $statement->execute($subseq); 
	
    
#Load contig table
		
		my $sth = $db->prepare("
        insert into contig(name, contig_id, dna_id, length, clone_id, embl_offset) 
        values(?, ?, ?, ?, ?, ?)
        "); 
	
		my $rv = $sth->execute(
				       $actmp,
				       $count,
				       $count,
				       $subseql,
				       $clone,
				       1,
				       );  

#Load the assembly table

	my $sth = $db->prepare("insert into assembly (chromosome_id,chr_start,chr_end,superctg_name,superctg_start,superctg_end,superctg_ori,contig_id,contig_start,contig_end,contig_ori,type) values (?,?,?,?,?,?,?,?,?,?,?,?)");
	$sth->execute(
		      $chr_id,
		      1,
		      $t,
		      "FPC_".$ac,
		      1,
		      $t,
		      1,
		      $count,
		      1,
		      $t,
		      1,
		      "FLYBASE"
		      
		      );

	    $count++;
	    $prev_end = $l+1;
	    $internal_count++;
	}
    


	if (($internal_count > 1) && ($internal_count < $div)) {
	    my $end = $prev_end + $l;
	    my $subseq = $seq->subseq($prev_end,$end);
	    my $subseql = length($subseq);
	
	    my $t = $l + 1;

	    $total_length = $total_length + $subseql;

	    my $actmp = $ac."_".$count;
	
	    #print STDERR  "AC: $actmp\nAC_CONTIG: $count\nDIV: $count\n";
	    #print STDERR  "$actmp\t$prev_end\t$length\n";


	       #Load DNA table
		my $statement = $db->prepare("
        insert into dna(sequence,created) 
        values(?, NOW())
        "); 
	    
		my $rv = $statement->execute($subseq); 
	
    
#Load contig table
		
		my $sth = $db->prepare("
        insert into contig(name, contig_id, dna_id, length, clone_id, embl_offset) 
        values(?, ?, ?, ?, ?, ?)
        "); 
	
		my $rv = $sth->execute(
				       $actmp,
				       $count,
				       $count,
				       $subseql,
				       $clone,
				       1,
				       ); 
	    
	    print STDERR "SUB: $subseql\tL: $l\n";

	    my $sth = $db->prepare("insert into assembly (chromosome_id,chr_start,chr_end,superctg_name,superctg_start,superctg_end,superctg_ori,contig_id,contig_start,contig_end,contig_ori,type) values (?,?,?,?,?,?,?,?,?,?,?,?)");
	$sth->execute(
		      $chr_id,
		      $prev_end,
		      $end,
		      "FPC_".$ac,
		      $prev_end,
		      $end,
		      1,
		      $count,
		      1,
		      $t,
		      1,
		      "FLYBASE"
		      
		      );

	    

	    $prev_end = $end+1;
	    $internal_count++;
	    $count++;

	}
	
	if ($internal_count == $div) {
	    
	    my $actmp = $ac."_".$count;
	    my $subseq = $seq->subseq($prev_end,$length);
	    my $subseql = length($subseq);
	    
	    my $t = $l + 1;
	   
	    $total_length = $total_length + $subseql;

#	    print STDERR  "AC: $actmp\nAC_CONTIG: $count\nDIV: $count\n";
#	    print STDERR  "$actmp\t$prev_end\t$length\n";
		
	       #Load DNA table
		my $statement = $db->prepare("
        insert into dna(sequence,created) 
        values(?, NOW())
        "); 
	    
		my $rv = $statement->execute($subseq); 
	
    
#Load contig table
		
		my $sth = $db->prepare("
        insert into contig(name, contig_id, dna_id, length, clone_id, embl_offset) 
        values(?, ?, ?, ?, ?, ?)
        "); 
	
		my $rv = $sth->execute(
				       $actmp,
				       $count,
				       $count,
				       $subseql,
				       $clone,
				       1,
				       ); 

	    my $sth = $db->prepare("insert into assembly (chromosome_id,chr_start,chr_end,superctg_name,superctg_start,superctg_end,superctg_ori,contig_id,contig_start,contig_end,contig_ori,type) values (?,?,?,?,?,?,?,?,?,?,?,?)");
	$sth->execute(
		      $chr_id,
		      $prev_end,
		      $length,
		      "FPC_".$ac,
		      $prev_end,
		      $length,
		      1,
		      $count,
		      1,
		      $subseql,
		      1,
		      "FLYBASE"
		      
		      );

	    $count++;
	    $internal_count++;
	    if ($total_length =! $l) {
		print STDERR  "TOTAL: $total_length\nLENGTH: $length\n";
		die;
	    }
	}
    }
}
    




