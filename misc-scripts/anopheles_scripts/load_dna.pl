## the chunking capability of this script has been script has ben replicated by new script split_fasta_in_subslices_anopheles.pl
## the loading into db capability of this script - which was partial anyway - is no longer functional.  Use peipleine scripts load_seq_region.pl and load_agp.pl instead


use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my $dbhost;
my $dbuser;
my $dbname;
my $dbpass;
my $input;
my $assembly;


GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
	   'dbuser:s'    => \$dbuser,
	   'dbpass:s'    => \$dbpass,
	   'input:s'     => \$input,
	   'assembly:s'  => \$assembly
	   );

print STDERR "Connecting to $dbhost, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    );


if (!defined $input) {
    die "Please provide a fasta file\n";
}

my $in  = Bio::SeqIO->new(-file => $input, '-format' =>'Fasta');
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
		  $ac,
		  1,
		  0,
		  3
		  );

#set length of chunks

    my $div;
    if ($length < 25000) {
      $div = 1;
    }
    else {
      $div = int ($length/25000);
    }

    my $l = int ($length/$div);
    
    print STDERR "AC: $ac\tDIV: $div\tL: $l\n";
    
    my $prev_end;

    while ($internal_count <= $div) {

	my $total_length;

#first chunk

	if ($internal_count == 1) {
	    
	    my $actmp = $ac."_1";
	    my $t = $l;
	    
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
		      "$assembly"
		      
		      );

	    $count++;
	    $prev_end = $l+1;
	    $internal_count++;
	}
    
#chunks between first and last

	if (($internal_count > 1) && ($internal_count < $div)) {
	    my $end = $prev_end + $l;
	    my $subseq = $seq->subseq($prev_end,$end);
	    my $subseql = length($subseq);
	
	    my $t = $l + 1;

	    $total_length = $total_length + $subseql;

	    my $actmp = $ac."_".$internal_count;
	
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
				       $prev_end,
				       ); 
	    
	    print STDERR "SUB: $subseql\tL: $l\n";


#Load the assembly table

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
		      "$assembly"
		      
		      );

	    

	    $prev_end = $end+1;
	    $internal_count++;
	    $count++;

	}

#last chunk
	
	if ($internal_count == $div) {
	    
	    my $actmp = $ac."_".$internal_count;
	    my $subseq = $seq->subseq($prev_end,$length);
	    my $subseql = length($subseq);
	    
	    my $t = $l + 1;
	   
	    $total_length = $total_length + $subseql;

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
				       $prev_end,
				       ); 


#Load the assembly table

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
		      "$assembly"
		      
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
    




