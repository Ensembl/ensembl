#!/usr/local/bin/perl --    # -*-Perl-*-
#
# Copyright (c) 2003 Tim Hubbard (th@sanger.ac.uk)
# Sanger Centre, Wellcome Trust Genome Campus, Cambs, UK
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation
#
# $Header$

# This is a driver script to populate and test the experimental dnac
# (compressed dna) table of ensembl.  Use -T for pod based tutorial

use strict;
use Getopt::Std;
use vars qw($opt_h $opt_T $opt_i 
	    $opt_U $opt_D $opt_P $opt_H $opt_p
	    $opt_s $opt_c $opt_C $opt_l $opt_d $opt_r $opt_n);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

getopts("hTi:U:D:P:H:p:s:cCl:dr:n:");

$|=1;

# specify defaults
my $def_U='ensadmin';
my $def_D='testdna';
my $def_P='3307';
my $def_H='127.0.0.1';
my $def_r=5;

if($opt_h){
    &help;
}elsif($opt_T){
    &help2;
}

sub help {
    print <<ENDHELP;
dna_compress.pl

  -h       for help
  -s char  string of DNA (for writing)
  -l char  label of clone, contig (for writing)
  -c       convert
  -C       compare
  -r num   number of randomisation subsequence selections [$def_r]
  -i id    clone_dbid (for comparing, converting, extracting DNA)
  -d       delete clone
  -n num   process N items

  -U user  ensembl db user [$def_U]
  -D db    ensembl db      [$def_D]
  -P port  port of mysql   [$def_P]
  -H host  host of mysql   [$def_H]
  -p pass  passwd for mysqluser
ENDHELP

    exit 0;
}

sub help2 {
    exec('perldoc', $0);
}

# defaults or options
$opt_U=$def_U unless $opt_U;
$opt_D=$def_D unless $opt_D;
$opt_P=$def_P unless $opt_P;
$opt_H=$def_H unless $opt_H;
$opt_r=$def_r unless $opt_r;

# db connection
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $opt_H,
					    -user   => $opt_U,
					    -port   => $opt_P,
					    -dbname => $opt_D,
					    -pass   => $opt_p);

my $seq_apt=$db->get_SequenceAdaptor;
my $contig_apt=$db->get_RawContigAdaptor;
my $clone_apt=$db->get_CloneAdaptor;

if($opt_C || $opt_c){
  
  # convert dna of sequences (loop over clone->contig->dna table,
  # inserting into dnac table when no entry is already present)
  # OR
  # compare raw and compressed sequences

  my $clones=$clone_apt->fetch_all();
  my $nd=0;
  my $ns=0;
  my $nc=0;
  foreach my $clone (@$clones){
    my $clone_dbid=$clone->dbID;
    next if($opt_i && $opt_i!=$clone_dbid);
    my $clone_id=$clone->id;
    my $created=$clone->created();
    print "Processing Clone $clone_id ($clone_dbid) $created\n";
    my $contigs=$clone->get_all_Contigs();
    foreach my $contig (@$contigs){
      my $contig_dbid=$contig->dbID;
      my $contig_id=$clone->id;
      my $seq=$seq_apt->fetch_by_RawContig_start_end_strand($contig,1,-1,1);
      my $len=length($seq);
      
      my $dna_id;
      # custom sql is required to get dna_id and check if there is a compressed id
      my $query="SELECT c.dna_id FROM contig c WHERE c.contig_id = ?";
      my $sth=$contig_apt->prepare($query);
      $sth->execute($contig_dbid);
      if(my $aref=$sth->fetchrow_arrayref()){
	($dna_id)=@$aref;
      }else{
	$contig_apt->thrown("No dna_id for Contig $contig_id");
      }
      
      print "  Processing Contig $contig_id ($contig_dbid) [$len] $dna_id\n";
      
      my $flag_compressed_exists=0;
      $query="SELECT d.dna_id FROM dnac d WHERE d.dna_id = ?";
      $sth=$seq_apt->prepare($query);
      $sth->execute($dna_id);
      if(my $aref=$sth->fetchrow_arrayref()){
	$flag_compressed_exists=1;
      }
      # convert
      if($opt_c){
	if($flag_compressed_exists){
	  print "    Compressed DNA entry already exists for dna_id $dna_id\n";
	  next;
	}
	# write compressed DNA record
	my $id=$seq_apt->store_compressed($seq,$created,$dna_id);
	print "    Created compressed DNA entry\n";
      }elsif($opt_C){
	if(!$flag_compressed_exists){
	  print "    Compressed DNA entry missing for dna_id $dna_id\n";
	  next;
	}
	my $seq2=$seq_apt->fetch_by_RawContig_start_end_strand($contig,1,-1,1,1);
	if($seq ne $seq2){
	  print "DIFFERENT\n";
	  $nd++;
	}else{
	  print "Same\n";
	  $ns++;
	}
	my $len2=length($seq2);
	my $lenx=30;
	if($len<$lenx){
	  $lenx=$len;
	  print substr($seq,0,$lenx)."\n".substr($seq2,0,$lenx)."\n\n";
	}else{
	  print substr($seq,0,$lenx)." .. ".substr($seq,$len-$lenx,$len)."\n";
	  print substr($seq2,0,$lenx)." .. ".substr($seq2,$len-$lenx,$len2)."\n\n";
	}

	# now test a subsequence
	my $lenr=$len-1;
	for(my $i=0;$i<$opt_r;$i++){
	  my $r=rand;
	  my $st=int($r*$lenr)+1;
	  $r=rand;
	  my $ed=int($r*($lenr-$st))+$st+1;
	  print "$lenr: $st-$ed\n";
	  $seq=$seq_apt->fetch_by_RawContig_start_end_strand($contig,$st,$ed,1);
	  $len=length($seq);
	  $seq2=$seq_apt->fetch_by_RawContig_start_end_strand($contig,$st,$ed,1,1);
	  if($seq ne $seq2){
	    print "DIFFERENT\n";
	    $nd++;
	  }else{
	    print "Same\n";
	    $ns++;
	  }
	  my $len2=length($seq2);
	  my $lenx=30;
	  if($len<$lenx){
	    $lenx=$len;
	    print substr($seq,0,$lenx)."\n".substr($seq2,0,$lenx)."\n\n";
	  }else{
	    print substr($seq,0,$lenx)." .. ".substr($seq,$len-$lenx,$len)."\n";
	    print substr($seq2,0,$lenx)." .. ".substr($seq2,$len-$lenx,$len2)."\n\n";
	  }
	}
      }
      
    }
    $nc++;
    last if($opt_n && $nc>$opt_n);
  }
  print "Same: $ns; Different: $nd\n" if $opt_C;
  exit 0;

}elsif($opt_s && $opt_l){

  # sequence and label defined, so write into db;

  my $acc=$opt_l;
  my $seq=$opt_s;
  my $ver=1;
  # clone object
  my $clone=new Bio::EnsEMBL::Clone;
  $clone->id($acc);
  $clone->htg_phase(4);
  $clone->embl_id($acc);
  $clone->version($ver);
  $clone->embl_version($ver);
  my $now = time;
  $clone->created($now);
  $clone->modified($now);

  # contig object
  my $contig = new Bio::EnsEMBL::RawContig;
  my $length = length($seq);
  $contig->name("$acc.$ver.1.$length");
  $contig->embl_offset(1);
  $contig->length($length);
  $contig->seq($seq);
  $clone->add_Contig($contig);

  my $dbclone;
  eval {
    $dbclone = $clone_apt->fetch_by_accession_version($acc,$ver);
  };
  if($dbclone){
    #$dbclone->delete_by_dbID;
    $clone_apt->remove($dbclone);
    print "remove old version of clone\n";
  }
  my $id=$clone_apt->store($clone);

  print "Clone '$acc' stored under dbid $id\n";
  
  exit 0;

}elsif($opt_d && $opt_i){

  # remove clone

  my $clone=$clone_apt->fetch_by_dbID($opt_i);
  if($clone ne ''){
    $clone_apt->remove($clone);
    print "Clone $opt_i removed\n";
  }else{
    print "Clone not found\n";
  }

}elsif($opt_i){

  foreach my $id (split(/,/,$opt_i)){
    my $clone=$clone_apt->fetch_by_dbID($id);
    my $contigs=$clone->get_all_Contigs;
    foreach my $contig (@$contigs){
      my $cid=$contig->dbID;
      print "$cid: ".$contig->seq."\n";
    }
  }
  exit 0;

}else{

    print "WARN - provide a valid contig_id (reading) or sequence to write\n";

}

__END__

=pod

=head1 NAME - name

dna_compress.pl

=head1 DESCRIPTION

Populates and tests the experimental dnac (compressed dna) table of ensembl.

=head1 SYNOPSIS

    

=head1 EXAMPLES

Create a clone/clontig/dna entry

    dna_compress.pl -s 'ACGTTTGGAANANGCCGTTNNNNACCGNGTGCTAAGCC' -l test

Clone 'test' stored under dbid 42


Create a compressed version of this entry

    dna_compress.pl -c -i 42


Compare compressed and uncompressed versions of this entry

    dna_compress.pl -C -i 42

(carries out a fetch of the entire contig and then 5 different random parts from it)


Extract sequence from DB, autosensing when only compressed is present.

    dna_compress.pl -i 42

28: ACGTTTGGAANANGCCGTTNNNNACCGNGTGCTAAGCC

    mysql> delete from dna where dna_id=42;

    dna_compress.pl -i 42

Switched to compressed DNA reading mode
28: ACGTTTGGAANANGCCGTTNNNNACCGNGTGCTAAGCC

    

=head1 FLAGS

=over 4

=item -h

Displays short help

=item -T

Displays this help message

=item -s <string>

String of DNA to write into database

=item -l <string>

Label for clonename, contigname to write into database

=item -c

Make compressed DNA entry for all clones in DB, or just those listed by dbid (-i)

=item -C

Compare uncompressed and compressed DNA entries where present, and/or listed by dbid (-i)

=item -i <num[,num]>

List of clone dbid

=back

=head1 VERSION HISTORY

=over 4

=item 14-Jul-2003

B<th> initial release

=back

=head1 BUGS

=head1 AUTHOR

B<Tim Hubbard> Email th@sanger.ac.uk

=cut
