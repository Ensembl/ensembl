#!/usr/local/bin/perl --    # -*-Perl-*-
#
# Copyright (c) 1994,1995 Tim Hubbard (th@mrc-lmb.cam.ac.uk)
# Centre for Protein Engineering, MRC Centre, Cambridge, UK
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation
#
# Function:
# T: 
# D: 
#
# Known bugs:
#

use strict;
use vars qw(
	    $opt_h $opt_d $opt_n $opt_s $opt_f $opt_i $opt_c $opt_l $opt_S
	    $opt_D $opt_r $opt_R $opt_w $opt_B
	    );
require "getopts.pl";
&Getopts('hdns:f:i:cl:SD:r:Rw:B');

#use DB_File;
#use SDBM_File;
#use GDBM_File;
#use NDBM_File;

my $sel=0;
if($opt_d){
    $sel++;
}
if($opt_n){
    $sel++;
}
if($opt_s){
    $sel++;
}
if($sel > 1){
    print "only one option allowed [-s] [-d] [-n]\n";
    &help;
}

if($opt_h || !($opt_f || $opt_B)){
    &help;
}

sub help {
    print <<ENDOFTEXT;
set_dbm.pl -f file [-s] [-d] [-n] [-i id]
carry out various operations on DBM files
  -f file   specify file [required]
  -i ID     ID to set or lookup
  -s xx     set ID to value [-i required]
  -d        delete ID
  -n        create file
  -h        for help
  -c        count entries
  -l num    list <num> entries
  -D file   duplicate to this file
  -r file   reference DBM for keys
  -R        reverse read of reference
  -w file   write out contents in flat file
  -B        build all DBM files from txt versions
ENDOFTEXT
    exit 0;
}

if($opt_B){
    opendir(DIR,".") || die "cannot open directory";
    my @files=readdir(DIR);
    closedir(DIR);
    foreach my $file (@files){
	if($file=~/^(.*\.dbm)\.txt$/){
	    my $n=0;
	    my %id;
	    dbmopen(%id,$file,0666) || die "Cannot open $file";
	    open(IN,$file) || die "cannot open $file";
	    while(<IN>){
		chomp;
		my($key,$val)=split(/\t/,$_);
		$id{$key}=$val;
		$n++;
	    }
	    close(IN);
	    dbmclose(%id);
	    print "Build DBM for $file ($n entries)\n";
	}
    }
    exit 0;
}

my $file=$opt_f.".dbm";

# can only create a file with -n
unless(-e "$file.dir"){
    if(-e "$opt_f.dir"){
	$file=$opt_f;
    }elsif($opt_i && $opt_n){
	print "*warn $file doesn't exist - created\n";
    }else{
	print "*warn $file doesn't exist\n";
	exit 0;
    }
}

# open write for create
my %id;
if($opt_i && ($opt_s || $opt_d || $opt_n)){
#    tie(%id,'SDBM_File',$file,O_RDWR|O_CREAT,0666) || die "Cannot open $file";
    dbmopen(%id,$file,0666) || die "Cannot open $file";
    if($opt_s){
	$id{$opt_i}=$opt_s;
	print "ID $opt_i set to $opt_s\n";
    }elsif($opt_d){
	delete($id{$opt_i});
	print "ID $opt_i deleted\n";
    }
#    untie(%id);
    dbmclose(%id);
}elsif($opt_c || $opt_l || $opt_i || $opt_w){
#    tie(%id,'SDBM_File',$file,O_RDWR,0666) || die "Cannot open $file\n";
#    tie(%id,'GDBM_File',$file,&GDBM_WRCREAT,0666) || die "Cannot open $file\n";
#    tie(%id,'NDBM_File',$file,O_RDWR|O_CREAT,0666) || die "Cannot open $file\n";
    dbmopen(%id,$file,0666) || die "Cannot open $file";
    if($opt_i){
	if($id{$opt_i}){
	    print "$opt_i listed as ".$id{$opt_i}."\n";
	}else{
	    print "$opt_i NOT listed\n";
	}
    }else{
	my $n=0;
	my $ns=0;
	my($key,$value,$maxlen);
	if($opt_r && $opt_D){
	    my %id3;
#	    tie(%id3,'DB_File',$opt_r) || die "Cannot open $opt_r";
	    dbmopen(%id3,$opt_r,0444) || die "Cannot open $opt_r";
	    if($opt_R){
		foreach my $key (keys %id3){
		    $value=$id{$key};
		    &_action($key,$value,\$n,\$ns,\$maxlen);
		}
	    }else{
		foreach my $key (reverse keys %id3){
		    $value=$id{$key};
		    &_action($key,$value,\$n,\$ns,\$maxlen);
		}
	    }
#	    untie(%id3);
	    dbmclose(%id3);
	}else{
	    open(OUT,">$opt_w") || die "cannot open $opt_w" if $opt_w;
	    while(($key,$value)=each %id){
		print OUT "$key\t$value\n" if $opt_w;
		&_action($key,$value,\$n,\$ns,\$maxlen);
	    }
	    close(OUT);
	}
	print "$n entries found [max length is $maxlen]\n";
	dbmclose(%id) || "cannot close DBM file ok";
	#untie(%id) || "cannot close DBM file ok";
    }
}

sub _action{
    my($key,$value,$rn,$rns,$rmaxlen)=@_;
    $$rn++;
    my $len=length($value);
    if($len>$$rmaxlen){$$rmaxlen=$len;};
    my %id2;
#    tie(%id2,'DB_File',$opt_D) || die "Cannot open $opt_D" if $opt_D;
    dbmopen(%id2,$opt_D,0666);
    if($opt_D){
	if(!$value && $opt_r){
	    print "[$$rns/$$rn] $key not found - ignored\n";
	}elsif($id2{$key}){
	    if($id2{$key} ne $value){
		print "[$$rns/$$rn] $key seen again: different\n";
	    }else{
		print "[$$rns/$$rn] $key seen again: same\n";
	    }
	}else{
	    if($value=~/(\S+)/){
		$value=$1;
	    }
	    if($key=~/(\S+)/){
		$key=$1;
	    }
	    $id2{$key}=$value;
	    $$rns++;
	    print "[$$rns/$$rn] $key new\n";
	}
    }
    if($opt_l && $$rn<$opt_l){
	if($value=~/\t/){
	    print "$$rn: $key listed as:\n";
	    foreach my $value2 (split(/\t/,$value)){
		print "  $value2\n";
	    }
	}else{
	    # undocumented option
	    # oneoff to fix em->Em bug
	    if($opt_S){
		my $key2=$key;
		$key=~s/^[Ee][Mm]:/Em:/;
		if($key ne $key2){
		    delete $id{$key2};
		    $id{$key}=$value;
		    print "$$rn: $key2->$key\n";
		}
	    }else{
		print "$$rn: $key listed as $value\n";
	    }
	}
    }
    untie(%id2) || "cannot close DBM file ok";
    dbmclose(%id2);
    if($opt_D){
	my @tmp=(stat "$opt_D.pag");
	print $tmp[7]."\n";
    }
}
