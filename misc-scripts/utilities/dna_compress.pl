#!/usr/local/ensembl/bin/perl --    # -*-Perl-*-
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
use vars qw($opt_H $opt_T
	    $opt_u $opt_d $opt_P $opt_h $opt_p $opt_C $opt_U);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

getopts("u:HTd:p:P:h:CU");

$|=1;

# specify defaults
my $def_u='ensadmin';
my $def_P='3306';
my $def_h='127.0.0.1';

if($opt_H){
    &help;
}elsif($opt_T){
    &help2;
}

sub help {
    print <<ENDHELP;

usage:
dna_compress.pl [options]

  -H       for help
  -T       tutorial
  -d       database
  -u user  ensembl db user [$def_u]
  -d db    ensembl db
  -P port  port of mysql   [$def_P]
  -h host  host of mysql   [$def_h]
  -p pass  passwd for mysqluser
  -C       compress
  -U       uncompress
ENDHELP

    exit 0;
}

sub help2 {
    exec('perldoc', $0);
}

# defaults or options
$opt_u ||= $def_u;
$opt_h ||= $def_h;
$opt_P ||= $def_P unless $opt_P;


if(!$opt_d) {
  print STDERR "ensembl db (-d) argument is required\n";
  help();
}

if(!$opt_C && !$opt_U) {
  print STDERR "either compress (-C) or uncompress (-U) is required\n";
  help();
}


# db connection
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $opt_h,
                                            -user   => $opt_u,
                                            -port   => $opt_P,
                                            -dbname => $opt_d,
                                            -pass   => $opt_p);



my $meta_container = $db->get_MetaContainer();
my $comp_seq_adaptor = $db->get_CompressedSequenceAdaptor();

my ($compressed) = 
  @{$meta_container->list_value_by_key('sequence.compression')};

if($opt_C) {
  #compress the dna in the database
  if($compressed) {
    throw("The database meta table indicates this database already contains"
     ." sequence compression.\nIf this is not the case do the following:\n"
     ."   mysql>  delete from meta where meta_key = 'sequence.compression'\n");
  }

  #one more sanity check - make sure that dna table is populated and dnac table
  #is not
  if(get_dna_count($db) == 0) {
    print STDERR "There is no uncompressed dna in this database to compress.";
    exit(0);
  }
  if(get_dnac_count($db) != 0) {
    throw("Unexpected: the dnac table already contains some rows.\n");
  }
  
  print STDERR "Compressing sequence and writing to dnac table.\n";

  #get every slice that contains sequence
  my $slice_adaptor = $db->get_SliceAdaptor();
  foreach my $slice (@{$slice_adaptor->fetch_all('seqlevel')}) {
    my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
    $comp_seq_adaptor->store($seq_region_id, $slice->seq);
    print STDERR ".";
  }

  print STDERR "\nSetting sequence.compression meta table entry.\n";
  $meta_container->delete_key('sequence.compression');
  $meta_container->store_key_value('sequence.compression', 1);

  print STDERR "Deleting uncompressed sequence.\n";
  my $sth = $db->prepare("delete from dna");
  $sth->execute();

  print STDERR "All done.\n";

} else {
  #uncompress the dna in the database
  if(!$compressed) {
    throw("The database meta table indicates that the sequence in this "
     ." database is already uncompressed.\nIf this is not the case do "
     ." the following:\n"
     ."   mysql>  delete from meta where meta_key = 'sequence.compression';\n"
     ."   mysql>  insert into meta (meta_key,meta_value) "
     . "values('sequence.compression', 1);\n");
  }

  if(get_dnac_count($db) == 0) {
    print STDERR "There is no compressed dna in this database to uncompress.";
    exit(0);
  }
  if(get_dna_count($db) != 0) {
    throw("Unexpected: the dna table already contains some rows.\n");
  }

  print STDERR "Setting sequence.compression meta table entry.\n";
  $meta_container->delete_key('sequence.compression');
  $meta_container->store_key_value('sequence.compression', 0);

  print STDERR "Uncompressing sequence and writing to dna table.\n";

  #get every slice that contains sequence
  my $slice_adaptor = $db->get_SliceAdaptor();
  my $seq_adaptor   = $db->get_SequenceAdaptor();
  
  if($seq_adaptor->isa('CompressedSequenceAdaptor')) {
    throw("Tried to get SequenceAdaptor but got CompressedSequenceAdaptor.");
  }

  foreach my $slice (@{$slice_adaptor->fetch_all('seqlevel')}) {
    my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
    my $seq = $comp_seq_adaptor->fetch_by_Slice_start_end_strand($slice);
    $seq_adaptor->store($seq_region_id, $$seq);
    print STDERR ".";
  }

  print STDERR "\nDeleting compressed sequence.\n";
  my $sth = $db->prepare("delete from dnac");
  $sth->execute();

  print STDERR "All done.\n";
}


sub get_dna_count {
  my $db = shift;

  my $sth = $db->prepare("SELECT count(*) from dna");
  $sth->execute();
  my ($count) = $sth->fetchrow_array();
  $sth->finish();

  return $count;
}


sub get_dnac_count {
  my $db = shift;

  my $sth = $db->prepare("SELECT count(*) from dnac");
  $sth->execute();
  my ($count) = $sth->fetchrow_array();
  $sth->finish();
  
  return $count;
}



__END__

=pod

=head1 NAME - name

dna_compress.pl

=head1 DESCRIPTION

Converts compressed sequence in an ensembl database to uncompressed sequence
or converts uncompressed sequence in an ensembl database to compressed
sequence.

=head1 SYNOPSIS

    

=head1 EXAMPLES

Compress the contents of the dna table in an ensembl database:

    dna_compress.pl -C -u ensadmin -p password -h host -d homo_sapiens_core_18_34 

Uncompress the contents of the dnac table in an ensembl database

    dna_compress.pl -U -u ensadmin -p password -h host -d homo_sapiens_core_18_34 


=head1 FLAGS

=over 4

=item -h

Displays short help

=item -T

Displays this help message

=item -C

Compress the dna in the database

=item -U

Uncompress the dna in the database

=item -h

The database host

=item -p

The database password

=item -P

The database port

=item -u

The database user

=back

=head1 VERSION HISTORY

=over 4

=item 14-Jul-2003

B<th> initial release

=item 08-Oct-2003

B<mcvicker> rewrote for new schema and newly implemented sequence compression

=back

=head1 BUGS

=head1 AUTHOR

B<Tim Hubbard> Email th@sanger.ac.uk

=cut
