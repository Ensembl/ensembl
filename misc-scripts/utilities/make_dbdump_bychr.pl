#!/usr/local/bin/perl

# This script generates a full dump of an EnsEMBL core database for 
# a particular chromosome. Useful to create a small but fully functional
# EnsEMBL db

=head1 NAME

make_dbdump_bychr

=head1 SYNOPSIS

  1)Needs to be called within a new directory where you want
  all the files to be written
  2)with a user that is allowed to use mysqldump
  3)needs to be run on the host that runs the daemon

  make_dbdump_bychr -chr 'chr22'

=head1 DESCRIPTION

This script generates a full dump of the EnsEMBL database for 
a particular chromosome. Useful to create a small but fully functional
EnsEMBL db (e.g. laptop mini-mirror)
=cut

use strict;

use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $workdir = `pwd`; chomp($workdir);
my $host = "localhost";
my $port   = '';
my $dbname = undef; # force users to specify e.g. 'homo_sapiens_core_110';
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $chr = 'chr21';                      # smaller than chr22
my $lim;
# my $mysqldump = '/nfs/croc/michele/mysql/bin/mysqldump';
#                  /mysql/current/bin/mysqldump
my $mysqldump = 'mysqldump'; # in $PATH we trust

&GetOptions( 
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'chr:s'      => \$chr,
	     'workdir:s'  => \$workdir,
	     'limit:n'    => \$lim
	     );

die "no database specified; use -dbname a-core(like)-database" unless $dbname;
die "chromosome names should start with 'chr'" unless $chr =~ /^chr/;

my $pass_arg=""; $pass_arg="-p$dbpass" if $dbpass;

$workdir= "$workdir/$dbname";           # following the import scheme

my $limit;
if ($lim) {
    $limit = "limit $lim";
}

my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
$db->{RaiseError}++;                    # carp as soon as something wrong

unless (-d $workdir) {
    my $cmd = "mkdir -p $workdir";
    system($cmd) && die "``$cmd'' exited with exit status $?";
}

print STDERR "Dumping data from small tables needed full:\n";

#interpro and interpro_description in theory could be trimmed based on the hid of the protein_feature, but the join would take long, and the table is tiny!
my @small_tables =  qw(analysis analysisprocess chromosome externalDB meta
                       species interpro interpro_description);
$"=' ';
my $command = "$mysqldump -u $dbuser $pass_arg -T $workdir $dbname @small_tables";
system ($command) && die "``$command'' exited with exit status $?";

$command = "rm $workdir/*.sql";
system ($command) && die "``$command'' exited with exit status $?";

#Dump schema
print STDERR "Dumping database schema to $workdir/$dbname.sql\n";
$command = "$mysqldump -u $dbuser $pass_arg -d $dbname > $workdir/$dbname.sql";
system($command) && die "``$command'' exited with exit status $?";

my $sth = $db->prepare("select * from karyotype where chr_name = '$chr' $limit into outfile '$workdir/karyotype.txt'");
$sth->execute;

$sth = $db->prepare("select chromosome_id from chromosome where name = '$chr'");
$sth->execute;
my ($chrom) = $sth->fetchrow_array;

# my $sth = $db->prepare("select * from map_density where chromosome_id = $chrom into outfile '$workdir/map_density.txt'");
$sth = $db->prepare("select * from map_density where chrname = '$chr' into outfile '$workdir/map_density.txt'");
$sth->execute;

warn "Finding markers for chromosome $chr\n";
my $q = "select * from contig_landmarkMarker where chr_name = '$chr'";
$sth= $db->prepare($q);
$sth->execute;
open (FILE,">$workdir/contig_landmarkMarker.txt") || die "";
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    print FILE join("\t",@array)."\n";
}
close (FILE);


print STDERR "Finding golden path contigs for chromosome $chr\n";
my $golden_path_q = "select * from static_golden_path where chr_name = '$chr' $limit";
$sth = $db->prepare($golden_path_q);

$sth->execute;
my @contig_ids;
open (FILE,">$workdir/static_golden_path.txt");
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on raw_id being 3rd element in the table!
    push (@contig_ids,$array[2]);
    print FILE join("\t",@array)."\n";
}
close (FILE);
die "no contigs found for ``$golden_path_q''" unless @contig_ids;

my $contig_list = &get_inlist(0,@contig_ids);

$sth = $db->prepare("select * from contig where internal_id in $contig_list");
$sth->execute;

my @clone_ids;
my @dna_ids;
open (FILE,">$workdir/contig.txt");
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on exon id being 1st element in the table!
    push (@dna_ids,$array[6]);
    push (@clone_ids,$array[2]);
    print FILE join("\t",@array)."\n";
}
close (FILE);

my $dna_list = &get_inlist(0,@dna_ids);
@clone_ids = &unique(@clone_ids);
my $clone_list = &get_inlist(0,@clone_ids);

$sth = $db->prepare("select * from dna where id in $dna_list into outfile '$workdir/dna.txt'");
$sth->execute;

# $sth = $db->prepare("select * from contig_landmarkMarker where contig in
# $contig_list into outfile '$workdir/contig_landmarkMarker.txt'");
# $sth->execute;

$sth = $db->prepare("select * from clone where internal_id in $clone_list into outfile '$workdir/clone.txt'");
$sth->execute;

print STDERR "Getting all features...\n";
$sth = $db->prepare("select * from repeat_feature where contig in $contig_list into outfile '$workdir/repeat_feature.txt'");
$sth->execute;

$sth = $db->prepare("select * from feature where contig in $contig_list");
$sth->execute;

my @feature_ids;
open (FILE,">$workdir/feature.txt");
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on the feature id being 1st element in the table!
    push (@feature_ids,$array[0]);
    print FILE join("\t",@array)."\n";
}
close (FILE);

my $feature_list = &get_inlist(0,@feature_ids);

#I had to comment this one out because the feature list is too long for a chromsome
#, MySQL doesn't like it
#my $sth = $db->prepare("select * from fset_feature where feature in $feature_list");
#$sth->execute;
my @fset_ids;
open (FILE,">$workdir/fset_feature.txt") || die $!;
foreach my $feature (@feature_ids) {
    
    my $sth = $db->prepare("select * from fset_feature where feature = $feature");
    $sth->execute;
    
    while( (my $arr = $sth->fetchrow_arrayref()) ) {
	my @array = @$arr;
	
	#Relying on the fset id being 2nd element in the table!
	push (@fset_ids,$array[1]);
        print FILE join("\t",@array)."\n";
    }
}
close (FILE);

@fset_ids = &unique(@fset_ids);
my $fset_list = &get_inlist(0,@fset_ids);

$sth = $db->prepare("select * from fset where id in $fset_list into outfile '$workdir/fset.txt'");
$sth->execute;

#Now get partial dumps of all tables that relate to contigs
#using the @contigs array
print STDERR "Dumping gene data from chromosome $chr\n";

#Start from gene data, get all exons on these contigs, keep their ids
$sth = $db->prepare("select * from exon where contig in $contig_list");
$sth->execute;
my @exon_ids;
open (FILE,">$workdir/exon.txt") || die $!;
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on exon id being 1st element in the table!
    push (@exon_ids,$array[0]);
    print FILE join("\t",@array)."\n";
}
close (FILE);

@exon_ids = &unique(@exon_ids);
my $exon_list = &get_inlist(1,@exon_ids);
$sth = $db->prepare("select * from exon_transcript where exon in $exon_list");
$sth->execute;
my @transcript_ids;
open (FILE,">$workdir/exon_transcript.txt") || die $!;
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on transcript id being 2nd element in the table!
    push (@transcript_ids,$array[1]);
    print FILE join("\t",@array)."\n";
}
close (FILE);

@transcript_ids = &unique(@transcript_ids);
my $transcript_list = &get_inlist(1,@transcript_ids);
$sth = $db->prepare("select * from transcript where id in $transcript_list");
$sth->execute;
my @gene_ids;
my @translation_ids;
open (FILE,">$workdir/transcript.txt") || die $!;
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on gene id being 3rd element in the table!
    push (@gene_ids,$array[2]);
    push (@translation_ids,$array[3]);
    print FILE join("\t",@array)."\n";
}
close (FILE);

@gene_ids = &unique(@gene_ids);
my $gene_list = &get_inlist(1,@gene_ids);
$sth = $db->prepare("select * from gene where id in $gene_list into outfile '$workdir/gene.txt'");
$sth->execute;

$sth = $db->prepare("select * from genetype where gene_id in $gene_list into outfile '$workdir/genetype.txt'");
$sth->execute;

$sth = $db->prepare("select * from gene_description where gene_id in $gene_list into outfile '$workdir/gene_description.txt'");
$sth->execute;

my $translation_list = &get_inlist(1,@translation_ids);
$sth = $db->prepare("select * from translation where id in $translation_list into outfile '$workdir/translation.txt'");
$sth->execute;

$sth = $db->prepare("select * from protein_feature where translation in $translation_list into outfile '$workdir/protein_feature.txt'");
$sth->execute;

$sth = $db->prepare("select * from objectXref where ensembl_id in $translation_list");
$sth->execute;
my @oxref_ids;
my @xref_ids;
open (FILE,">$workdir/objectXref.txt") || die $!;
while( (my $arr = $sth->fetchrow_arrayref()) ) {
    my @array = @$arr;
    
    #Relying on xref id being 1st element in the table!
    push(@xref_ids,$array[3]);
    push(@oxref_ids,$array[0]);
    print FILE join("\t",@array)."\n";
}
close (FILE);

@xref_ids = &unique(@xref_ids);
my $xref_list = &get_inlist(0,@xref_ids);
my $oxref_list = &get_inlist(0,@oxref_ids);
$sth = $db->prepare("select * from Xref where xrefId in $xref_list into outfile '$workdir/Xref.txt'");
$sth->execute;

$sth = $db->prepare("select * from externalSynonym where xrefId in $xref_list into outfile '$workdir/externalSynonym.txt'");
$sth->execute;
$sth = $db->prepare("select * from identityXref where objectxrefId in $oxref_list into outfile '$workdir/identityXref.txt'");
$sth->execute;

$sth = $db->prepare("select * from supporting_feature where exon in $exon_list into outfile '$workdir/supporting_feature.txt'");
$sth->execute;

sub unique {
    
    my @unique;
    my %seen = ();
    foreach my $item (@_) {
	push(@unique,$item) unless $seen{$item}++;
    }
    return @unique;
}

sub get_inlist {
    my $string_flag = shift (@_);
    my $string;
    foreach my $element (@_) {
	if ($string_flag) {
	    $string .= "'$element',";
	}
	else {
	    $string .= "$element,";
	}
    }
    $string =~ s/,$//;
    return "($string)";
} 



