#! /usr/local/ensembl/bin/perl -w

# human_clone_loader.pl
#
# Read clones from AGP file and write any missing or new clones to 
# reference database.  Written specifically for the human genome.  
# Script to read one part of the two-part NCBI AGP for the human
# genome.
#
# Usually, the NCBI distribute the assembly from a human genome build
# as two files.  One is called seq_contig.md and maps NT contigs to
# genomic coordinates/chromosome.  The second file is called something
# like allcontig.agp.hs.build34 and maps clone sequences to NT contig.
# The second file is needed to determine which clones are used to
# build the full genome assembly - all these clones need to be the 
# reference database before the assembly can be loaded.
#
# USAGE:  Set the global variables at the top of this script
# to reflect your database setup.  Additionally set the global variable 
# to indicate the location of the NT contigs to clones file. Run.
#
# $DBNAME - reference database name
# $HOST   - database host
# $USER   - database user (needs write permission)
# $PORT   - database port
# $PASS   - password
#
# $INFILE - full path to the NCBI agp file (probably called 
#           something like allcontig.agp.hs.build34)
# 
# NOTE:  This script will also delete any non-golden clones.  Beware
# (or comment out the lines below that do something like:
#     print "Removing clone " . $clone->id . " [version:". $clone->embl_version . "]\n";
#    $clone->remove;
# )
#
# NOTE II:  This script assumes that NCBI agp file format has not changed.
# The expected format looks like the example below.  If the file currently
# looks different you will need to alter the portion of the script that reads
# the file - this is the very first thing the scrip does.
# NT_077402.1	1	616	1	F	AP006221.1	36116	36731	-
# NT_077402.1	617	167280	2	F	AL627309.15	241	166904	+
# NT_077911.1	1	40302	1	F	AP006222.1	1	40302	+
# NT_077912.1	1	153649	1	F	AL732372.15	1	153649	+
# NT_034471.3	1	111549	1	F	AC114498.2	1	111549	+
# NT_034471.3	111550	291116	2	F	AL669831.13	1	179567	+
# NT_077913.2	1	130979	1	F	AL390719.47	1	130979	+
# NT_077913.2	130980	241138	2	F	AL162741.44	2001	112159	+
# NT_077913.2	241139	323583	3	F	AL139287.24	2001	84445	+
# NT_077914.2	1	67923	1	F	AL391244.11	1	67923	-
# NT_077914.2	67924	106917	2	F	AL157945.21	101	39094	+
# NT_077914.2	106918	192022	3	F	AL645728.31	1	85105	-
# NT_077915.1	1	110608	1	F	AL031282.1	1	110608	+


use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;

# Global variables:

my $DBNAME = 'dan_human_apr04_ref';
my $HOST   = 'ecs4';
my $USER   = 'ensadmin';
my $PORT   = '3352';
my $PASS   = 'ensembl';

my $INFILE = 'allcontig.agp.hs.build34.altered290104';

# Read input information - determine all the clones in our golden path
# and suck up the information about how clones are broken into contigs.

open INFILE,  "< $INFILE"  or die "Can't open file $INFILE";

my %clone_data;

while (<INFILE>) {
  next if /fragment/;

  my ($seq_phase, $clone_id, $contig_start, $contig_end)
    = (split)[4, 5, 6, 7];

  $clone_data{$clone_id} = [] 
    unless $clone_data{$clone_id};

  my %contig = (clone_id     => $clone_id,
		seq_phase    => $seq_phase,
		contig_start => $contig_start,
		contig_end   => $contig_end);

  push @{$clone_data{$clone_id}}, \%contig;
}

close INFILE;

print "The golden path refers to " . scalar (keys %clone_data) . " clones.\n";

# Determine a list of clones and contigs (if any) presently in the database.

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $DBNAME,
					     -host   => $HOST,
					     -user   => $USER,
					     -pass   => $PASS,
					     -port   => $PORT);

my $clone_adaptor = $db->get_CloneAdaptor;
my $clones = $clone_adaptor->fetch_all;

print "There are " . scalar @$clones . " clones already in the database.\n";

# If any of the clones are not in the golden path, delete them now.

my %existing_clones;

foreach my $clone (@$clones) {

  my $faked_clone_id = $clone->embl_id . '.' . $clone->embl_version;

  $existing_clones{$faked_clone_id}++;

  unless ($clone_data{$faked_clone_id}) {
    print "Removing clone " . $clone->id . " [version:". $clone->embl_version . "]\n";
    $clone->remove;
  }
}

# Load any clones/contigs/dna that are in the AGP but not in the database.

my %htg_phase_enum = ('F' => 4,
		      'W' => -1,
		      'D' => -1,
		      'N' => -1);

my @clones_to_store;

my $pfetcher 
  = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new(-options => '-A');  
    # Some of these clones wont be in the current archive - need to 
    # use the -A option to scan the archives too.

my $fallback_pfetcher 
  = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new;  
    # But sometimes the -A option causes everything to come unstuck,
    # so we will use a vanilla version as a fallback.

foreach my $agp_clone (keys %clone_data){
  # Skip clone if it already exists in the database.

  next if $existing_clones{$agp_clone};

  # Make clone.

  my ($name, $embl_version) = split /\./, $agp_clone;
  print "Clone : $name $embl_version\n";

  my $htg_phase;

  if ($htg_phase_enum{$clone_data{$agp_clone}->[0]->{seq_phase}}) {
    $htg_phase = $htg_phase_enum{$clone_data{$agp_clone}->[0]->{seq_phase}}
  } else {
    $htg_phase = -1;
  }

  my $created = time;

  my $new_clone = Bio::EnsEMBL::Clone->new($db,
					   undef,
					   $name,
					   $name,
					   $embl_version,
					   $embl_version,
					   $htg_phase,
					   $created,
					   $created);

  # Fetch clone sequence

  my $clone_seq;

  eval{
    $clone_seq = $pfetcher->get_Seq_by_acc($agp_clone);
  };

  $clone_seq = $fallback_pfetcher->get_Seq_by_acc($agp_clone)
    unless $clone_seq;

  # Here we include a small fudge to expand the first and last 
  # contigs in each clone out to the clone edges (providing these 
  # are not just Ns.  If there are internal gaps these should remain 
  # between the contigs.

  my @agp_contigs 
    = sort {$a->{contig_start} <=> $b->{contig_start}} 
      @{$clone_data{$agp_clone}};

  if ($agp_contigs[0]->{contig_start} != 1) {

    unless ($clone_seq->subseq(1, 
			       ($agp_contigs[0]->{contig_start} - 1)) 
	    =~ /NNNNNNNNNN/i) {
      $agp_contigs[0]->{contig_start} = 1;
    }
  }

  if ($agp_contigs[-1]->{contig_end} != $clone_seq->length) {
    unless ($clone_seq->subseq($agp_contigs[-1]->{contig_end},
			   ($clone_seq->length - 1)) 
	=~ /NNNNNNNNNN/i) {
      $agp_contigs[-1]->{contig_end} = $clone_seq->length;
    }
  }


  # Make contig(s)

  foreach my $agp_contig (@agp_contigs){

    my $contig_id = $agp_contig->{clone_id} . '.' . $agp_contig->{contig_start} . 
      '.' . $agp_contig->{contig_end};

    my $length = ($agp_contig->{contig_end} - $agp_contig->{contig_start}) + 1;

    my $new_contig = Bio::EnsEMBL::RawContig->new(undef,
						  undef,
						  $contig_id,
						  undef,
						  $length,
						  $new_clone,
						  $agp_contig->{contig_start});


    # Attach clone subseq to contig
    $new_contig->seq($clone_seq->subseq($agp_contig->{contig_start}, 
					$agp_contig->{contig_end}));

    # Attach contig to clone
    $new_clone->add_Contig($new_contig);
  }

  # Add to array of clones that will soon be stored.

  push @clones_to_store, $new_clone;

}

# If all has gone well to this point, store the new clones.

print "About to add " . scalar @clones_to_store . " clones to the database.\n";

foreach my $clone (@clones_to_store) {
  print "Storing clone " . $clone->id . "\n";
  $clone_adaptor->store($clone);
}


print "Done.\n";
