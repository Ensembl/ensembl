
#!/usr/local/ensembl/bin/perl

#script to set canonical transcripts for each gene
#the rule is if the cluster contains translated transcripts take the one with the 
#longest cds otherwise take the one with the longest cdna

#an example commandline is
#
# perl set_canonical_transcripts.pl -dbhost host -dbuser user -dbpass *** -dbname 
# dbname -dbport 3306 -coord_system toplevel -write
#

#To check the script has run correctly you can run the CoreForeignKeys healthcheck
# e.g ./run-healthcheck.sh -d dbname -output problem CoreForeignKeys

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

my $host;
my $port;
my $dbname;
my $user;
my $pass;
my $coord_system;
my $seq_region_name;
my $write = 0;
my $include_non_ref = 1;
my $verbose = 0;

&GetOptions( 
            'dbhost:s'      => \$host,
            'dbport:n'      => \$port,
            'dbname:s'    => \$dbname,
            'dbuser:s'    => \$user,
            'dbpass:s'      => \$pass,
            'coord_system_name:s' => \$coord_system,
            'seq_region_name:s' => \$seq_region_name,
            'write!' => \$write,
            'include_non_ref!' => \$include_non_ref,
            'verbose!' => \$verbose,
           );

unless ( $write ) {  
   print " you've not used the -write option so results will not be written into the database\n" ; 
} 

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                              -host   => $host,
                                              -user   => $user,
                                              -port   => $port,
                                              -dbname => $dbname,
                                              -pass => $pass,
                                             );


my $sa = $db->get_SliceAdaptor;
my $slices;

if($seq_region_name){
  my $slice = $sa->fetch_by_region($coord_system, $seq_region_name, $include_non_ref);
  push(@$slices, $slice);
}else{
  $slices = $sa->fetch_all($coord_system, '', $include_non_ref);
}

my $gene_update_sql = "update gene set canonical_transcript_id = ? where gene_id = ?";
my $sth = $db->dbc->prepare($gene_update_sql);
SLICE:foreach my $slice(@$slices){
  print "Getting genes for ".$slice->name."\n" if($verbose);
  my $genes = $slice->get_all_Genes(undef, undef, 1);
  my %canonical;
 GENE:foreach my $gene(@$genes){
    my $transcripts = $gene->get_all_Transcripts;
    if(@$transcripts == 1){
      $canonical{$gene->dbID} = $transcripts->[0]->dbID;
      next GENE;
    }
    my $has_translation =0;
    my $count = 0;
    my @with_translation;
    my @no_translation;
    foreach my $transcript(@$transcripts){
      if($transcript->translation && ($gene->biotype ne 'processed_transcript') 
         && ($gene->biotype ne 'pseudogene')){
        push(@with_translation, $transcript)
      }else{
        push(@no_translation, $transcript);
      }
    }

    my @sorted;
    if(@with_translation){
      my @len_and_trans;
      foreach my $trans (@with_translation) {
        my $h = { trans => $trans, len => $trans->translate->length };
        push @len_and_trans,$h;
      }
      my @tmp_sorted = sort { $b->{len} <=> $a->{len} } @len_and_trans;

      foreach my $h (@tmp_sorted) {
        #print "Adding to sorted " . $h->{trans}->dbID . "\n";
        push @sorted,$h->{trans};
      }
    }else{
      @sorted = sort {$b->length <=> $a->length} @no_translation;
    }
    $canonical{$gene->dbID} = $sorted[0]->dbID;
  }
  foreach my $gene_id (keys(%canonical)){
    my $transcript_id = $canonical{$gene_id};
    $sth->execute($transcript_id, $gene_id) if($write);
  }
}

