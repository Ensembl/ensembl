=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefParser::RGDParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

sub run {


  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $source_sql = "select source_id from source where name = 'RGD' and priority_description = 'direct_xref'";
  my $sth = $dbi->prepare($source_sql);
  $sth->execute();
  my ($direct_source_id);
  $sth->bind_columns(\$direct_source_id);
  $sth->fetch();
  $sth->finish();

  my $file = @{$files}[0];

  my (%refseq) = %{$self->get_valid_codes("refseq",$species_id, $dbi)};

  my $rgd_io = $self->get_filehandle($file);

  if ( !defined $rgd_io ) {
    print "ERROR: Could not open $file\n";
    return 1;
  }
  my $line;
  my $found =0;
  while((!$found) and ($line = $rgd_io->getline())){ # ignore comments
    if(!($line =~ /^#/)){
      $found = 1;
    }
  };
  chomp $line;
  my @linearr = split(/\t/,$line);

  #
  #warn if sanity check fails
  #

  if($linearr[0] =~ /GENE_RDB_ID/){
   die ($linearr[0]."!= GENE_RDB_ID is not the first element in the header\n$line\n");
  }
  if($linearr[1] ne "SYMBOL"){
    die ("SYMBOL is not the second element in the header\n$line\n");
  }
  if($linearr[2] ne "NAME"){
    die ("NAME is not the third element in the header\n$line\n");
  }
  if($linearr[23] ne "GENBANK_NUCLEOTIDE"){
    die ("GENBANK_NUCLEOTIDE is not the twentysixth element in the header but ".$linearr[23]." is.\\n");
  }
  if($linearr[29] ne "OLD_SYMBOL"){
    die ("OLD_SYMBOL is not the 30th element in the header\n$line\n");
  }  
  if($linearr[37] ne "ENSEMBL_ID"){
    die ("ENSEMBL_ID is not the 38th element in the header\n$line\n");
  }

  my $sql = "insert ignore into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);    
  
  my $count= 0;
  my $ensembl_count = 0;
  my $mismatch = 0;
  my $syn_count = 0;
  while ( $line = $rgd_io->getline() ) {
    chomp $line;
    my ($rgd, $symbol, $name, $refseq,$old_name, $ensembl_id) = (split (/\t/,$line))[0,1,2,23,29, 37];
    my @nucs = split(/\;/,$refseq);
    my $done = 0;
    my $failed_list ="";
    foreach my $nuc (reverse @nucs){
      if(!$done){
	my $xref=undef; 
	if(defined($refseq{$nuc})){
	  foreach my $xref (@{$refseq{$nuc}}){
	    $done = 1;
	    my $xref_id = $self->add_dependent_xref({ master_xref_id => $xref,
						      acc            => $rgd,
						      label          => $symbol,
						      desc           => $name,
						      source_id      => $source_id,
                                                      dbi            => $dbi,
						      species_id     => $species_id} );
	    $count++;
	    my @syns  = split(/\;/,$old_name);
	    foreach my $syn(@syns){
	      $add_syn_sth->execute($xref_id, $syn);
	      $syn_count++;
	    }
	  }
	}
        if ($ensembl_id) {
          my @ensembl_ids = split(/\;/, $ensembl_id);
          $done = 1;
          foreach my $id (@ensembl_ids) {
            $ensembl_count++;
            $self->add_to_direct_xrefs({ stable_id => $id,
                                         type => 'gene',
                                         acc => $rgd,
                                         label => $symbol,
                                         desc => $name,
                                         dbi  => $dbi,
                                         source_id => $direct_source_id,
                                         species_id => $species_id} );
            my $xref_id = $self->get_xref($rgd, $direct_source_id, $species_id, $dbi);
            my @syns = split(/\;/, $old_name);
            foreach my $syn(@syns) {
              $add_syn_sth->execute($xref_id, $syn);
            }
          }
        }
	else{
	  $failed_list .= " $nuc";
	}
      }
    }

    if(!$done){
      $self->add_xref({ acc        => $rgd,
			label      => $symbol,
			desc       => $name,
			source_id  => $source_id,
			species_id => $species_id,
                        dbi        => $dbi,
			info_type  => "MISC"} );
      $mismatch++;
    }

  }

  $rgd_io->close();

  if($verbose){
    print "\t$count xrefs succesfully loaded and dependent on refseq\n";
    print "\t$mismatch xrefs added but with NO dependencies\n";
    print "\t$ensembl_count direct xrefs successfully loaded\n";
    print "added $syn_count synonyms\n";
  }
  return 0;
}

1;
