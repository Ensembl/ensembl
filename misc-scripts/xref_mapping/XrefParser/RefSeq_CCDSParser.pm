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

package XrefParser::RefSeq_CCDSParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;

# Parse file of Refseq records and assign direct xrefs

sub run_script {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $user = "ensro";
  my $host;
  my $port;
  my $dbname;
  my $pass;

  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }

  my $mrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_mRNA','ccds', $dbi);
  my $pred_mrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_mRNA_predicted','ccds', $dbi);

  if($verbose){
     print "RefSeq_mRNA source ID = $mrna_source_id\n";
     print "RefSeq_mRNA_predicted source ID = $pred_mrna_source_id\n" ;
   }

  # becouse the direct mapping have no descriptions etc
  # we have to steal these from the previous Refseq parser.

  my %label;
  my %version;
  my %description;

  my $sql =(<<'RSS');
SELECT xref.accession, xref.label, xref.version,  xref.description 
  FROM xref, source 
    WHERE xref.source_id = source.source_id AND
          source.name = ?
RSS

  my $sth = $dbi->prepare($sql);
  foreach my $refseq (qw(RefSeq_mRNA RefSeq_ncRNA)){
    $sth->execute($refseq);
    my ($acc, $lab, $ver, $desc);
    $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
    while (my @row = $sth->fetchrow_array()) {
      $label{$acc} = $lab;
      $version{$acc} = $ver;
      $description{$acc} = $desc;
    }
  }
  $sth->finish;
 



  my $x_sql =(<<'XSL');
SELECT x.accession, x.xref_id, d.ensembl_stable_id, "Transcript"
  FROM xref x, transcript_direct_xref d, source s
    WHERE s.source_id = x.source_id AND
          x.xref_id = d.general_xref_id AND 
          s.name like "CCDS"
XSL

  $sth = $dbi->prepare($x_sql);
  $sth->execute();
  my ($access, $old_xref_id, $stable_id, $type);
  $sth->bind_columns(\$access, \$old_xref_id, \$stable_id, \$type);
  my %ensembl_stable_id;
  my %ensembl_type;
  my %old_xref;
  while (my @row = $sth->fetchrow_array()) {
      push @{$ensembl_stable_id{$access}}, $stable_id;
      $ensembl_type{$access} = $type;
      $old_xref{$access} = $old_xref_id; 
  }
  $sth->finish;


  my $line_count = 0;
  my $xref_count = 0;
  my $direct_count = 0;
  my %seen;
  my %old_to_new;

  #
  # dbi2 is the ccds database
  #
  my $ccds_db =  XrefParser::Database->new({ host   => $host,
					port   => $port,
					user   => $user,
					dbname => $dbname,
					pass   => $pass});

  my $dbi2 = $ccds_db->dbi();

  if(!defined($dbi2)){
    return 1;
  }


##############################NEW#########################################

  # get ccds -> xref transcript_id                 ensembl_stable_id{CCDS1} = ENST00001
  # get ccds -> internal transcript_id             ccds_to_internal_id(CCDS1} = 12345

  my $ccds_sql =(<<CCDS);
SELECT x.dbprimary_acc, ox.ensembl_id
  FROM xref x, object_xref ox, external_db e
    WHERE x.xref_id = ox.xref_id AND
          x.external_db_id = e.external_db_id AND
	  ensembl_object_type = 'Transcript' AND
          e.db_name like ?
      ORDER BY x.version
CCDS

# order by version added so that the hash gets overwritten with the latest version.


  # calculate internal_id -> xref transcript_id
  my %internal_to_stable_id;
  my ($acc, $internal_id);

  $sth = $dbi2->prepare($ccds_sql);
  $sth->execute("CCDS") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $acc = $row[0];
    my $internal_id = $row[1];
    if(defined($ensembl_stable_id{$acc})){
      $internal_to_stable_id{$internal_id} =  $ensembl_stable_id{$acc};
    }
    else{
      print "$acc not found in ccds database????\n";
    }
  }    

  # for each object_xref for refseq_mRNA change internal_id to xref transcript_id
  $sth->execute("RefSeq_mRNA") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $refseq = $row[0];
    my $internal_id = $row[1];

    if(defined($internal_to_stable_id{$internal_id})){

	#If CCDS points to multiple ensembl transcripts, we don't want to store the direct RefSeq xrefs for it.
	#This would produce a many to many relationship between RefSeqs and ensembl transcripts.
	#They will be sequence mapped.
	if ( scalar(@{$internal_to_stable_id{$internal_id}}) > 1) {
	    next;
	}
    }
    else{
      print "Problem no internal_to_stable_id for $internal_id\n"; 
      next;
    }
  
    $line_count++;
    if(!defined($seen{$refseq})){
      $seen{$refseq} = 1;
      my $new_source_id;
      if ($refseq =~ /^XM_/ ){
	$new_source_id = $pred_mrna_source_id;
      } 
      elsif( $refseq =~ /^NM/) {
	$new_source_id = $mrna_source_id;
      } else {
	croak "refseq $refseq does not start with XM_ or NM_";
      }

      my $xref_id = $self->add_xref({ acc        => $refseq,
				      version    => $version{$refseq} ,
				      label      => $label{$refseq}||$refseq ,
				      desc       => $description{$refseq},
				      source_id  => $new_source_id,
				      species_id => $species_id,
                                      dbi        => $dbi,
				      info_type  => "DIRECT"} );


      foreach my $stable_id (@{$internal_to_stable_id{$internal_id}}){
	$self->add_direct_xref($xref_id, $stable_id, "Transcript", "", $dbi);
	$direct_count++;
      }

      if(defined $old_xref{$refseq}){
	$old_to_new{$old_xref{$refseq}} = $xref_id;
      }
      $xref_count++;
    }
  }
############################END NEW######################################

  #for each one seen get all its dependent xrefs and load them fro the new one too;

  my $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
  my $get_dependent_xref_sth = $dbi->prepare("SELECT dependent_xref_id, linkage_annotation "
					    .  "FROM  dependent_xref where master_xref_id = ?");

  foreach my $old_xref (keys %old_to_new){
      my $linkage;
      my $dependent_id;
      $get_dependent_xref_sth->execute($old_xref);
      $get_dependent_xref_sth->bind_columns(\$dependent_id, \$linkage);
      while(my @row = $get_dependent_xref_sth->fetchrow_array()){
	  $add_dependent_xref_sth->execute($old_to_new{$old_xref}, $dependent_id, $linkage, $source_id); 
      }   
  }


  print "Parsed $line_count RefSeq_mRNA identifiers from $file, added $xref_count xrefs and $direct_count direct_xrefs  from $line_count lines.\n" if ($verbose);


  return 0;

}

1;
