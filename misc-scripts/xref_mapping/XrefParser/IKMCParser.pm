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

package XrefParser::IKMCParser;

use strict;
use warnings;
use Carp;
use LWP::UserAgent;

use base qw( XrefParser::BaseParser );

sub run_script {

  my ($self, $ref_arg) = @_;
  my $orig_source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $orig_source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my ($type, $my_args) = split(/:/,$file);

  my %type2id;

  foreach my $t ("No products available yet", "Vector available", "ES cells available", "Mice available"){
      my $ikmc = "IKMCs_".$t;
      $ikmc =~ s/ /_/g;
      $type2id{$t}  = $self->get_source_id_for_source_name($ikmc);
#      print $ikmc."\t".$type2id{$t}."\n";
      if(!defined( $type2id{$t})){
	die  "Could not get source id for $ikmc\n";
      }
    }

  my $xml = (<<'XXML');
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
  
  <Dataset name = "dcc" interface = "default" >
    <Attribute name = "mgi_accession_id" />
    <Attribute name = "marker_symbol" />
    <Attribute name = "vector_available" />
    <Attribute name = "escell_available" />
    <Attribute name = "mouse_available" />
    <Attribute name = "ensembl_gene_id" />
    </Dataset>
</Query>
XXML


#  print $xml."\nYO\n";
    
  my %symbols;
  my %ensembl_ids;
  my %status;
    
  my $path="http://www.i-dcc.org/biomart/martservice?";
  my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
  my $ua = LWP::UserAgent->new;

#  print "getting data from url\n";
  my $line_count=0;
  my $old_data="";
  my $chunks = 0;
  $ua->request($request,
	       sub{
		 my($data, $response) = @_;
		 if ($response->is_success) {
		   chomp $data;
		   if($data =~ /^MGI:/ and $chunks){
		     $old_data .= "\n";
		   }
		   my $data_line= $old_data.$data;
		   my @lines = split(/\n/,$data_line);
		   if(length($lines[-1]) == 0){
		     pop @lines;
		   }
		   $old_data = "";
		   my $count=0;
		   $chunks++;
		   my $max= scalar(@lines);
		   foreach my $entry (@lines){
		     $count++;
		     my @fields = split(/\t/,$entry);
		     next if (!length($entry));
		     if($count == $max){ # possible incomplete line
		       $old_data = $entry;
		       next;
		     }
		     elsif($count > $max){
		       croak "What the celery is going on here";
		     }
		     else{
		       $line_count++;
		       my $mgi_id = $fields[0];
		       foreach my $index (1..5){
			 if((!defined $fields[$index]) or ($fields[$index] eq "")){
			   $fields[$index] = 0;
			 }
		       }
#		       print {*STDERR} join(", ",@fields) ."\n";
		       if(!($mgi_id =~ /MGI:/)){
			 print "PROB1:$data_line\n";
			 print "PROB2:".join(', ',@fields)."\n";
		       }
		       $symbols{$mgi_id}=$fields[1];
		       $ensembl_ids{$mgi_id}=$fields[5];
		       if ((!defined $status{$mgi_id}) or  ($status{$mgi_id} eq '') ){
			 $status{$mgi_id} = 1;
		       }
		       if ($status{$mgi_id} < 4 && defined $fields[4] && $fields[4] == 1){
			 $status{$mgi_id} = 4;
		       }
		       elsif ($status{$mgi_id} < 3 && defined $fields[3] && $fields[3] == 1){
			 $status{$mgi_id} = 3;
		       }
		       elsif ($status{$mgi_id} < 2 && defined $fields[2] && $fields[2] == 1){
			 $status{$mgi_id} = 2;
		       }
		     }
		   }
		 }
		 else {
		   carp ("Problems with the web server: ".$response->status_line);
		   return 1;
		 }
	       },1000);

#  print "Number of chunks is $chunks\n";
  if($old_data){
    my @fields = split(/\t/,$old_data);

    $line_count++;
    #		     chop $line[5];
    my $mgi_id = $fields[0];
    if(!($mgi_id =~ /MGI:/)){
      print "PROB3:$old_data\n";
      print "PROB4:".join(', ',@fields)."\n";
    }
    $symbols{$mgi_id}=$fields[1];
    $ensembl_ids{$mgi_id}=$fields[5];
    $status{$mgi_id} = 1 if ( (!defined($status{$mgi_id}) or ($status{$mgi_id} eq '')) );
    if ($status{$mgi_id} < 4 && $fields[4] == 1){
      $status{$mgi_id} = 4;
    }
    elsif ($status{$mgi_id} < 3 && $fields[3] == 1){
      $status{$mgi_id} = 3;
    }
    elsif ($status{$mgi_id} < 2 && $fields[2] == 1){
      $status{$mgi_id} = 2;
    }
  }
#  print "obtained $line_count lines\n";

  my $parsed_count = 0;
  my $direct_count = 0;
  foreach my $acc (keys %symbols){
    my $source_id;
    $source_id = $type2id{'No products available yet'} if $status{$acc} == 1;
    $source_id = $type2id{'Vector available'} if $status{$acc} == 2;
    $source_id = $type2id{'ES cells available'} if $status{$acc} == 3;
    $source_id = $type2id{'Mice available'} if $status{$acc} == 4;

    my $label = $symbols{$acc} || $acc;
    my $ensembl_id = $ensembl_ids{$acc};
    #    print OUT "$acc\t$symbols{$acc}\t$description\t$ensembl_ids{$acc}\n";
    my $ensembl_type        = 'gene';

    ++$parsed_count;

    my $xref_id =
      $self->get_xref( $acc, $source_id, $species_id );

    if ( !defined($xref_id) || $xref_id eq '' ) {
      $xref_id =
	$self->add_xref({ acc        => $acc,
			  label      => $label,
			  source_id  => $source_id,
			  species_id => $species_id,
			  info_type  => "DIRECT"} );
    }
    next if(!defined($ensembl_ids{$acc}));
    if($ensembl_id){
      $self->add_direct_xref( $xref_id, $ensembl_id, $ensembl_type, $acc );
      $direct_count++;
    }
  }
  printf( "%d  xrefs succesfully parsed and %d direct xrefs added\n", $parsed_count, $direct_count );
  
  return 0;
} ## end sub run

1;
