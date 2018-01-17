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

package XrefParser::MIMParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $general_source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $general_source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];


  my %old_to_new;
  my %removed;
  my $source_id;
  my @sources;

  push @sources, $general_source_id;

  my $gene_source_id = $self->get_source_id_for_source_name("MIM_GENE", undef, $dbi);
  push @sources, $gene_source_id;
  my $morbid_source_id =  $self->get_source_id_for_source_name("MIM_MORBID", undef, $dbi);
  push @sources, $morbid_source_id;

  print "sources are:- ".join(", ",@sources)."\n" if($verbose);

  local $/ = "*RECORD*";

  my $mim_io = $self->get_filehandle($file);

  if ( !defined $mim_io ) {
    print "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $gene = 0;
  my $phenotype = 0;
  my $removed_count =0;

  $mim_io->getline();    # first record is empty with *RECORD* as the
                         # record seperator

  while ( $_ = $mim_io->getline() ) {
    #get the MIM number
    my $number = 0;
    my $label = undef;
    my $long_desc;
    my $is_morbid = 0;
    my $type =undef;
    if(/\*FIELD\*\s+NO\n(\d+)/){
      $number = $1;
      $source_id = $gene_source_id;      
      # if(/\*FIELD\*\sTI\n([\^\#\%\+\*]*)\d+(.*)\n(.*)\n\*/){
      # 	$label =$2; # taken from description as acc is meaning less
      # 	$long_desc = $2;
      # 	$long_desc .= $3 if defined $3;
      # 	$type = $1;
      # 	$label =~ s/\;\s[A-Z0-9]+$//; # strip gene name at end
      # 	$label = substr($label,0,35)." [".$type.$number."]";

      if(/\*FIELD\*\sTI(.+)\*FIELD\*\sTX/s) { # grab the whole TI field
	my $ti = $1;
        $ti =~ s/\n//g; # Remove return carriages
        # extract the 'type' and the whole description
	$ti =~ /([\^\#\%\+\*]*)\d+(.+)/s;
	my $type = $1;
	my $long_desc = $2;
        $long_desc =~ s/^\s//; # Remove white space at the start
        my @fields = split(";;", $long_desc);

        # Use the first block of text as description
	my $label = $fields[0] . " [" . $type . $number . "]";

	if($type eq "*"){ # gene only
	  $gene++;
	  $self->add_xref({ acc        => $number,
			    label      => $label,
			    desc       => $long_desc,
			    source_id  => $gene_source_id,
			    species_id => $species_id,
                            dbi        => $dbi,
			    info_type  => "DEPENDENT"} );
	}
	elsif((!defined $type) or ($type eq "") or ($type eq "#") or ($type eq "%")){ #phenotype only
	  $phenotype++;
	  $self->add_xref({ acc        => $number,
			    label      => $label,
			    desc       => $long_desc,
			    source_id  => $morbid_source_id,
			    species_id => $species_id,
                            dbi        => $dbi,
			    info_type  => "DEPENDENT"} );
	}
	elsif($type eq "+"){ # both
	  $gene++;
 	  $phenotype++;
	  $self->add_xref({ acc        => $number,
			    label      => $label,
			    desc       => $long_desc,
			    source_id  => $gene_source_id,
			    species_id => $species_id,
                            dbi        => $dbi,
			    info_type  => "DEPENDENT"} );

	  $self->add_xref({ acc        => $number,
			    label      => $label,
			    desc       => $long_desc,
			    source_id  => $morbid_source_id,
			    species_id => $species_id,
                            dbi        => $dbi,
			    info_type  => "DEPENDENT"} );
	}
	elsif($type eq "^"){
	  if(/\*FIELD\*\sTI\n[\^]\d+ MOVED TO (\d+)/){
            if ($1 eq $number) { next; }
	    $old_to_new{$number} = $1;
	  }
	  else{
	    $removed{$number} = 1;
	    $removed_count++;
	  }
	
	}
      }
    }
  }

  $mim_io->close();

  my $syn_count =0;
  foreach my $mim (keys %old_to_new){
    my $old= $mim;
    my $new= $old_to_new{$old};
    while(defined($old_to_new{$new})){
      $new = $old_to_new{$new};
    }
    if(!defined($removed{$new})){
      $self->add_to_syn_for_mult_sources($new, \@sources, $old, $species_id, $dbi);
      $syn_count++;
    }
  }
  print "$gene genemap and $phenotype phenotype MIM xrefs added\n" if($verbose);
  print "added $syn_count synonyms (defined by MOVED TO)\n" if($verbose);
  return 0; #successful
}

1;
