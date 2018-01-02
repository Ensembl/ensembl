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

package XrefParser::WormbaseDirectParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use XrefParser::BaseParser;

use base qw( XrefParser::BaseParser );


sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files)){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my $wormbasegene_src_id = $self->get_source_id_for_source_name('wormbase_gene');
  my $wormbasegseq_src_id = $self->get_source_id_for_source_name('wormbase_gseqname');
  my $wormbaselocus_src_id = $self->get_source_id_for_source_name('wormbase_locus');
  my $wormbasetran_src_id = $self->get_source_id_for_source_name('wormbase_transcript');
  my $wormpep_src_id = $self->get_source_id_for_source_name('wormpep_id');

  my $xref_wgene_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbasegene_src_id AND species_id=$species_id");
  my $xref_gseq_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbasegseq_src_id AND species_id=$species_id");
  my $xref_wloc_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbaselocus_src_id AND species_id=$species_id");
  my $xref_wtran_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbasetran_src_id AND species_id=$species_id");
  my $xref_wpep_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormpep_src_id AND species_id=$species_id");

  my $pep_io = $self->get_filehandle($file);

  if ( !defined $pep_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my ($x_count, $d_count);

  my (%wbgene2seqid, %wbgene2loc, %tran2wbtran, %tran2wpep);

  while ( $_ = $pep_io->getline() ) {
    next if /^\/\//;
    
    my ($gseqid, $wbgeneid, $locus, $wbtranscript, $wormpep) = split(/\t/, $_);

    # Each WBGeneid should have only one sequence name and (optionally) one locus name
    $wbgene2seqid{$wbgeneid} = $gseqid;
    $wbgene2loc{$wbgeneid} = $locus if $locus ne '.';

    $tran2wbtran{$wbtranscript} = 1;
    $tran2wpep{$wbtranscript} = $wormpep if $wormpep ne '.';

  }
  $pep_io->close();

  foreach my $wbgid (keys %wbgene2seqid) {
    # reuse or create xref
    $xref_wgene_sth->execute($wbgid);
    my $xref_id = ($xref_wgene_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref({ acc        => $wbgid,
                                   label      => $wbgid,
                                   source_id  => $wormbasegene_src_id,
                                   species_id => $species_id,
                                   info_type  => "DIRECT"} );
      $x_count++;
    }
    $self->add_direct_xref($xref_id, $wbgid, "gene", "");
    $d_count++;
    
    my $gseqname = $wbgene2seqid{$wbgid};

    $xref_gseq_sth->execute($wbgid);
    $xref_id = ($xref_gseq_sth->fetchrow_array())[0];
    if (not $xref_id) {
      $xref_id = $self->add_xref({ acc        => $wbgid,
                                   label      => $gseqname,
                                   source_id  => $wormbasegseq_src_id,
                                   species_id => $species_id,
                                   info_type  => "DIRECT"} );
      $x_count++;
    }
    $self->add_direct_xref($xref_id, $wbgid, "gene", "");
    $d_count++;


    if (exists $wbgene2loc{$wbgid}) {
      my $loc_sym = $wbgene2loc{$wbgid};

      $xref_wloc_sth->execute($wbgid);    
      $xref_id = ($xref_wloc_sth->fetchrow_array())[0];
      if (!$xref_id) {
        $xref_id = $self->add_xref({ acc        => $wbgid,
                                     label      => $loc_sym,
                                     source_id  => $wormbaselocus_src_id,
                                     species_id => $species_id,
                                     info_type  => "DIRECT"} );
        $x_count++;
      }
    }
    
    # and direct xref
    $self->add_direct_xref($xref_id, $wbgid, "gene", "");
    $d_count++;
  }
  

  foreach my $tid (keys %tran2wbtran) {
    $xref_wtran_sth->execute($tid);      
    my $xref_id = ($xref_wtran_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref({ acc        => $tid,
                                   label      => $tid,
                                   source_id  => $wormbasetran_src_id,
                                   species_id => $species_id,
                                   info_type  => "DIRECT"} );
      $x_count++;
    }
    
    # and direct xref
    $self->add_direct_xref($xref_id, $tid, "transcript", "");
    $d_count++;
  }

  foreach my $tid (keys %tran2wpep) {
    my $wpep = $tran2wpep{$tid};

    $xref_wpep_sth->execute($wpep);
      
    my $xref_id = ($xref_wpep_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref({ acc        => $wpep,
                                   label      => $wpep,
                                   source_id  => $wormpep_src_id,
                                   species_id => $species_id,
                                   info_type  => "DIRECT"} );
      $x_count++;
    }

    # and direct xref
    $self->add_direct_xref($xref_id, $tid, "translation", "");
    $d_count++;
  }

  print "Added $d_count direct xrefs and $x_count xrefs\n" if($verbose);
  return 0;
}

1;
