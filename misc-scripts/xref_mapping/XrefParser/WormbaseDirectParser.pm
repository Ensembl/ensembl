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
  my $wormbaselocus_src_id = $self->get_source_id_for_source_name('wormbase_locus');
  my $wormbasetran_src_id = $self->get_source_id_for_source_name('wormbase_transcript');
  my $wormpep_src_id = $self->get_source_id_for_source_name('wormpep_id');

  my $xref_wgene_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbasegene_src_id AND species_id=$species_id");
  my $xref_wloc_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbaselocus_src_id AND species_id=$species_id");
  my $xref_wtran_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormbasetran_src_id AND species_id=$species_id");
  my $xref_wpep_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=$wormpep_src_id AND species_id=$species_id");

  my $pep_io = $self->get_filehandle($file);

  if ( !defined $pep_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 error
  }

  my ($x_count, $d_count);

  my (%gene2wbgene, %gene2wbloc, %tran2wbtran, %tran2wpep);

  while ( $_ = $pep_io->getline() ) {
    my ($gid, $wbgeneid, $locus, $wbtranscript, $wormpep) = split(/\t/, $_);

    $gene2wbgene{$gid}->{$wbgeneid} = 1;
    $tran2wbtran{$wbtranscript}->{$wbtranscript} = 1;

    $gene2wbloc{$gid}->{$locus} = 1 if $locus ne '.';
    $tran2wpep{$wbtranscript}->{$wormpep} = 1 if $wormpep ne '.';

  }
  $pep_io->close();

  foreach my $gid (keys %gene2wbgene) {
    # reuse or create xref
    foreach my $wbgid (keys %{$gene2wbgene{$gid}}) {
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
      
      # and direct xref
      $self->add_direct_xref($xref_id, $gid, "gene", "");
      $d_count++;
    }
    
    my @locs = (exists $gene2wbloc{$gid}) ?  keys %{$gene2wbloc{$gid}} : ($gid);
    
    foreach my $wbloc (@locs) {
      $xref_wloc_sth->execute($wbloc);
      
      my $xref_id = ($xref_wloc_sth->fetchrow_array())[0];
      if (!$xref_id) {
        $xref_id = $self->add_xref({ acc        => $gid,
                                     label      => $wbloc,
                                     source_id  => $wormbaselocus_src_id,
                                     species_id => $species_id,
                                     info_type  => "DIRECT"} );
        $x_count++;
      }
      
      # and direct xref
      $self->add_direct_xref($xref_id, $gid, "gene", "");
      $d_count++;
    }
  }

  foreach my $tid (keys %tran2wbtran) {
    foreach my $wbtran (keys %{$tran2wbtran{$tid}}) {
      $xref_wtran_sth->execute($wbtran);
      
      my $xref_id = ($xref_wtran_sth->fetchrow_array())[0];
      if (!$xref_id) {
        $xref_id = $self->add_xref({ acc        => $wbtran,
                                     label      => $wbtran,
                                     source_id  => $wormbasetran_src_id,
                                     species_id => $species_id,
                                     info_type  => "DIRECT"} );
        $x_count++;
      }

      # and direct xref
      $self->add_direct_xref($xref_id, $tid, "transcript", "");
      $d_count++;
    }
  }

  foreach my $tid (keys %tran2wpep) {
    foreach my $wpep (keys %{$tran2wpep{$tid}}) {
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
  }

  print "Added $d_count direct xrefs and $x_count xrefs\n" if($verbose);
  return 0;
}

1;
