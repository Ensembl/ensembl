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

package AssemblyMapper::BlastzAligner;

=head1 NAME

BlastzAligner.pm - module to do a whole genome alignment between two closely
related assemblies and create assembly entries from it.

=head1 SYNOPSIS

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);
my $aligner = BlastzAligner->new(-SUPPORT => $support);

# create a tempdir for storing input and output
$aligner->create_tempdir;

# write sequences to fasta and nib files
$aligner->write_sequence(
    $alt_slice,
    $support->param('altassembly'), 
    "alt_sequence.1"
);
$aligner->write_sequence(
    $ref_slice,
    $support->param('assembly'), 
    "ref_sequence.1"
);

# run blastz
$aligner->run_blastz("alt_sequence.1", "ref_sequence.1");

=head1 DESCRIPTION

This modules contains helper functions to generate a whole genome alignment
between two closely related assemblies using blastz from scratch. Alignments
are then stored in an Ensembl assembly table.

The module is part of a series of scripts to create a mapping between two
assemblies. See "Related scripts" below for an overview of the whole process.

=head1 RELATED SCRIPTS

The whole process of creating a whole genome alignment is done by these
scripts:

    ensembl/misc-scripts/assembly/load_alternative_assembly.pl
    ensembl/misc-scripts/assembly/align_by_clone_identity.pl
    ensembl/misc-scripts/assembly/align_nonident_regions.pl

See documention in the respective script for more information.


=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use File::Path;

use constant FMT1 => "%-30s%10.0f (%3.2f%%)\n";
use constant FMT2 => "%-30s%10.0f\n";
use constant FMT3 => "%-8s%-12s%-5s%-10s%-10s%-10s%-10s\n";
use constant FMT4 => "%-8s%-12s%-5s%8.0f  %8.0f  %8.0f  %8.0f\n";

=head2 new

  Arg[-SUPPORT] : a Bio::EnsEMBL::Utils::ConversionSupport object
  Example       : my $support = new Bio::EnsEMBL::Utils::ConversionSupport(
                    $SERVERROOT);
                  my $aligner = BlastzAligner->new(-SUPPORT => $support);
  Description   : object constructor method
  Return type   : a BlastzAligner object
  Exceptions    : none
  Caller        : general

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = {};
    bless $self, $class;

    my ($support) = rearrange([qw(SUPPORT)], @args);
    $self->support($support);

    # set bindir
    $self->bindir($self->support->param('bindir') || '/usr/local/ensembl/bin');

    return $self;
}

=head2 create_tempdir

  Arg[1]      : String $tmpdir - temporary directory name
  Example     : $aligner->create_tempdir;
  Description : Creates a temporary directory in /tmp with a semi-random name
                (username.timestamp.randomnumber). Alternatively, you can pass
                it the name of a directory to use.
  Return type : String - name of the tempdir created
  Exceptions  : Thrown if tempdir can't be created
  Caller      : general

=cut

sub create_tempdir {
    my $self = shift;
    my $tempdir = shift;

    if ($tempdir) {
      unless (-d $tempdir) {
        $self->support->log_error("Can't find tempdir $tempdir: $!");
      }

    } else {
      # create tmpdir to store input and output
      my $user = `whoami`;
      chomp $user;
      $tempdir = "/tmp/$user.".time.".".int(rand(100000));
      $self->support->log("Creating tmpdir $tempdir...\n");
      system("mkdir $tempdir") == 0 or
          $self->support->log_error("Can't create tmp dir $tempdir: $!\n");
      $self->support->log("Done.\n");
    }

    $self->tempdir($tempdir);
    return $tempdir;
}

=head2 remove_tempdir

  Example     : $aligner->remove_tempdir;
  Description : Removes a temporary directory created by $self->create_tempdir.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub remove_tempdir {
    my $self = shift;
    rmtree($self->tempdir) or $self->support->log_warning("Could not delete ".$self->tempdir.": $!\n");
}

=head2 write_sequence

  Arg[1]      : Bio::EnsEMBL::Slice $slice - slice for which to write sequence
  Arg[2]      : String $assembly - assembly name
  Arg[3]      : String $basename1 - basename of single sequence file
  Arg[4]      : (optional) String $basename2 - basename of multiple sequence
                file
  Example     : $aligner->write_sequence($slice, 'NCBI35', 'e_seq.1');
  Description : Writes a slice's sequence to a fasta file and converts it to nib
                format. Optionally appends the sequence to another
                multi-sequence fasta file.
  Return type : none
  Exceptions  : thrown if faToNib or file appending fails
  Caller      : general

=cut

sub write_sequence {
    my ($self, $slice, $assembly, $basename1, $basename2) = @_;

    my $tmpdir = $self->tempdir;

    unless (-e "$tmpdir/$basename1.fa") {
      my $fh = $self->support->filehandle('>', "$tmpdir/$basename1.fa");
      my $cs_name = $slice->coord_system_name;
      print $fh join(':', ">$basename1 dna:chromfrag $cs_name",
                            $assembly,
                            $slice->start,
                            $slice->end,
                            $slice->strand
                      ), "\n";
      print $fh $slice->get_repeatmasked_seq(undef, 1)->seq, "\n";
      close($fh);
    }
    
    # convert fasta to nib (needed for lavToAxt)
    unless (-e "$tmpdir/$basename1.nib") {
      system($self->bindir."/faToNib $tmpdir/$basename1.fa $tmpdir/$basename1.nib") == 0
          or $self->support->log_error("Can't run faToNib: $!\n");
    }

    if ($basename2) {
        system("cat $tmpdir/$basename1.fa >> $tmpdir/$basename2.fa") == 0 or
            $self->support->log_error("Can't concat fasta files: $!\n");
    }
}

=head2 run_blastz

  Arg[1]      : String $A_basename - basename of alternative fasta file
  Arg[2]      : String $R_basename - basename of reference fasta file
  Example     : $aligner->run_blastz('alt_seq.1', 'ref_seq.1');
  Description : Runs blastz between an alternative and multiple reference
                sequences.
  Return type : none
  Exceptions  : thrown if blastz fails
  Caller      : general

=cut

sub run_blastz {
    my ($self, $A_basename, $R_basename) = @_;

    my $tmpdir = $self->tempdir;
    my $bindir = $self->bindir;
    my $id = $self->id;

    my $blastz_cmd = qq($bindir/blastz $tmpdir/$A_basename.fa $tmpdir/$R_basename.fa Q=blastz_matrix.txt T=0 L=10000 H=2200 Y=3400 > $tmpdir/blastz.$id.lav);
    
    unless (-e "$tmpdir/blastz.$id.lav") {
      system($blastz_cmd) == 0 or
        $self->support->log_error("Can't run blastz: $!\n");
    }
}

=head2 lav_to_axt

  Example     : $aligner->lav_to_axt;
  Description : Converts blastz output from lav to axt format. Target and query
                sequences must be present in nib format in the temporary
                directory for this to work.
  Return type : none
  Exceptions  : thrown if lavToAxt fails
  Caller      : general

=cut

sub lav_to_axt {
    my $self = shift;

    my $tmpdir = $self->tempdir;
    my $id = $self->id;

    unless (-e "$tmpdir/blastz.$id.axt") {
      system($self->bindir."/lavToAxt $tmpdir/blastz.$id.lav $tmpdir $tmpdir $tmpdir/blastz.$id.axt") == 0 or $self->support->log_error("Can't run lavToAxt: $!\n");
    }
}

=head2 find_best_alignment

  Example     : $aligner->find_best_alignment;
  Description : Finds the best set of non-overlapping alignments by running
                axtBest.
  Return type : none
  Exceptions  : thrown if axtBest fails
  Caller      : general

=cut

sub find_best_alignment {
    my $self = shift;

    my $tmpdir = $self->tempdir;
    my $id = $self->id;

    unless (-e "$tmpdir/blastz.$id.best.axt") {
      system($self->bindir."/axtBest $tmpdir/blastz.$id.axt all $tmpdir/blastz.$id.best.axt") == 0 or $self->support->log_error("Can't run axtBest: $!\n");
    }
}

=head2 parse_blastz_output

  Example     : $aligner->parse_blastz_output;
  Description : Reads a blastz alignment result from an axt file and creates
                a datastructure containing ungapped alignments from it. Note
                that mismatches are allowed, but separate stats will be
                collected for them.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub parse_blastz_output {
    my $self = shift;

    # read file
    my $tmpdir = $self->tempdir;
    my $id = $self->id;
    my $fh = $self->support->filehandle('<', "$tmpdir/blastz.$id.best.axt");

    # initialize stats
    $self->init_stats(qw(match mismatch gap alignments bp));
    
    my $i = 1;
    my ($header, $A_seq, $R_seq);

    while (my $line = <$fh>) {
        # there are blocks of 4 lines, where line 1 is the header, line 2 is
        # A_seq, line3 is R_seq
        $header = $line unless (($i-1) % 4);
        $A_seq = $line unless (($i-2) % 4);
        chomp $A_seq;
        my @A_arr = split(//, $A_seq);
        $R_seq = $line unless (($i-3) % 4);
        chomp $R_seq;
        my @R_arr = split(//, $R_seq);

        # compare sequences letter by letter
        if ($i % 4 == 0) {
            my $match_flag = 0;
            $self->init_stats(qw(A_gap R_gap));
            my %coords;
            @coords{'R_id', 'A_start', 'A_end', 'R_start', 'R_end', 'strand'} =
                (split(/ /, $header))[4, 2, 3, 5, 6, 7];
            $coords{'R_id'} =~ s/ref_seq\.(.*)/$1/;
            ($coords{'strand'} eq '+') ? ($coords{'strand'} = 1) :
                                         ($coords{'strand'} = -1);
            for (my $j = 0; $j < scalar(@A_arr); $j++) {
                # gap
                if ($A_arr[$j] eq '-' or $R_arr[$j] eq '-') {
                    $self->stats_incr('gap', 1);
                    $self->stats_incr('A_gap', 1) if ($A_arr[$j] eq '-');
                    $self->stats_incr('R_gap', 1) if ($R_arr[$j] eq '-');
                    $match_flag = 0;
                } else {
                    $self->found_match($match_flag, $j, \%coords);
                    $match_flag = 1;

                    # match
                    if ($A_arr[$j] eq $R_arr[$j]) {
                        $self->stats_incr('match', 1);
                    # mismatch
                    } else {
                        $self->stats_incr('mismatch', 1);
                    }
                }
            }
            $self->stats_incr('bp', scalar(@A_arr));
            $self->stats_incr('alignments', 1);
        }
        
        $i++;
    }
}

=head2 cleanup_tmpfiles

  Arg[1-N]    : (optional) list @files - additional tmp files to delete
  Example     : $self->cleanup_tmpfiles('e_seq.fa', 'v_seq.fa');
  Description : deletes temporary files
  Return type : none
  Exceptions  : Warning if file cannot be deleted
  Caller      : general
  Status      : stable

=cut

sub cleanup_tmpfiles {
  my $self = shift;
  my @files = @_;
  
  my $tmpdir = $self->tempdir;
  my $id = $self->id;

  push @files,
    "blastz.$id.lav",
    "blastz.$id.axt",
    "blastz.$id.best.axt";

  foreach my $file (@files) {
    unlink("$tmpdir/$file") or $self->support->log_warning("Couldn't delete file $file: $!");
  }
}

=head2 found_match

  Arg[1]      : Boolean $match_flag - flag indicating if last bp was a match
  Arg[2]      : Int $j - current bp position in the alignment
  Arg[3]      : Hashref $coords - alignment coordinates and strand from blastz
                output
  Description : Populates a datastructure describing blocks of alignment.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_match {
    my ($self, $match_flag, $j, $coords) = @_;

    my $id = $self->id;
    my $align = $self->get_stats('alignments');
    my $R_chr = $self->seq_region_name;

    # last position was a match
    if ($match_flag) {
        
        # adjust align block end
        if ($self->{'_match'}->{$R_chr}->{$id}->[$align]) {
            my $c = scalar(@{ $self->{'_match'}->{$R_chr}->{$id}->[$align] }) - 1;
            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[1] =
                $coords->{'A_start'} + $j - $self->get_stats('A_gap');
            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3] =
                $coords->{'R_start'} + $j - $self->get_stats('R_gap');
        }
    
    # last position was a non-match
    } else {
    
        # start a new align block
        push @{ $self->{'_match'}->{$R_chr}->{$id}->[$align] }, [
            $coords->{'A_start'} + $j - $self->get_stats('A_gap'),
            $coords->{'A_start'} + $j - $self->get_stats('A_gap'),
            $coords->{'R_start'} + $j - $self->get_stats('R_gap'),
            $coords->{'R_start'} + $j - $self->get_stats('R_gap'),
            $coords->{'strand'},
            $coords->{'R_id'},
        ];
    }
}


=head2 adjust_coords

  Arg[1]      : Int $A_start - start of alternatvie block in chromosomal coords
  Arg[2]      : Int $A_end - end of alternative block in chromosomal coords
  Arg[3]      : Arrayref $R_coords - list of start/end pairs of reference
                blocks in chromosomal coords
  Example     : my $R_coords = [ [1, 1000], [3000, 5000] ];
                $aligner->adjust_coords(1, 30000000, $R_coords);
  Description : Adjusts coordinates of blastz alignments to chromosomal coords.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub adjust_coords {
    my ($self, $A_start, $A_end, $R_coords) = @_;
    my $R_chr = $self->seq_region_name;
    my $id = $self->id;

    for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$R_chr}->{$id} || [] }); $align++) {
        for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$R_chr}->{$id}->[$align] || []}); $c++) {
            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[0] += $A_start - 1;
            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[1] += $A_start - 1;

            # forward strand match
            if ($self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[4] == 1) {
                $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[2] += $R_coords->{$self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[5]}->[0] - 1;
                $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3] += $R_coords->{$self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[5]}->[0] - 1;
            
            # reverse strand match
            } else {
                my $tmp_start = $R_coords->{$self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[5]}->[1] - $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3] + 1;

                $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3] = $R_coords->{$self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[5]}->[1] - $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[2] + 1;

                $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[2] = $tmp_start;
            }

            # sanity check: aligned region pairs must have same length
            my $A_len = $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[1] - $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[0];
            my $R_len = $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3] - $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[2];
            $self->support->log_warning("Length mismatch: $A_len <> $R_len in block $id, alignment $align, stretch $c\n", 2) unless ($A_len == $R_len);
        }
    }
}

=head2 filter_overlaps

  Description : DEPRECATED. This filtering algorithm didn't work well. Please
                run the separate script fix_overlaps.pl.

=cut

sub filter_overlaps {
    my $self = shift;
    my $id = $self->id;
    
    foreach my $R_chr (sort keys %{ $self->{'_match'} }) {
        my $coord_check = [];
        # rearrange the datastructure so that we can find overlaps
        foreach my $id (keys %{ $self->{'_match'}->{$R_chr} }) {
            for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$R_chr}->{$id} }); $align++) {
                for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$R_chr}->{$id}->[$align] || []}); $c++) {
                    push @{ $coord_check }, [
                        $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[0],
                        $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[1],
                        $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[2],
                        $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3],
                        $id,
                        $align,
                        $c,
                    ];
                }
            }
        }

        my @A_sort = sort { $a->[0] <=> $b->[0] } @{ $coord_check };
        my @R_sort = sort { $a->[2] <=> $b->[2] } @{ $coord_check };

        # sanity check: alternative alignments must not overlap (axtBest should
        # guarantee that)
        my $last;
        foreach my $c (@A_sort) {
            $self->support->log_warning("Overlapping alternative alignment at ".join(':', $R_chr, $c->[0], $c->[1])." (last_end ".$last->[1].")\n", 1) if ($last and $c->[0] <= $last->[1]);
            $last = $c;
        }

        # now filter reference overlaps
        my @seen;
        $last = undef;
        foreach my $c (@R_sort) {
            if ($last and $c->[2] <= $last->[3]) {
                $self->support->log_verbose("Overlapping reference alignment at ".join(':', $R_chr, $c->[2], $c->[3])." (last_end ".$last->[3].")\n", 1);

                # if last alignment was longer, delete this one
                if ($last->[3]-$last->[2] > $c->[3]-$c->[2]) {
                    undef $self->{'_match'}->{$R_chr}->{$c->[4]}->[$c->[5]]->[$c->[6]];

                # if last alignment was shorter, trace back and delete all
                # overlapping shorter alignments
                } else {
                    foreach my $s (@seen) {
                        # earlier alignment still overlapping
                        if ($c->[2] <= $s->[3]) {
                            # earlier alignment shorter -> delete it
                            if ($s->[3]-$s->[2] < $c->[3]-$c->[2]) {
                                undef $self->{'_match'}->{$R_chr}->{$s->[4]}->[$s->[5]]->[$s->[6]];

                            # this alignment shorter -> delete it
                            } else {
                                undef $self->{'_match'}->{$R_chr}->{$c->[4]}->[$c->[5]]->[$c->[6]];
                                $last = $s;
                                last;
                            }
                        } else {
                            $last = $s;
                            last;
                        }
                    }                
                    
                    $last = $c;
                }
            }
            unshift @seen, $c;
            $last = $c unless ($last);
        }
    }
}

=head2 write_assembly

  Arg[1]      : Bio::*::DBSQL::DBAdaptor $R_dba - reference database adaptor
  Example     : my $R_dba = $support->get_database('ensembl');
                $aligner->write_assembly($R_dba);
  Description : Writes assembly entries for blastz alignments to the database.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub write_assembly {
    my ($self, $R_dba) = @_;

    my $R_dbh = $R_dba->dbc->db_handle;
    my $R_sa = $R_dba->get_SliceAdaptor;

    my $sth = $R_dbh->prepare(qq(
        INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
            asm_start,asm_end, cmp_start, cmp_end, ori)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ));

    $self->support->log("Adding assembly entries for alignments...\n");

    my $i;
    foreach my $R_chr (sort _by_chr_num keys %{ $self->{'_match'} }) {
        # get seq_region_id for alternative and reference chromosome
        my $A_chr = $R_chr;
        my $R_slice = $R_sa->fetch_by_region('toplevel', $R_chr, undef, undef, undef, $self->support->param('assembly'));
        my $R_sid = $R_sa->get_seq_region_id($R_slice);
        # we need to fetch the alternative slice from the reference db
        # explicitely by coord_system, since toplevel attribute is not set
        # there
        my $cs_name = $R_slice->coord_system_name;
        my $A_sid = $R_sa->get_seq_region_id($R_sa->fetch_by_region($cs_name, $A_chr, undef, undef, undef, $self->support->param('altassembly')));

        foreach my $id (sort { $a <=> $b } keys %{ $self->{'_match'}->{$R_chr} }) {
            for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$R_chr}->{$id} }); $align++) {
                for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$R_chr}->{$id}->[$align] }); $c++) {
                    if ($self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]) {
                        $sth->execute(
                            $R_sid,
                            $A_sid,
                            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[2],
                            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[3],
                            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[0],
                            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[1],
                            $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[4],
                        );
                        $i++;
                    }
                }
            }
        }
    }
    $self->support->log("Done inserting $i entries.\n");
}

=head2 stats_incr

  Arg[1]      : String $code - stats code
  Arg[2]      : Int $incr - number by which stats should be incremented
  Example     : $aligner->stats_incr('total_alignments', 1);
  Description : Increments stats.
  Return type : Int - stat value after increment
  Exceptions  : none
  Caller      : general

=cut

sub stats_incr {
    my ($self, $code, $incr) = @_;
    $self->{'_stats'}->{$code} += $incr;
    return $self->{'_stats'}->{$code};
}

=head2 init_stats

  Arg[1]      : Array @codes - stats codes to initialise
  Example     : $aligner->init_stats('match', 'mismatch', 'alignments');
  Description : Initialises stats (i.e. sets to 0).
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub init_stats {
    my ($self, @codes) = @_;
    foreach my $code (@codes) {
        $self->{'_stats'}->{$code} = 0;
    }
}

=head2 get_stats

  Arg[1]      : String $code - stats code
  Example     : my $num_mismatches = $aligner->get_stats('mismatch');
  Description : Stats getter.
  Return type : Int - stats value
  Exceptions  : none
  Caller      : general

=cut

sub get_stats {
    my ($self, $code) = @_;
    return $self->{'_stats'}->{$code};
}

=head2 log_block_stats

  Arg[1]      : Int $indent - indentation level
  Example     : $aligner->log_block_stats(3);
  Description : Logs stats for an alignment block.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub log_block_stats {
    my ($self, $indent) = @_;

    $self->support->log("Blastz alignment stats:\n", $indent);
    $self->support->log(sprintf(FMT2, "Alignments:", $self->get_stats('alignments')), $indent+1);
    if ($self->get_stats('alignments')) {
        $self->support->log(sprintf(FMT1,
            "Matches:",
            $self->get_stats('match'),
            $self->get_stats('match')/$self->get_stats('bp')*100),
        $indent+1);
        $self->support->log(sprintf(FMT1,
            "Mismatches:",
            $self->get_stats('mismatch'),
            $self->get_stats('mismatch')/$self->get_stats('bp')*100),
        $indent+1);
        $self->support->log(sprintf(FMT1,
            "Gaps:",
            $self->get_stats('gap'),
            $self->get_stats('gap')/$self->get_stats('bp')*100),
        $indent+1);
    }
    map { $self->stats_incr($_.'_total', $self->get_stats($_)) }
        qw(match mismatch gap bp);
}

=head2 log_overall_stats

  Example     : $aligner->log_overall_stats;
  Description : Logs overall alignment stats.
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub log_overall_stats {
    my $self = shift;
    
    # blastz
    $self->support->log("\nOverall blastz alignment stats:\n");

    unless ($self->get_stats('alignments')) {
      $self->support->log("No alignments found.\n", 1);
      return;
    }
    
    $self->support->log(sprintf(FMT1,
        "Matches:",
        $self->get_stats('match_total'),
        $self->get_stats('match_total')/$self->get_stats('bp_total')*100),
    1);
    $self->support->log(sprintf(FMT1,
        "Mismatches:",
        $self->get_stats('mismatch_total'),
        $self->get_stats('mismatch_total')/$self->get_stats('bp_total')*100),
    1);
    $self->support->log(sprintf(FMT1,
        "Gaps:",
        $self->get_stats('gap_total'),
        $self->get_stats('gap_total')/$self->get_stats('bp_total')*100),
    1);

    # alignments to be written to assembly table
    $self->support->log_verbose("\nAlignments that will be written to assembly table:\n");
    $self->support->log_verbose(sprintf(FMT3,
        qw(CHR BLOCK ALIGNMENT ALT_START ALT_END REF_START REF_END)),
    1);
    $self->support->log_verbose(('-'x63)."\n", 1);
    foreach my $R_chr (sort _by_chr_num keys %{ $self->{'_match'} }) {
        foreach my $id (sort { $a <=> $b } keys %{ $self->{'_match'}->{$R_chr} }) {
            for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$R_chr}->{$id} }); $align++) {
                for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$R_chr}->{$id}->[$align] }); $c++) {
                    if ($self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]) {
                        $self->support->log_verbose(sprintf(FMT4,
                            $R_chr,
                            $id,
                            $align+1,
                            @{ $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c] }),
                        1);
                        my $l = $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[1] - $self->{'_match'}->{$R_chr}->{$id}->[$align]->[$c]->[0];
                        $self->stats_incr('alignments_total', 1);
                        $self->stats_incr('short1_10_total', 1) if ($l < 11);
                        $self->stats_incr('short11_100_total', 1) if ($l > 10 and $l < 101);
                    }
                }
            }
        }
    }
    $self->support->log("\nAssembly entry stats:\n");
    $self->support->log(sprintf(FMT2,
        "Total alignment blocks:",
        $self->get_stats('alignments_total')),
    1);
    $self->support->log(sprintf(FMT2,
        "Alignments 1-10 bp:",
        $self->get_stats('short1_10_total')),
    1);
    $self->support->log(sprintf(FMT2,
        "Alignments 11-100 bp:",
        $self->get_stats('short11_100_total')),
    1);
}

=head2 _by_chr_num

  Example     : my @sorted = sort _by_chr_num qw(X, 6-COX, 14, 7);
  Description : Subroutine to use in sort for sorting chromosomes. Sorts
                numerically, then alphabetically
  Return type : values to be used by sort
  Exceptions  : none
  Caller      : internal

=cut

sub _by_chr_num {
    my @awords = split /-/, $a;
    my @bwords = split /-/, $b;

    my $anum = $awords[0];
    my $bnum = $bwords[0];

    if ($anum !~ /^[0-9]*$/) {
        if ($bnum !~ /^[0-9]*$/) {
            return $anum cmp $bnum;
        } else {
            return 1;
        }
    }
    if ($bnum !~ /^[0-9]*$/) {
        return -1;
    }

    if ($anum <=> $bnum) {
        return $anum <=> $bnum;
    } else {
        if ($#awords == 0) {
            return -1;
        } elsif ($#bwords == 0) {
            return 1;
        } else {
            return $awords[1] cmp $bwords[1];
        }
    }
}

=head2 AUTOLOAD

  Arg[1]      : (optional) String/Object - attribute to set
  Example     : # setting a attribute
                $self->attr($val);
                # getting the attribute
                $self->attr;
                # undefining an attribute
                $self->attr(undef);
  Description : lazy function generator for getters/setters
  Return type : String/Object
  Exceptions  : none
  Caller      : general

=cut

sub AUTOLOAD {
    my $self = shift;
    my $attr = our $AUTOLOAD;
    $attr =~ s/.*:://;
    return unless $attr =~ /[^A-Z]/;
    no strict 'refs';
    *{$AUTOLOAD} = sub {
        $_[0]->{'_data'}->{$attr} = $_[1] if (@_ > 1);
        return $_[0]->{'_data'}->{$attr};
    };
    $self->{'_data'}->{$attr} = shift if (@_);
    return $self->{'_data'}->{$attr};
}

1;

