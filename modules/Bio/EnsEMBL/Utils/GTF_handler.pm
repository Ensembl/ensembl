#
# BioPerl module for Bio::EnsEMBL::Utils::GTF_handler
#
# Cared for by Elia Stupka <elia@sanger.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::GTF_handler - Utilities for loading and dumping genes
                                   from/to gtf files.

=head1 SYNOPSIS

    GTF_handler->new();

    # To parse a GTF file, and create genes
    my @genes = GTF_handler->parse_file($file);

    # To dump genes in GTF format
    GTF_handler->dump_genes($file,@genes);

=head1 DESCRIPTION

This module is a handler for the GTF format, i.e. a GTF parser and dumper.

Dumping genes is done simply by calling the dump_genes method and passing an
array of genes and a filhandle.

Parsing GTF files is done in three steps, first parsing a file, then mapping
golden path coordinates to ensembl coordinates, and finally writing genes to a database.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Utils::GTF_handler;

use strict;
use vars qw(@ISA);
use Bio::Root::RootI;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;


@ISA = qw(Bio::Root::RootI);

=head2 new

 Title   : new
 Usage   : GTF_handler->new()
 Function: Constructor
 Example : 
 Returns : Reference to an object
 Args    : 

=cut

sub new {
    my ($class, @args) = @_;
    my $self = bless {}, $class;
    $self->reset_exon_phases(1);        # not used, set to zero
    return $self;
}

=head2 seq_name_parser

 Title   : seq_name_parser
 Usage   : GTF_handler->seq_name_parser(CODE)
 Function: Get or set a custom parser for the sequence name field.
           The value returned from the subroutine will be used as
           the contig_id in the Ensembl genes created by parse_file.
           This subroutine could do complicated lookups in an Ensembl
           db to get the correct contig ID.
           The default is to keep the sequence name as is.
 Example : GTF_handler->seq_name_parser(\&my_parser)
 Returns : A ref to a subroutine, or undef
 Args    : Ref to a subrotine

=cut

sub seq_name_parser {
    my( $self, $parser ) = @_;
    
    if ($parser) {
        $self->throw("Not a code reference '$parser'")
            unless ref($parser) eq 'CODE';
        $self->{'_seq_name_parser'} = $parser;
    }
    return $self->{'_seq_name_parser'};
}

=head2 reset_exon_phases

 Title   : reset_exon_phases
 Usage   : GTF_handler->reset_exon_phases(BOOLEAN)
 Function: Flags whether to reset the phases of exons
           given in the CDS lines in the GTF file, or
           to trust them.
           Defaults to TRUE.
 Example : GTF_handler->reset_exon_phases(0); # Keep the values in the GTF file
 Returns : TRUE or FALSE
 Args    : TRUE or FALSE

=cut

# not used, all set to 0
sub reset_exon_phases {
    my( $self, $flag ) = @_;
    
    if (defined $flag) {
        $self->{'_reset_exon_phases'} = $flag;
    }
    return $self->{'_reset_exon_phases'};
}

=head2 parse_file

 Title   : parse_file
 Usage   : GTF_handler->parse_file(filename)
 Function: Parses a GTF file, reading in features
 Example : GTF_handler->parse_file(gens.gtf)
 Returns : array of Bio::EnsEMBL::Gene objects
 Args    : name of the file

=cut

sub parse_file {
    my( $self, $fh ) = @_;
    
    $self->warn("Do not use for dumping or translating; all phases are set to 0 !\n");

    $self->throw("No filehandle supplied") unless $fh;
    
    my %valid_gtf_feature = map {$_, 1} qw{ exon cds start_codon stop_codon };
    my( %gtf,           # Big hash that stores an entire parsed GTF file.
        %gene_type,     # Just maps gene names to types.
        );
    my $seq_name_parser = $self->seq_name_parser;
    
  GTF_LINE: while (<$fh>) {
        next if /^#/;
        next if /^\s*$/;
        chomp;
        
        my $gtf_line = $_;
        my ($seq_name,    $source, $feature,
            $start,  $end,    $score,
            $strand, $phase,  $group_field) = split /\s+/, $gtf_line, 9;
        
        $feature = lc $feature;
        $seq_name = &$seq_name_parser($seq_name) if $seq_name_parser;
        
        # Convert strand to Ensembl convention
        if ($strand eq '+') {
            $strand = 1;
        }
        elsif ($strand eq '-') {
            $strand = -1;
        }
        else {
            $strand = 0;
        }
        
        # It is isn't a legal GTF file without a group field
        unless ($group_field) {
            warn("Skipping unparseable line : '$gtf_line'\n");
            next;
        }

        # Extract the extra information from the final field of the GTF line.
        my ($gene_name, $gene_id, $transcript_id, $exon_num, $exon_id) =
            $self->parse_group_field($group_field);

        unless ($gene_id) {
            warn("Skipping line with no gene_id: '$_'\n");
	    next GTF_LINE;
        }
        unless ($transcript_id) {
            warn("Skipping line with no transcript_id: '$_'\n");
	    next GTF_LINE;
        }
        
        # Put information from the feature into the appropriate
        # place in the monster GTF hash.
        if ($valid_gtf_feature{$feature}) {
            
            unless (defined($exon_num)) { # typically for start/stop_codons
                #warn("Skipping line with no exon_number: '$gtf_line'\n");
                next GTF_LINE;
            }
            
            # Create an entry in the %gtf hash if isn't there already
            $gtf{$gene_id}{$transcript_id}{$exon_num} ||= {};
            my $exon = $gtf{$gene_id}{$transcript_id}{$exon_num};
            
            if ($feature    eq 'exon'
                or $feature eq 'stop_codon'
                or $feature eq 'start_codon'
                or $feature eq 'cds') {
                $exon->{$feature} = [$seq_name, $start, $end, $strand, $phase];
            }

            # Record the exon_id if we have it
            if ($exon_id and defined($exon_num)) {
                $exon->{'exon_id'} = $exon_id;
            }
        }
        else {
            #warn "Ignoring '$feature' feature: '$gtf_line'\n";
            next GTF_LINE;
        }
        
        # Record the source that reported this gene
        $gene_type{$gene_id} = $source if $source;

        # keep track of igi here? 
    }
    
    # Transform the ugly, great, multi-level hash we've
    # created into nice Ensembl gene objects.
    my( @genes );
    foreach my $gene_id (keys %gtf) {
        push(@genes, $self->_make_gene_from_gtf_hash($gene_id, 
                                                     $gene_type{$gene_id}, 
                                                     $gtf{$gene_id}));
    }
    return @genes;
}                                       # parse_file

sub parse_group_field {
    my( $self, $group_field ) = @_;
    
    my ($gene_name, $gene_id, $transcript_id, $exon_num, $exon_id);

    # Parse the group field
    foreach my $tag_val (split /;/, $group_field) {

        # Trim trailing and leading spaces
        $tag_val =~ s/^\s+|\s+$//g;

        my($tag, $value) = split /\s+/, $tag_val, 2;

        # Remove quotes from the value
        $value =~ s/^"|"$//g;
        $tag = lc $tag;

        if ($tag eq 'gene_name') {
            $gene_name = $value;
        }
        elsif ($tag eq 'gene_id') {
            $gene_id = $value;
        }
        elsif ($tag eq 'transcript_id') {
            $transcript_id = $value;
        }
        elsif ($tag eq 'exon_number') {
            $exon_num = $value;
        }
        elsif ($tag eq 'exon_id') {
            $exon_id = $value;
        }
        else {
            #warn "Ignoring group field element: '$tag_val'\n";
        }
    }
    
    return($gene_name, $gene_id, $transcript_id, $exon_num, $exon_id);
}                                       # parse_group_field

sub _make_gene_from_gtf_hash {
    my( $self, $gene_id, $source, $gene_gtf ) = @_;
    
    # Make a new gene object
    my $gene = Bio::EnsEMBL::Gene->new();
    
    # NO IDEA WHAT TO DO WITH STABLE IDs
    #$gene->id($gene_id);
    $gene->type($source);
    
    # The %gene_exons hash is used so each exon is only made once,
    # even if it occurs in multiple transcripts.  The keys of the
    # hash are "<seqname>-<start>-<end>-<strand>".
    my( %gene_exons );
    foreach my $trans_id (keys %$gene_gtf) {
        my $trans_gtf = $gene_gtf->{$trans_id};
        eval {
            $self->_make_transcript($gene, $trans_id, $trans_gtf, \%gene_exons);
        };
        warn $@ if $@;
    }
    $self->throw("Failed to make any transcripts")
      unless $gene->each_Transcript;
    return $gene;
}                                       # _make_gene_from_gtf_hash

sub _make_transcript {
    my( $self, $gene, $trans_id, $trans_gtf, $gene_exons ) = @_;
    
    # Make a new transcript object
    my $transcript = Bio::EnsEMBL::Transcript->new;

    # NO IDEA WHAT TO DO WITH STABLE IDS
    #$transcript->id($trans_id);

    my @exon_num = sort {$a <=> $b} keys %$trans_gtf;
    my $translation_start = []; # For recording the exon number and sequence position
    my $translation_end   = []; # of the start and end of the translation.
    foreach my $n (@exon_num) {

        # Try to make a new exon
        my( $exon, $t_start, $t_end );
        eval {
            ($exon, $t_start, $t_end) =
                $self->_make_exon($gene, $n, $trans_gtf->{$n});
        };
        if ($@) {
            # Provide a nice exception if something is wrong with the data
            my $exon_error_string = "$@\nError making exon in transcript '$trans_id' from data:\n";
            $exon_error_string .= $self->_pretty_exon_data_string($n, $trans_gtf->{$n});
            $self->throw($exon_error_string);
        }

        # Use an existing exon instead if we already have it
        my $exon_hash_key = join('-',
            $exon->contig_id,
            $exon->start,
            $exon->end,
            $exon->strand,
            );
        if ($gene_exons->{$exon_hash_key}) {
            $exon = $gene_exons->{$exon_hash_key};
        } else {
            $gene_exons->{$exon_hash_key} = $exon;
        }
        
        # Record translation starts and ends if they are in this exon
        if ($t_start) {
            $translation_start = [$t_start, $exon];
        }
        if ($t_end) {
            $translation_end = [$t_end, $exon];
        }

        # Add this exon to the transcript
        $transcript->add_Exon($exon);
    }
    
    $self->_make_translation($transcript, $translation_start, $translation_end);

    ## all phases are 0, now. 
    # $self->_correct_exon_phases($transcript)
    #    if $self->reset_exon_phases;
    
    $gene->add_Transcript($transcript);
}                                       # _make_transcript

# returns a new Exon, and translation start/stop if it can find it.
sub _make_exon {
    my( $self, $gene, $exon_num, $exon_gtf ) = @_;
        
    my( $seq_name, $start, $end, $strand, $phase );

    # Could check that the sequence names are all the same
    foreach my $feature (qw{ exon cds }) {
        last if $seq_name = $exon_gtf->{$feature}[0];
    }
    die "No sequence name" unless $seq_name;
    
    foreach my $feature (qw{ exon cds }) {
        last if $start = $exon_gtf->{$feature}[1];
    }
    die "No exon start" unless $start;

    foreach my $feature (qw{ exon cds }) {
        last if $end = $exon_gtf->{$feature}[2];
    }
    die "No exon end" unless $end;
    
    # Could check that all the strands are the same here
    foreach my $feature (qw{ exon cds }) {
        last if $strand = $exon_gtf->{$feature}[3];
    }
    die "No strand" unless $strand;

    $phase = 0; # ignore it altogether!
#    $phase = $exon_gtf->{'cds'}[4] if $exon_gtf->{'cds'};
    
    # Make an exon object
    my $exon = Bio::EnsEMBL::Exon->new($start, $end, $strand);


    $exon->contig_id($seq_name);
    $exon->seqname($seq_name);  ### What is seqname for ???
    $exon->phase($phase);
    $exon->sticky_rank(1);

    
    my $exon_id = $exon_gtf->{'exon_id'};
    # Make a, hopefull unique, exon ID, it there wasn't one in the file
    unless ($exon_id) {
        my $gene_id = "$gene";
	$gene_id = s/HASH//g;
        my $type    = $gene->type;
        $exon_id = "$type-$gene_id-$exon_num";
    }
    # No idea what to do with stable ids
    #$exon->id($exon_id);
    
    # Fill in the translation start and ends if we have them
    my( $t_start, $t_end );
    if ($exon_gtf->{'start_codon'}) {
        if ($strand == 1) {
            $t_start = $exon_gtf->{'start_codon'}[1];
        } else {
            # Strand is -1
            $t_start = $exon_gtf->{'start_codon'}[2];
        }
    }
    
    if ($exon_gtf->{'stop_codon'}) {
        if ($strand == 1) {
            $t_end = $exon_gtf->{'stop_codon'}[1] - 1;
        } else {
            # Strand is -1
            $t_end = $exon_gtf->{'stop_codon'}[2] + 1;
        }
    }

    ## HA! translation start/end is wrong; EnsEMBL code expects it 
    ## in exon coordinates, not in the fpc coords we're given:
    
    
    return($exon, $t_start, $t_end);
}                                       # _make_exon

sub _pretty_exon_data_string {
    my( $self, $exon_num, $exon_gtf ) = @_;
    
    my $pretty_string = "$exon_num => {\n";
    foreach my $key (sort keys %$exon_gtf) {
        my $field = $exon_gtf->{$key};
        if (ref($field) eq 'ARRAY') {
            $pretty_string .= "  $key => [";
            $pretty_string .= join(', ', map "'$_'", @{$field});
            $pretty_string .= "]\n";
        } else {
            $pretty_string .= "  $key => '$field'\n";
        }
    }
    $pretty_string .= "}\n";
}

sub _make_translation {
    my( $self, $transcript, $translation_start, $translation_end ) = @_;

    # Set the translation start to the start of the first exon
    # if it isn't explicitly stated
    if (@$translation_start == 0) {
        my $exon = $transcript->start_exon;
        my( $start );
        if ($exon->strand == 1) {
            $start = $exon->start;
        } else {
            $start = $exon->end;
        }
        $translation_start = [$start, $exon];
    }

    # Set the translation end to the end of the last exon
    # if it isn't explicitly stated
    if (@$translation_end == 0) {
        my $exon = $transcript->end_exon;
        my( $end );
        if ($exon->strand == 1) {
            $end = $exon->end;
        } else {
            $end = $exon->start;
        }
        $translation_end = [$end, $exon];
    }

    ### these start/ends are in fpc contig coords; we now have to map it
    ### back to exon coords (the new convention); 
    {   my ($start, $exonid) = @$translation_start;
        $start -= $transcript->start_exon->start; 
        $start +=1;
        $translation_start = [ $start, $exonid ];
    }
    # and same for end:
    {   my ($end, $exonid) = @$translation_end;
        $end -= $transcript->end_exon->start; 
        $end +=1;
        $translation_end = [ $end, $exonid ] ;
    }

    # New translation object
    my $translation = Bio::EnsEMBL::Translation->new;
#    my $tl_id = $transcript->id;
#    $tl_id =~ s/ENST/ENSQ/g;            # invent an id?
    ### note: re-using transcript_ids for translation_id's, for simplicity!

    # NO IDEA WHAT TO DO with STABLE IDs
    #$translation->id($transcript->id);
    
    # Translation start
    $translation->start        ($translation_start->[0]);
    $translation->start_exon($translation_start->[1]);

    # Translation end
    $translation->end        ($translation_end->[0]);
    $translation->end_exon($translation_end->[1]);
    
    $transcript->translation($translation);
}                                       # _make_translation


## not used anymore ...
sub _correct_exon_phases {
    my( $self, $transcript ) = @_;
    
    my $trans_id = $transcript->id;
    my $translation = $transcript->translation;
    my $t_start       = $translation->start;
    my $start_exon_id = $translation->start_exon_id;
    my $t_end         = $translation->end;
    my $end_exon_id   = $translation->end_exon_id;
    
    my( $phase );
    foreach my $exon ($transcript->each_Exon) {
        my $exon_id    = $exon->id;
        my $exon_phase = $exon->phase;
        if ($exon_id eq $start_exon_id) {
            # First coding exon
            if (defined($exon_phase) and $exon_phase != 0) {
                #warn "Resetting phase for first coding exon '$exon_id' from '$exon_phase' to '0'\n";
            }
            $phase = 0;
            $exon->phase($phase);
            if ($exon->strand == 1) {
                $phase = $self->calculate_end_phase($phase, $t_start, $exon->end, $exon->strand);
            } else {
                $phase = $self->calculate_end_phase($phase, $exon->start, $t_start, $exon->strand);
            }
        }
        elsif (defined $phase) {
            if (defined($exon_phase) and $exon_phase != $phase) {
                #warn "Resetting phase for exon '$exon_id' from '$exon_phase' to '$phase'\n";
            }
            $exon->phase($phase);
            $phase = $self->calculate_end_phase($phase, $exon->start, $exon->end, $exon->strand);
        }
        else {
            # We haven't got to the first coding exon yet
            next;
        }
        
        if ($exon->id eq $end_exon_id) {
            # Check that we end in frame
            if ($exon->strand == 1) {
                $phase = $self->calculate_end_phase($phase, $exon->start, $t_end, $exon->strand);
            } else {
                $phase = $self->calculate_end_phase($phase, $t_start, $exon->end, $exon->strand);
            }
            #warn "Transcript '$trans_id' ends in phase '$phase', not '0' in exon '$exon_id'\n"
                #unless $phase == 0;
            
            # We've reached the last coding exon
            last;
        }
    }
}                                       # _correct_exon_phases

# not used anymore, all set to 0
sub calculate_end_phase {
    my( $self, $phase, $start, $end, $strand ) = @_;
    
    if ($strand == -1) {
        ($start, $end) = ($end, $start);
    }
    my $length = $end - $start + 1;
    return ($phase + $length) % 3;
}

=head2 dump_genes

 Title   : dump_genes
 Usage   : GTF_handler->dump_genes($file,@genes)
 Function: Dumps genes to a GTF file
 Example : 
 Returns : 
 Args    : array of Bio::EnsEMBL::Gene objects

=cut

sub dump_genes {
    my ($self, $fh, @genes) = @_;
    
    print $fh 
        "##gff-version 2\n",
        "##source-version EnsEMBL-Gff 1.0\n",
        "##date ".gff_date()."\n";

  GENE: foreach my $gene (@genes) {
        my $gene_id = $gene->id;
        print STDERR "Dumping gene '$gene_id'\n";
        
        my $type = $gene->type || 'ensembl';
	
      TRANSCRIPT: foreach my $trans ($gene->each_Transcript) {
            my $trans_id = $trans->id;
	    print STDERR "Dumping transcript '$trans_id'\n";
            
            # Extract the data needed from the translation object
            my $skipping_string = "Skipping transcript '$trans' : no translation";
            
            my $translation = $trans->translation
                or $self->warn("$skipping_string object")
                and next TRANSCRIPT;
	    
            my $translation_start = $translation->start
		or $self->warn("$skipping_string start")
		and next TRANSCRIPT;

            #tania's fix to use the stable_id becaure start_exon_id is deprecated
            my $start_exon_id = $translation->start_exon->stable_id
		or $self->warn("$skipping_string start_exon_id")
		and next TRANSCRIPT;

            my $translation_end = $translation->end
		or $self->warn("$skipping_string end")
		and next TRANSCRIPT;

            #tania's fix to use the stable_id becaure end_exon is deprecated
            my $end_exon_id = $translation->end_exon->stable_id
		or $self->warn("$skipping_string end_exon_id")
		and next TRANSCRIPT;

            my @exons = $trans->each_Exon;
            
            # Find the start and end exons
            my( $start_exon, $end_exon );
            foreach my $ex (@exons) {
                
                $start_exon = $ex if ($ex->id eq $start_exon_id);
                $end_exon   = $ex if ($ex->id eq $end_exon_id);
            }

	    my $transcript_string = '';
            my @group_fields = (
                qq{gene_id "$gene_id"}, 
                qq{transcript_id "$trans_id"}
                );
	    #Loop through all exons to find start and end
            my $exon_num = 0;
            my( $seen_start, $seen_end );
	    foreach my $exon (@exons) {
                my $exon_id    = $exon->id;
                my $exon_start = $exon->start;
                my $exon_end   = $exon->end;
                my $seq_name   = $exon->seqname;
                my $phase      = $exon->phase;
                my $score      = $exon->score || 0;
                my $strand    = ($exon->strand == 1) ? '+' : '-';

                # Make the group field for this exon
                $exon_num++;
                my $exon_num_field = qq{exon_number $exon_num};
                my $exon_id_field  = qq{exon_id "$exon_id"};
                my $group = join('; ', (@group_fields, $exon_num_field, $exon_id_field));
                
                # Is the start codon here?
                if ($exon_id eq $start_exon_id) {
                    $seen_start = 1;
                    my( $x, $y ) = ($translation_start + $exon_start-1, 
                                    $translation_start + $exon_start-1 + 2);
                    $transcript_string .= join("\t",
                        $seq_name, $type, 'start_codon',
                        $x, $y, $score, $strand,
                        '.', # phase
                        $group) ."\n";
                }
                
                # Add strings for the exon and CDS lines
                $transcript_string .= join("\t",
                    $seq_name, $type, 'exon',
                    $exon_start, $exon_end, $score, $strand,
                    '.', # exon lines don't have phase
                    $group) ."\n";
                
                $transcript_string .= join("\t",
                    $seq_name, $type, 'CDS',
                    $exon_start, $exon_end, $score, $strand,
                    $phase, # CDS lines do have phases
                    $group) ."\n";
                
                # Is the end codon here?
                if ($exon_id eq $end_exon_id) {
                    $seen_end = 1;
                    my( $x, $y ) = ($translation_end + $exon_start-1 -2, 
                                    $translation_end + $exon_start-1);

                    $transcript_string .= join("\t",
                        $seq_name, $type, 'stop_codon',
                        $x, $y, $score, $strand,
                        '.', # phase
                        $group) ."\n";
                }
	    }
            if ( ! $seen_start ) {
                $self->warn("Could not find start exon for transcript '$trans_id'");
            }
            elsif ( ! $seen_end ) {
                $self->warn("Could not find start end for transcript '$trans_id'");
            }
            else {
                print $fh $transcript_string
            }
	}
    }
}

sub gff_date {
    my $time = shift || time;

    # Get time info
    my ($mday, $mon, $year) = (localtime($time))[3,4,5];

    # Change numbers to double-digit format
    ($mon, $mday) = ('00'..'31')[($mon + 1), $mday];

    # Make year
    $year += 1900;

    return "$year-$mon-$mday";

}

1;

__END__
