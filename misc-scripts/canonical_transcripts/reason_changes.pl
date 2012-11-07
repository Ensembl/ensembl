#!/usr/bin/env perl

###################################################
#
# A script to investigate the results of canonical assignment pre-run log 
# output. We know of 3 normal scenarios for assignment
#
#   - Stable ID ordering
#   - Transcript length (longest)
#   - Translation length (longest)
#
# A not so normal scenario is the choice of an NMD over an Ensembl protein
# coder. To emit these in all their glory, add NMD to the command line. e.g.
#
# perl reason_changes.pl logfile NMD
#
# If it was not one of these we will complain bitterly and you should
# investigate. If it is a known reason after investigation, that is an
# allowed scenario, then code it into this script.
#
##################################################

use strict;
use warnings;
use JSON;

my $file = $ARGV[0];
my $blurt = $ARGV[1] ;

if(! $file) {
  die "Was not given a file. Command is $0 PATH_TO_FILE";
}
if(! -f $file) {
  die "Could not find file";
}

my %reason_count = (
  transcript_length => 0,
  id_ordering => 0,
  translation_length => 0,
  havana_merge_nmd_over_e_coding => 0,
  other => 0
);

sub compare {
  my ($old, $new) = @_;
  my (@old_contents) = @{$old}[1..4];
  my ($old_id, $old_length) = ($old->[-1], $old->[-2]);
  my (@new_contents) = @{$new}[1..4];
  my ($new_id, $new_length) = ($new->[-1], $new->[-2]); 
  
  my $old_key = join(q{'=!='}, @old_contents) if defined($old_contents[0]);
  my $new_key = join(q{'=!='}, @new_contents);
  
  my $reason_key;
  my $possible_reason;
  if($old_key && $old_key eq $new_key) {
    if($new_length > $old_length) {
      $reason_key = 'transcript_length';
    }
    elsif($new_length < $old_length) {
      $reason_key = 'other';
      $possible_reason = 'the old transcript length was longer than the new transcript';
    }
    else {
      if( ($new_id cmp $old_id) > 0 ) {
        $reason_key = 'other';
        $possible_reason = 'the old ID was lexographically lower than the new ID';
      }
      else {
        $reason_key = 'id_ordering';
      }
    }
  }
  else {
    #If all other keys are equal then we used translation length
    if( defined($old->[1]) && join(q{=!=},@{$old}[1..3]) eq join(q{=!=},@{$new}[1..3]) ) {
      $reason_key = 'translation_length';
    } 
    #If we have translatable genes and we prefered Havana merged NMD over E/H PC
    elsif(
      $old->[1] && $new->[1] && #We had translatable genes 
      $old->[2] == 3 &&         #Old was an Ensembl/Havana transcript
      $new->[2] == 2 &&         #New was E! Havana merged 
      $old->[3] == 1 &&         #Old was a Protein Coding
      $new->[3] == 2 &&         #New was a "NMD"
      $blurt !~ /nmd/i          #User suppressing this test     
      ) {
      $reason_key = 'havana_merge_nmd_over_e_coding';
    } elsif (! $old->[1] ) {
        # The old canonical transcript was not found.
        $reason_key = 'new';
    } else {
      $reason_key = 'other';
    }
  }
  
  if($reason_key eq 'other') {
    if($possible_reason) {
      printf("ERROR: %s -> %s change was wrong because %s. Please investigate\n", $old_id, $new_id, $possible_reason);
    }
    else {
      printf("ERROR: %s -> %s change not due to a known reason. Please investigate\n", $old_id, $new_id);
    }
    printf("\tURL: www.ensemblgenomes.org/id/%s\n", $old_id);
    printf("\tURL: www.ensemblgenomes.org/id/%s\n", $new_id);
    printf("\tENC: %s\n", encode_json($old));
    printf("\tENC: %s\n", encode_json($new));
    print("\tSee Bio::EnsEMBL::Utils::TranscriptSelector::sort_into_canonical_order() POD to decode\n");
  }
  
  $reason_count{$reason_key}++;
  
  return;
}

my $old;
my $new;

open my $fh, '<', $file or die "Could not open $file for reading: $!";
while(my $line = <$fh>) {
  chomp $line;
  if($line eq '//' && $old && $new) {
    compare($old, $new);
    ($old, $new) = (undef, undef);
    #Perform comparison
  }
  elsif($line =~ /^Old=(\[.+\])$/) {
    $old = eval $1;
  }
  elsif($line =~ /^New=(\[.+\])$/) {
    $new = eval $1;
  }
}

compare($old, $new) if $old && $new;

foreach my $reason (keys %reason_count) {
  printf('%s = %d'."\n", $reason, $reason_count{$reason});
}