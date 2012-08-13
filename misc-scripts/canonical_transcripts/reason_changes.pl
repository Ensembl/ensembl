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
# If it was not one of these we will complain bitterly and you should
# investigate. If it is a known reason after investigation, that is an
# allowed scenario, then code it into this script.
#
##################################################

use strict;
use warnings;
use JSON;

my $file = $ARGV[0];

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
  other => 0
);

sub compare {
  my ($old, $new) = @_;
  my (@old_contents) = @{$old}[1..4];
  my ($old_id, $old_length) = ($old->[-1], $old->[-2]);
  my (@new_contents) = @{$new}[1..4];
  my ($new_id, $new_length) = ($new->[-1], $new->[-2]); 
  
  my $old_key = join(q{'=!='}, @old_contents);
  my $new_key = join(q{'=!='}, @new_contents);
  
  my $reason_key;
  my $possible_reason;
  if($old_key eq $new_key) {
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
    if( join(q{=!=},@{$old}[1..3]) eq join(q{=!=},@{$new}[1..3]) ) {
      $reason_key = 'translation_length';
    } 
    else {
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