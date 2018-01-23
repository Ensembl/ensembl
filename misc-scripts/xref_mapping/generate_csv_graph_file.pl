#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
   
#if any parameter is passed then dump out accumulated data
# else just the nu,berof hits at each cutoff

my $acc;
my $temp = shift;
if(defined($temp)){
  if($temp eq "acc"){
    $acc = $temp
  }
  else{
    die "only acc allowed to set to dump acummulated passes\n";
  }
}


my  @q_bin;
my  @t_bin;
my  @m_bin;

for(my $i=0;$i<=100; $i++){
  $q_bin[$i] = 0;
  $t_bin[$i] = 0;
  $m_bin[$i] = 0;
}

# percentage_vals
while(<>){
  my ($q,$t,$acc) = split;

  $q_bin[$q]++;
  $t_bin[$t]++;
  if($q>$t){
    $m_bin[$q]++;
  }
  else{
    $m_bin[$t]++;
  }
}
if(!defined($acc)){
  for(my $i=0;$i<=100; $i++){
    print $i.",".$q_bin[$i].",".$t_bin[$i].",".$m_bin[$i]."\n";
  }
  exit;
}


my $q_accum=0;
my $t_accum=0;
my $m_accum=0;

for(my $i=100;$i>=0;$i--){
  $q_accum+=$q_bin[$i];
  $t_accum+=$t_bin[$i];
  $m_accum+=$m_bin[$i];
  
  print "$i,$q_accum,$t_accum,$m_accum\n";
}


