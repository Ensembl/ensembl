#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#
# gtfSerializer.t
#
# Test Bio::EnsEMBL::Utils::IO::GTFSerializer.
#

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Differences;
use IO::String;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO::GTFSerializer;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBEntry;

my $mtdb = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $mtdb->get_DBAdaptor("core");

{
  # Creating a temporary transcript to get a transcript which will have a Selenocysteinea,
  # 5' and 3' UTRs and a CCDS reference
  my $slice = $db->get_SliceAdaptor()->fetch_by_toplevel_location('20');
  my $exon = Bio::EnsEMBL::Exon->new(-START => 30274331, -END => 30274348, -STRAND => 1, -SLICE => $slice, -PHASE => 0, -STABLE_ID => 'e1');
  my $exon_two = Bio::EnsEMBL::Exon->new(-START => 30274401, -END => 30274404, -STRAND => 1, -SLICE => $slice, -PHASE => 0, -STABLE_ID => 'e2');
  my $transcript = Bio::EnsEMBL::Transcript->new(-EXONS => [$exon, $exon_two], -STABLE_ID => 'TRANS', -BIOTYPE => 'protein_coding');
  my $seq_edit = Bio::EnsEMBL::SeqEdit->new(-CODE => '_selenocysteine', -START => 2, -END => 2, -ALT_SEQ => 'U');
  my $translation = Bio::EnsEMBL::Translation->new(-START_EXON => $exon, -END_EXON => $exon, -SEQ_START => 4, -SEQ_END => 15, -STABLE_ID => 'PEP');
  $translation->add_Attributes($seq_edit->get_Attribute());
  $transcript->translation($translation);
  my $gene = Bio::EnsEMBL::Gene->new(-TRANSCRIPTS => [$transcript], -STABLE_ID => 'GENE', -BIOTYPE => 'protein_coding');
  $transcript->add_DBEntry(Bio::EnsEMBL::DBEntry->new(-PRIMARY_ID => 'CCDS.1', -DBNAME => 'CCDS'));
  
  # Stupid transcript code has a cache per DB if using it.
  $transcript->{dbentriesCCDS} = $transcript->{dbentries};
  
  my $fh = IO::String->new();
  my $gtf_serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($fh);
  $gtf_serializer->print_Gene($gene);
  my $gtf = <<GTF;
20\tprotein_coding\tgene\t30274331\t30274404\t.\t+\t.\tgene_id \"GENE\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\";
20\tprotein_coding\ttranscript\t30274331\t30274404\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; tag \"seleno\";
20\tprotein_coding\tSelenocysteine\t30274337\t30274339\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; tag \"seleno\";
20\tprotein_coding\texon\t30274331\t30274348\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; exon_number \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; exon_id \"e1\"; tag \"seleno\";
20\tprotein_coding\tCDS\t30274334\t30274345\t.\t+\t0\tgene_id \"GENE\"; transcript_id \"TRANS\"; exon_number \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; protein_id \"PEP\"; tag \"seleno\";
20\tprotein_coding\texon\t30274401\t30274404\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; exon_number \"2\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; exon_id \"e2\"; tag \"seleno\";
20\tprotein_coding\tUTR\t30274331\t30274333\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; tag \"seleno\";
20\tprotein_coding\tUTR\t30274346\t30274348\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; tag \"seleno\";
20\tprotein_coding\tUTR\t30274401\t30274404\t.\t+\t.\tgene_id \"GENE\"; transcript_id \"TRANS\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; tag \"CCDS\"; ccds_id \"CCDS.1\"; tag \"seleno\";
GTF
  eq_or_diff(${$fh->string_ref}, $gtf, 'Checking custom Gene object dumps UTRs, Selenocysteine, seleno tag and CCDS');
}


my $transcripts_gtf = 
  {
   ENST00000310998 => "#!genome-version NCBI33
20\tprotein_coding\ttranscript\t30274334\t30298904\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\";
20\tprotein_coding\texon\t30274334\t30274425\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"1\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; exon_id \"ENSE00001155821\";
20\tprotein_coding\tCDS\t30274334\t30274425\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"1\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; protein_id \"ENSP00000308980\";
20\tprotein_coding\texon\t30284451\t30284562\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"2\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; exon_id \"ENSE00000859878\";
20\tprotein_coding\tCDS\t30284451\t30284562\t.\t+\t1\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"2\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; protein_id \"ENSP00000308980\";
20\tprotein_coding\texon\t30285597\t30285782\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"3\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; exon_id \"ENSE00000661128\";
20\tprotein_coding\tCDS\t30285597\t30285782\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"3\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; protein_id \"ENSP00000308980\";
20\tprotein_coding\texon\t30295724\t30295792\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"4\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; exon_id \"ENSE00000991632\";
20\tprotein_coding\tCDS\t30295724\t30295792\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"4\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; protein_id \"ENSP00000308980\";
20\tprotein_coding\texon\t30296506\t30296579\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"5\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; exon_id \"ENSE00001155786\";
20\tprotein_coding\tCDS\t30296506\t30296579\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"5\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; protein_id \"ENSP00000308980\";
20\tprotein_coding\texon\t30298823\t30298904\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"6\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; exon_id \"ENSE00001180831\";
20\tprotein_coding\tCDS\t30298823\t30298904\t.\t+\t1\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000310998\"; exon_number \"6\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf125\"; transcript_source \"ensembl\"; protein_id \"ENSP00000308980\";
",

   ENST00000278995 => "#!genome-version NCBI33
20\tprotein_coding\ttranscript\t30285705\t30300924\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\";
20\tprotein_coding\texon\t30285705\t30285782\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"1\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; exon_id \"ENSE00000991635\";
20\tprotein_coding\tCDS\t30285705\t30285782\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"1\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; protein_id \"ENSP00000278995\";
20\tprotein_coding\texon\t30295724\t30295792\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"2\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; exon_id \"ENSE00000991632\";
20\tprotein_coding\tCDS\t30295724\t30295792\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"2\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; protein_id \"ENSP00000278995\";
20\tprotein_coding\texon\t30298823\t30298913\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"3\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; exon_id \"ENSE00000991636\";
20\tprotein_coding\tCDS\t30298823\t30298913\t.\t+\t0\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"3\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; protein_id \"ENSP00000278995\";
20\tprotein_coding\texon\t30300869\t30300924\t.\t+\t.\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"4\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; exon_id \"ENSE00000991637\";
20\tprotein_coding\tCDS\t30300869\t30300924\t.\t+\t2\tgene_id \"ENSG00000131044\"; transcript_id \"ENST00000278995\"; exon_number \"4\"; gene_name \"C20orf125\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9BR18\"; transcript_source \"ensembl\"; protein_id \"ENSP00000278995\";
",

   ENST00000252021 => "#!genome-version NCBI33
20\tprotein_coding\ttranscript\t30301733\t30318881\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\";
20\tprotein_coding\texon\t30301733\t30301887\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; exon_id \"ENSE00001155773\";
20\tprotein_coding\tCDS\t30301733\t30301887\t.\t+\t0\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; protein_id \"ENSP00000252021\";
20\tprotein_coding\texon\t30309589\t30309718\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"2\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; exon_id \"ENSE00001180789\";
20\tprotein_coding\tCDS\t30309589\t30309718\t.\t+\t1\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"2\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; protein_id \"ENSP00000252021\";
20\tprotein_coding\texon\t30310552\t30310748\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"3\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; exon_id \"ENSE00001048784\";
20\tprotein_coding\tCDS\t30310552\t30310748\t.\t+\t0\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"3\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; protein_id \"ENSP00000252021\";
20\tprotein_coding\texon\t30313256\t30313369\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"4\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; exon_id \"ENSE00001048782\";
20\tprotein_coding\tCDS\t30313256\t30313369\t.\t+\t1\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"4\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; protein_id \"ENSP00000252021\";
20\tprotein_coding\texon\t30315002\t30315126\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; exon_id \"ENSE00001155819\";
20\tprotein_coding\tCDS\t30315002\t30315126\t.\t+\t1\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; protein_id \"ENSP00000252021\";
20\tprotein_coding\texon\t30318805\t30318881\t.\t+\t.\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"6\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; exon_id \"ENSE00001180793\";
20\tprotein_coding\tCDS\t30318805\t30318878\t.\t+\t2\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"6\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\"; protein_id \"ENSP00000252021\";
20\tprotein_coding\tstop_codon\t30318879\t30318881\t.\t+\t0\tgene_id \"ENSG00000174873\"; transcript_id \"ENST00000252021\"; exon_number \"6\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_source \"ensembl\";
",

   ENST00000202017 => "#!genome-version NCBI33
20\tprotein_coding\ttranscript\t30320853\t30327869\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\";
20\tprotein_coding\texon\t30327735\t30327869\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"1\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; exon_id \"ENSE00001155739\";
20\tprotein_coding\tCDS\t30327735\t30327869\t.\t-\t0\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"1\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; protein_id \"ENSP00000202017\";
20\tprotein_coding\texon\t30326172\t30326247\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"2\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; exon_id \"ENSE00000661139\";
20\tprotein_coding\tCDS\t30326172\t30326247\t.\t-\t0\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"2\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; protein_id \"ENSP00000202017\";
20\tprotein_coding\texon\t30324668\t30324742\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"3\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; exon_id \"ENSE00000661138\";
20\tprotein_coding\tCDS\t30324668\t30324742\t.\t-\t2\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"3\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; protein_id \"ENSP00000202017\";
20\tprotein_coding\texon\t30322356\t30322436\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"4\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; exon_id \"ENSE00000661137\";
20\tprotein_coding\tCDS\t30322356\t30322436\t.\t-\t2\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"4\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; protein_id \"ENSP00000202017\";
20\tprotein_coding\texon\t30320853\t30321749\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"5\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; exon_id \"ENSE00001109504\";
20\tprotein_coding\tCDS\t30321670\t30321749\t.\t-\t2\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; exon_number \"5\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\"; protein_id \"ENSP00000202017\";
20\tprotein_coding\tUTR\t30320853\t30321669\t.\t-\t.\tgene_id \"ENSG00000088356\"; transcript_id \"ENST00000202017\"; gene_name \"C20orf126\"; gene_source \"vega\"; gene_biotype \"protein_coding\"; transcript_name \"C20orf126\"; transcript_source \"vega\";
",


   ENST00000246203 => "#!genome-version NCBI33
20\tprotein_coding\ttranscript\t30565065\t30566129\t.\t-\t.\tgene_id \"ENSG00000125979\"; transcript_id \"ENST00000246203\"; gene_name \"TSPYL3\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"TSPYL3\"; transcript_source \"ensembl\";
20\tprotein_coding\texon\t30565065\t30566129\t.\t-\t.\tgene_id \"ENSG00000125979\"; transcript_id \"ENST00000246203\"; exon_number \"1\"; gene_name \"TSPYL3\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"TSPYL3\"; transcript_source \"ensembl\"; exon_id \"ENSE00000859919\";
20\tprotein_coding\tCDS\t30565065\t30566129\t.\t-\t0\tgene_id \"ENSG00000125979\"; transcript_id \"ENST00000246203\"; exon_number \"1\"; gene_name \"TSPYL3\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"TSPYL3\"; transcript_source \"ensembl\"; protein_id \"ENSP00000246203\";
20\tprotein_coding\tstart_codon\t30566127\t30566129\t.\t-\t0\tgene_id \"ENSG00000125979\"; transcript_id \"ENST00000246203\"; exon_number \"1\"; gene_name \"TSPYL3\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"TSPYL3\"; transcript_source \"ensembl\";
",

   ENST00000201961 => "#!genome-version NCBI33
20\tprotein_coding\ttranscript\t30885729\t30911383\t.\t-\t.\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\";
20\tprotein_coding\texon\t30911297\t30911383\t.\t-\t.\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"1\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; exon_id \"ENSE00000661216\";
20\tprotein_coding\tCDS\t30911297\t30911383\t.\t-\t0\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"1\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; protein_id \"ENSP00000201961\";
20\tprotein_coding\tstart_codon\t30911381\t30911383\t.\t-\t0\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"1\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\";
20\tprotein_coding\texon\t30903618\t30903773\t.\t-\t.\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"2\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; exon_id \"ENSE00000661215\";
20\tprotein_coding\tCDS\t30903618\t30903773\t.\t-\t0\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"2\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; protein_id \"ENSP00000201961\";
20\tprotein_coding\texon\t30896671\t30896782\t.\t-\t.\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"3\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; exon_id \"ENSE00000560920\";
20\tprotein_coding\tCDS\t30896671\t30896782\t.\t-\t0\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"3\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; protein_id \"ENSP00000201961\";
20\tprotein_coding\texon\t30887207\t30887316\t.\t-\t.\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"4\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; exon_id \"ENSE00000560923\";
20\tprotein_coding\tCDS\t30887207\t30887316\t.\t-\t2\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"4\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; protein_id \"ENSP00000201961\";
20\tprotein_coding\texon\t30885729\t30885800\t.\t-\t.\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"5\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; exon_id \"ENSE00000661212\";
20\tprotein_coding\tCDS\t30885729\t30885800\t.\t-\t0\tgene_id \"ENSG00000088303\"; transcript_id \"ENST00000201961\"; exon_number \"5\"; gene_name \"Q9NQF5\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\"; transcript_name \"Q9NQF5\"; transcript_source \"ensembl\"; protein_id \"ENSP00000201961\";
"
};

my $transcript_adaptor = $db->get_TranscriptAdaptor;

my @keys = keys %{$transcripts_gtf};
foreach my $transcript_id (sort @keys) {
  my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
  defined $transcript or 
    $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);

  skip "Cannot retrieve transcript $transcript_id. Skipping test", 1
    unless defined $transcript;

  my $fh = IO::String->new();
  my $gtf_serializer = 
    Bio::EnsEMBL::Utils::IO::GTFSerializer->new($fh);
  $gtf_serializer->print_main_header($db);
  $gtf_serializer->print_feature($transcript);

  eq_or_diff(${$fh->string_ref()}, $transcripts_gtf->{$transcript_id}, "Transcript $transcript_id serialises to GTF as expected");
    
}

done_testing();

