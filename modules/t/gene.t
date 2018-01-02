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

use Test::More;
use Test::Warnings qw(warning);
use Test::Exception;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DnaDnaAlignFeature;

# switch on the debug prints
our $verbose = 0;

debug("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor("core");

debug("Test database instatiated");
ok($db);

my $ontology = Bio::EnsEMBL::Test::MultiTestDB->new('ontology');

my $odb = $ontology->get_DBAdaptor("ontology");
note("Ontology database instatiated");
ok($odb);
my $go_adaptor;
warning { $go_adaptor = $odb->get_OntologyTermAdaptor(); };

my $gene;
my $ga = $db->get_GeneAdaptor();

debug("Gene->list_dbIDs");
my $ids = $ga->list_dbIDs();
ok(@{$ids});

debug("Gene->list_stable_ids");
my $stable_ids = $ga->list_stable_ids();
ok(@{$stable_ids});

debug("Gene->list_seq_region_ids");
my $region_ids = $ga->list_seq_region_ids();
ok(@{$region_ids});

$gene = $ga->fetch_by_display_label("T9S4_HUMAN");
ok($gene && $gene->dbID() == 18262);

$gene = $ga->fetch_by_stable_id("ENSG00000171456");
debug("Gene->fetch_by_stable_id()");
ok($gene);

my @date_time = localtime($gene->created_date());
ok($date_time[3] == 6 && $date_time[4] == 11 && $date_time[5] == 104);

@date_time = localtime($gene->modified_date());
ok($date_time[3] == 6 && $date_time[4] == 11 && $date_time[5] == 104);

debug("Gene dbID: " . $gene->dbID());
ok($gene->dbID() == 18267);

debug("Gene start: " . $gene->start);
ok($gene->start() == 30735607);

debug("Gene end: " . $gene->end);
ok($gene->end() == 30815178);

debug("Gene external name: " . $gene->external_name);
ok($gene->external_name eq "Q9H466");

debug("Gene external dbname: " . $gene->external_db);
ok($gene->external_db eq "Uniprot/SPTREMBL");

debug("Gene display xref id: " . $gene->display_xref->dbID);
ok($gene->display_xref->dbID() == 128324);

# test the getters and setters
ok(test_getter_setter($gene, "external_name", "banana"));
ok(test_getter_setter($gene, "external_db",   "dummy"));
ok(test_getter_setter($gene, "display_xref",  42));
ok(test_getter_setter($gene, "created_date",  time()));
ok(test_getter_setter($gene, "modified_date", time()));

my $links = $gene->get_all_DBLinks();
debug("Links: " . scalar(@$links));

ok(scalar @$links == 6);

# now create a new gene ...

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region("chromosome", "20", 30_249_935, 31_254_640);

debug("Slice from SliceAdaptor");
ok($slice);

my $analysis   = $db->get_AnalysisAdaptor->fetch_by_logic_name("ensembl");
my $f_analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name("Vertrna");
debug("Analysis from AnalysisAdaptor");
ok($analysis);

$gene = Bio::EnsEMBL::Gene->new();

my $transcript1 = Bio::EnsEMBL::Transcript->new();
$transcript1->analysis($analysis);
my $transcript2 = Bio::EnsEMBL::Transcript->new();
$transcript2->analysis($analysis);

my $ex1 = Bio::EnsEMBL::Exon->new();
my $ex2 = Bio::EnsEMBL::Exon->new();
my $ex3 = Bio::EnsEMBL::Exon->new();

my $translation1 = Bio::EnsEMBL::Translation->new();
my $translation2 = Bio::EnsEMBL::Translation->new();

ok($gene);

is($gene->version,         1, 'Default gene version = 1');
is($transcript1->version,  1, 'Default transcript version = 1');
is($ex1->version,          1, 'Default exon version = 1');
is($translation1->version, 1, 'Default translation version = 1');

$ex1->start(13586);
$ex1->end(13735);
$ex1->phase(0);
$ex1->end_phase(0);
$ex1->slice($slice);
$ex1->strand(1);
$ex1->analysis($analysis);

my @feats;
my $fp = new Bio::EnsEMBL::FeaturePair;

$fp->start(13586);
$fp->end(13705);
$fp->strand(1);
$fp->score(10);
$fp->slice($slice);
$fp->hstart(100);
$fp->hend(219);
$fp->hstrand(1);
$fp->hseqname('dummy-hid');

push(@feats, $fp);

$fp = new Bio::EnsEMBL::FeaturePair;
$fp->start(13707);
$fp->end(13735);
$fp->strand(1);
$fp->score(10);
$fp->slice($slice);

$fp->hstart(220);
$fp->hend(248);
$fp->hstrand(1);
$fp->hseqname('dummy-hid');
push(@feats, $fp);

#
#
# 2 Test DnaDnaAlignFeature::new(-features)
#
my $dnaf;
warning { $dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@feats); };
$dnaf->analysis($f_analysis);

$ex1->add_supporting_features($dnaf);

$ex2->start(201372);
$ex2->end(201571);
$ex2->phase(0);
$ex2->end_phase(-1);
$ex2->slice($slice);
$ex2->strand(1);
$ex2->analysis($analysis);

@feats = ();
$fp    = new Bio::EnsEMBL::FeaturePair;

$fp->start(201372);
$fp->end(201471);
$fp->strand(1);
$fp->score(10);
$fp->slice($slice);
$fp->hstart(100);
$fp->hend(199);
$fp->hstrand(1);
$fp->hseqname('dummy-hid');

push(@feats, $fp);

$fp = new Bio::EnsEMBL::FeaturePair;
$fp->start(201472);
$fp->end(201571);
$fp->strand(1);
$fp->score(10);
$fp->slice($slice);

$fp->hstart(201);
$fp->hend(300);
$fp->hstrand(1);
$fp->hseqname('dummy-hid');
push(@feats, $fp);

#
#
# 2 Test DnaDnaAlignFeature::new(-features)
#
warning { $dnaf = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@feats); };
$dnaf->analysis($f_analysis);

$ex2->add_supporting_features($dnaf);

$ex3->start(210600);
$ex3->end(210800);
$ex3->phase(-1);
$ex3->end_phase(-1);
$ex3->slice($slice);
$ex3->strand(1);
$ex3->analysis($analysis);

$transcript1->add_Exon($ex1);
$transcript1->add_Exon($ex2);
$translation1->start_Exon($ex1);
$translation1->end_Exon($ex2);
$translation1->start(1);
$translation1->end(150);
$transcript1->translation($translation1);

$transcript2->add_Exon($ex1);
$transcript2->add_Exon($ex2);
$transcript2->add_Exon($ex3);
$translation2->start_Exon($ex1);
$translation2->end_Exon($ex2);
$translation2->start(1);
$translation2->end(180);
$transcript2->translation($translation2);

debug("Transcripts created");
ok($transcript1);

$gene->add_Transcript($transcript1);
$gene->add_Transcript($transcript2);

$gene->analysis($analysis);

debug("Getting all the Transcripts/Exons from new Gene");

my $count      = 0;
my $translates = 1;

foreach my $tr (@{$gene->get_all_Transcripts()}) {
  if ($tr->translate()->seq() =~ /\*./) {
	$translates = 0;
	debug("Translate failed.");
  }
  debug("Translation: " . $tr->translate()->seq());
  foreach my $exon (@{$tr->get_all_Exons()}) {
	debug("  Exon start: " . $exon->start());
	debug("  Exon end:   " . $exon->end());
	debug("  Exon strand " . $exon->strand());
	$count++;
  }
}

ok($count == 5);
ok($translates);

# Verify Transcript cache is not leaky
my $transcripts = $gene->get_all_Transcripts;
$count = @$transcripts;
pop @$transcripts;
$transcripts = $gene->get_all_Transcripts;
cmp_ok(scalar @$transcripts, '==', $count, "Gene's transcript cache is not modified by changing transcript lists in caller code");

ok(scalar(@{$gene->get_all_Exons()}) == 3);

$gene = $gene->transform("chromosome");

my $desc      = 'test description for a gene';
my $stable_id = 'ENSG00000171456';
$gene->description($desc);
$gene->stable_id($stable_id);

$gene->created_date( time());
$gene->modified_date(time());

$multi->hide("core", "meta_coord", "gene", "transcript", "exon", "exon_transcript", "translation", "supporting_feature", "dna_align_feature", 'xref', 'object_xref', 'identity_xref');

my $gene_ad = $db->get_GeneAdaptor();
debug("Storing the gene");
$gene_ad->store($gene);

ok(1);

# Cache all mappings needed for genes
$gene_ad->cache_gene_seq_mappings();

my $genes = $slice->get_all_Genes();

ok(scalar(@$genes) == 1);

my $gene_out = $genes->[0];

#make sure the stable_id was stored
ok($gene_out->stable_id eq $stable_id);

#make sure the version was stored
ok($gene_out->version == 1);

#make sure the description was stored
ok($gene_out->description eq $desc);

debug("gene_out created_date   = ", $gene_out->created_date);
debug("gene_out modified_date  = ", $gene_out->modified_date);

is($gene_out->created_date,  $gene->created_date,  'created_date roundtrips');
is($gene_out->modified_date, $gene->modified_date, 'modified_date roundtrips');

ok(scalar(@{$gene_out->get_all_Exons()}) == 3);

foreach my $tr (@{$gene_out->get_all_Transcripts()}) {
  debug("NewTranscript: " . $tr->dbID());
  foreach my $exon (@{$tr->get_all_Exons()}) {
	debug("  NewExon: " . $exon->start() . " " . $exon->end() . " " . $exon->strand());
  }
}

my $exons = $gene_out->get_all_Transcripts()->[0]->get_all_Exons();

ok($exons->[0]->start == 13586);

ok($exons->[1]->strand == 1);
ok($exons->[1]->phase == 0);

my $pep;
my $translate = 0;
foreach my $trans (@{$gene_out->get_all_Transcripts()}) {

  my $pep = $trans->translate();
  debug("Peptide: " . $pep->seq());

  if ($pep->seq !~ /\*./) {
	$translate = 1;
  } else {
	$translate = 0;
  }
}

ok($translate == 1);

my $t = $gene_out->get_all_Transcripts()->[1];

my $e  = $t->get_all_Exons()->[0];
my $se = $e->get_all_supporting_features();

debug("Got " . scalar(@$se) . " supporting features.");
ok(scalar(@$se) == 1);

my $se_start = $se->[0]->start();

my $se_end = $se->[0]->end();

debug("Supporting start $se_start, end $se_end");
debug("Exon start " . $e->start() . " end " . $e->end());

ok($se_start == $e->start());
ok($se_end == $e->end());

my $pep1 = $t->translate()->seq();

$e->phase(1);
my $pep2 = $t->translate()->seq();

debug("Pep phase 0: $pep1");
debug("Pep phase 1: $pep2");

ok($pep1 ne $pep2);
debug("checking external references");

$multi->restore();

$slice = $db->get_SliceAdaptor()->fetch_by_region("chromosome", "20", 30_252_000, 31_252_001);

$genes = $slice->get_all_Genes();

# try and count the genes on the slice
note 'Processing counts';
my $geneCount = $ga->count_all_by_Slice($slice);
is(scalar(@$genes), $geneCount, 'Gene count of array on Chr 20 subset');
$geneCount = $ga->count_all_by_Slice($slice, 'protein_coding');
is(scalar(@$genes), $geneCount, 'Protein Coding gene count of array on Chr 20 subset');
$geneCount = $ga->count_all_by_Slice($slice, 'banana');
is($geneCount, 0, 'Gene count on Chr 20 subset with bogus biotype');
$geneCount = $ga->count_all_by_Slice($slice, ['banana', 'protein_coding']);
is($geneCount, scalar(@$genes), 'Protein coding gene count matches array size on Chr20 subset');
$geneCount = $ga->count_all_by_Slice($slice, 'protein_coding', 'ensembl');
my $vega_geneCount = $ga->count_all_by_Slice($slice, 'protein_coding', 'vega');
is($geneCount, scalar(@$genes) - $vega_geneCount, "Almost all genes are of source ensembl");

# Time to do more complex counts involving slice projections
{
  my $hap_slice = $db->get_SliceAdaptor()->fetch_by_region('chromosome', '20_HAP1', 4_000_000);
  my $hap_gene_array = $hap_slice->get_all_Genes();
  is(scalar(@{$hap_gene_array}), 14, 'Expected amount of genes when fetching all from a HAP region');
  my $gene_count = $ga->count_all_by_Slice($hap_slice);
  is($gene_count, 14, 'Counts should span slices when given patches/haps/pars');
}


#save contents of gene table
$multi->save('core', 'gene');

# tests for update method
# go get a fresh gene again
$gene = $ga->fetch_by_stable_id("ENSG00000171456");

# the first update should no effect
$ga->update($gene);

my $newgene = $ga->fetch_by_stable_id("ENSG00000171456");
ok($newgene->display_xref->dbID() == 128324);
ok($newgene->biotype eq 'protein_coding');

# now change the original gene and update it
my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

$gene->display_xref($dbEntryAdaptor->fetch_by_dbID(614));
$gene->biotype('dummy');
$ga->update($gene);

$newgene = $ga->fetch_by_stable_id("ENSG00000171456");
ok($newgene->display_xref->dbID() == 614);
ok($newgene->biotype eq 'dummy');

$multi->restore('core', 'gene');

# test update_coords method
# old coords: 30735607 - 30815178
$gene = $ga->fetch_by_stable_id("ENSG00000171456");
$gene->start(30730000);
$gene->end(30815178);
$ga->update_coords($gene);

my $new_gene = $ga->fetch_by_stable_id("ENSG00000171456");
cmp_ok($new_gene->start(), '==', 30735607, 'Updated gene start');
cmp_ok($new_gene->end(), '==', 30815178, 'Updated gene end');

#
# test GeneAdaptor::fetch_all_by_domain
#
my @genes = @{$ga->fetch_all_by_domain('IPR000010')};

debug("Fetch by domain 'IPR000010'");

ok(@genes == 2);
debug("Got " . scalar(@genes) . " genes");
ok(($genes[0]->stable_id() eq 'ENSG00000131044') || ($genes[1]->stable_id() eq 'ENSG00000131044'));
ok(($genes[0]->stable_id() eq 'ENSG00000174873') || ($genes[1]->stable_id() eq 'ENSG00000174873'));

#
# test GeneAdaptor::fetch_all_by_external_name
#

#Q15691
($gene) = @{$ga->fetch_all_by_external_name('MAE1_HUMAN')};
debug($gene->stable_id);
ok($gene->stable_id() eq 'ENSG00000101367');

#
# test GeneAdaptor::fetch_all_by_Slice
#
@genes = @{ $ga->fetch_all_by_Slice($slice) };
is(scalar(@genes), 19, "Found 19 genes with fetch_all_by_Slice");
is($genes[0]->stable_id(), 'ENSG00000131044', "First gene is ENSG00000131044");
is($genes[1]->stable_id(), 'ENSG00000174873', "Second gene is ENSG00000174873");

#
# test GeneAdaptor::fetch_all_by_Slice_and_external_dbname_link
#
@genes = @{ $ga->fetch_all_by_Slice_and_external_dbname_link($slice, undef, 0, "HUGO") };
is(scalar(@genes), 13, "Found 13 genes with fetch_all_by_Slice_and_external_dbname_link");
is($genes[0]->stable_id(), 'ENSG00000131044', "First gene is ENSG00000131044");
is($genes[1]->stable_id(), 'ENSG00000088356', "Second gene is ENSG00000088356");

warning { @genes = @{ $ga->fetch_all_by_Slice_and_external_dbname_link($slice, undef, 0, "random") }; };
is(scalar(@genes), 0, "No genes with db random");

#
# test GeneAdaptor::fetch_all_by_display_label
#
@genes = @{ $ga->fetch_all_by_display_label("C20orf125") };
is(scalar(@genes), 1, "Found 1 genes with fetch_all_by_display_label");
is($genes[0]->stable_id(), 'ENSG00000131044', "First gene is ENSG00000131044");

#
# test GeneAdaptor::fetch_all_by_transcript_supporting_evidence
#
@genes = @{ $ga->fetch_all_by_transcript_supporting_evidence('Q9NUG5', 'protein_align_feature') };
is(scalar(@genes), 0, "Found 0 genes with fetch_all_by_transcript_supporting_evidence");


#
# test GeneAdaptor::fetch_all_by_exon_supporting_evidence
#
@genes = @{ $ga->fetch_all_by_exon_supporting_evidence('BF346221.1', 'dna_align_feature') };
is(scalar(@genes), 1, "Found 1 genes with fetch_all_by_exon_supporting_evidence");
my $aa = $db->get_AnalysisAdaptor();
$analysis = $aa->fetch_by_logic_name('RepeatMask');
@genes = @{ $ga->fetch_all_by_exon_supporting_evidence('BF346221.1', 'dna_align_feature', $analysis) };
is(scalar(@genes), 0, "No genes for random analysis");

#
# test GeneAdaptor::fetch_all_by_logic_name
#
@genes = @{ $ga->fetch_all_by_logic_name('ensembl') };
is(scalar(@genes), 21, "All genes of analysis ensembl");


#
# test GeneAdaptor::fetch_all_by_GOTerm
#
my $go_term = $go_adaptor->fetch_by_accession('GO:0003677');
@genes = @{ $ga->fetch_all_by_GOTerm($go_term) };
is(scalar(@genes), 2, "Found 2 genes with fetch_all_by_GOTerm");


#
# test GeneAdaptor::fetch_all_by_GOTerm_accession
#
#@genes = @{ $ga->fetch_all_by_GOTerm_accession("GO:0004835") };
#is(scalar(@genes), 13, "Found 13 genes with fetch_all_by_GOTerm_accession");
#is($genes[0]->stable_id(), 'ENSG00000131044', "First gene is ENSG00000131044");

#
# test fetch_all_by_description
#

my $gene_list = $ga->fetch_all_by_description('%APC-BINDING PROTEIN EB1%');
ok(scalar(@$gene_list) == 1, "Get by description");
ok($gene_list->[0]->stable_id eq "ENSG00000101367", "check we got the right one by description");

#
# test fetch_all_by_external_name with wildcard restrictions
#
(@genes) = @{$ga->fetch_all_by_external_name('AF_%')};
# Should = 0 because _ is auto-escaped.
debug('Genes found under external_name AF_%: ' . scalar(@genes));
ok(scalar(@genes) == 0);
(@genes) = @{$ga->fetch_all_by_external_name('AF_%', undef, 'override')};
debug('Genes found under external_name AF_% with override on: ' . scalar(@genes));
debug($genes[0]->stable_id());
debug($genes[1]->stable_id());
debug($genes[2]->stable_id());
debug($genes[3]->stable_id());
# Note that 9 AF_% xrefs correspond to 4 unique ensembl IDs.
ok(scalar(@genes) == 4);
#
# test fetch_all_by_external_name with wildcard matching
#
@genes = @{$ga->fetch_all_by_external_name('MAE__HUMAN')};
debug("Wildcard test:" . $genes[0]->stable_id);
ok($genes[0]->stable_id() eq 'ENSG00000101367');

@genes = @{$ga->fetch_all_by_external_name('M_%')};
debug("Wildcard test:" . $genes[0]->stable_id());
debug(scalar @genes . " genes found");
ok(scalar @genes == 2);

# Test performance protection (very vague queries return no hits)
debug("Testing vague query protection");
{
  my $warnings = q{};
  local $SIG{'__WARN__'} = sub {
	$warnings .= $_[0];
  };
  ok(scalar(@{$ga->fetch_all_by_external_name('M%')}) == 0);
  ok(scalar(@{$ga->fetch_all_by_external_name('%')}) == 0);
  like($warnings, qr/is too vague and will monopolise database/, 'Checking for warnings being emitted by the above methods');
}

#
# test GeneAdaptor::get_Interpro_by_geneid
#
debug("Test get_Interpro_by_geneid");
my @interpro = @{$ga->get_Interpro_by_geneid('ENSG00000174873')};
ok(@interpro == 1);
debug($interpro[0]);

ok($gene->display_id eq $gene->stable_id);

#
# test GeneAdaptor::fetch_all_by_biotype
#
debug("Test fetch_all_by_biotype");
@genes = @{$ga->fetch_all_by_biotype('protein_coding')};
ok(@genes == 21, "Gene count for protein coding");
@genes = @{$ga->fetch_all_by_biotype(['protein_coding', 'sRNA'])};
ok(@genes == 21, "Gene count for protein coding and sRNA");
$geneCount = $ga->count_all_by_biotype('protein_coding');
ok($geneCount == 21, "Gene count via method call for protein coders");
$geneCount = $ga->count_all_by_biotype(['protein_coding', 'sRNA']);
ok($geneCount == 21, "Gene count via method call for protein coding and sRNA");

#
# test GeneAdaptor::fetch_all_by_source
#
note("Test fetch_all_by_source");
@genes = @{$ga->fetch_all_by_source('ensembl')};
note "Got ".scalar(@genes)." ensembl genes\n";
ok(@genes == 19);
$geneCount = $ga->count_all_by_source('ensembl');
ok($geneCount == 19);
@genes = @{$ga->fetch_all_by_source(['havana','vega'])};
note "Got ".scalar(@genes)." (havana, vega) transcripts\n";
ok(@genes == 2);
$geneCount = $ga->count_all_by_source(['havana', 'vega']);
ok($geneCount == 2);


#
# Test Gene: get_all_Introns
#

#
# Test get_all_Introns by joining Exons and introns
# and comparing it to the original
#

foreach my $stable_id (qw(ENSG00000174873 ENSG00000101367)){ #test both strands

  my $gene = $ga->fetch_by_stable_id($stable_id);

  my @exons = sort { $a->start <=> $b->start } (@{$gene->get_all_Exons()});
  my @introns = sort { $a->start <=> $b->start } (@{$gene->get_all_Introns()});

  my $orig_seq = $gene->slice->subseq(
                 $gene->start(),
                 $gene->end(),
                 $gene->strand());

  my $idl=0;
  my $new_seq = $exons[0]->seq()->seq();
  foreach my $intron (@introns){
    $new_seq .= $intron->seq;
    $new_seq .= $exons[$idl+1]->seq->seq();
    $idl++;

  }

  is ($gene->slice->subseq($gene->start(), $gene->start() + 9, $gene->strand()), substr($exons[0]->seq()->seq(), 0, 10), "Initial sequences match");

  is (length($orig_seq), length($new_seq), "Sequence length match");
  is($orig_seq, $new_seq, 'Correct new origin seq');

}

#
# test Gene: get_all_alt_alleles
#

is($gene->is_reference, 1, "If no alt allele, gene is reference");
$gene = $ga->fetch_by_dbID(18256);
my $alt_genes = $gene->get_all_alt_alleles();

ok($gene->is_reference == 1);
# expect the following alleles
is_deeply(
  [sort {$a <=> $b} map {$_->dbID()} @{$alt_genes}], 
  [18257,18258,18259], 
  'Checking retrieved alt allele IDs are expected'
) or diag explain $alt_genes;

#
# test storing a new allele group
#

$multi->hide('core', 'alt_allele');

#TODO Fix. Code current raises this warning. Don't use alt genes on a ref slice or change the code
#
#-------------------- WARNING ----------------------
#MSG: More than one alternative allele on the reference sequence (gene ids: 18270,18271,18272). Ignoring.
#FILE: EnsEMBL/DBSQL/GeneAdaptor.pm LINE: 1073
#CALLED BY: modules/t/gene.t  LINE: 514
#Ensembl API version = 67
#---------------------------------------------------

my @alt_genes;
push(@alt_genes, $ga->fetch_by_dbID(18270));
push(@alt_genes, $ga->fetch_by_dbID(18271));
push(@alt_genes, $ga->fetch_by_dbID(18272));

note 'Do not use GeneAdaptor::store_alt_alleles() for storing alt alleles. Ensuring we still warn';
warns_like {
  $ga->store_alt_alleles(\@alt_genes);
} qr/Unsupported.+alternative.+reference\ssequence.+Ignoring/is, 'Checking we are still warning about multiple alt_alleles on refs';
$gene = $ga->fetch_by_dbID(18270);
is(scalar(@{$gene->get_all_alt_alleles()}), 0, 'Checking we have no alleles retrieved because we had too many representatives (too many ref genes)');

#
# Gene remove test
#

warning { $multi->save("core", "gene", "transcript", "translation", "protein_feature", "exon", "exon_transcript", "supporting_feature", "object_xref", "ontology_xref", "identity_xref", "dna_align_feature", "protein_align_feature", 'meta_coord'); };

$gene = $ga->fetch_by_stable_id("ENSG00000171456");

my $gene_count  = count_rows($db, "gene");
my $exon_count  = count_rows($db, "exon");
my $trans_count = count_rows($db, "transcript");
my $tl_count    = count_rows($db, "translation");

my $tminus = scalar(@{$gene->get_all_Transcripts()});
my $eminus = scalar(@{$gene->get_all_Exons()});

debug("Genes before " . $gene_count);
debug("Exons before " . $exon_count);
debug("Transcripts before " . $trans_count);
debug("Translations before " . $tl_count);
debug("Gene has " . $tminus . " transcripts");
debug("Gene has " . $eminus . " exons");

$ga->remove($gene);

ok(count_rows($db, "gene") ==       ($gene_count - 1));
ok(count_rows($db, "transcript") == ($trans_count - $tminus));
ok(count_rows($db, "exon") ==       ($exon_count - $eminus));

ok(!defined($gene->dbID()));
ok(!defined($gene->adaptor()));

$multi->restore('core');

#
# regression test - test the recalculation of coords
# in the Gene.  This was setting the end incorrectly
# before.
#
$gene = Bio::EnsEMBL::Gene->new();

$gene->slice($slice);

my $first_ex = Bio::EnsEMBL::Exon->new(-START     => 10,
									   -END       => 100,
									   -STRAND    => 1,
									   -PHASE     => 0,
									   -END_PHASE => 1,
									   -SLICE     => $slice);

my $second_ex = Bio::EnsEMBL::Exon->new(-START     => 200,
										-END       => 400,
										-STRAND    => 1,
										-PHASE     => 1,
										-END_PHASE => 0,
										-SLICE     => $slice);

$transcript1 = Bio::EnsEMBL::Transcript->new(-EXONS => [$first_ex, $second_ex]);
$transcript1->analysis($analysis);

$transcript2 = Bio::EnsEMBL::Transcript->new(-EXONS => [$first_ex]);
$transcript2->analysis($analysis);

$gene->add_Transcript($transcript1);
$gene->add_Transcript($transcript2);

$gene->recalculate_coordinates();

ok($gene->start() == 10);
ok($gene->end() == 400);

{
  # Test transformming genes over a gapped alignment, when gaps
  # only occur in introns
  # target_slice = sliceAdaptor->fetch_by_region( "alt_chrom", "gap_map_test" );

  my $gene = $ga->fetch_by_dbID(18274);
  debug(":::: Gene: $gene");
  my $new_gene = $gene->transform("alt_chrom");
  debug(":::: New Gene $new_gene");
  my $trans_orig   = $gene->get_all_Transcripts()->[0];
  my $trans_mapped = $new_gene->get_all_Transcripts()->[0];
  ok($trans_orig->spliced_seq() eq $trans_mapped->spliced_seq());

  # the assembly is rigged to have no gaps between the first and second exon
  ok($trans_mapped->get_all_Exons()->[0]->end() + 1 == $trans_mapped->get_all_Exons()->[1]->start());
}

#
# test that the display_xref_id is set for the gene and its transcript
#
$multi->hide("core", "gene", "transcript", "exon", 'xref', 'object_xref', "exon_transcript", "translation", 'meta_coord');

$gene->analysis($analysis);
my $analysis_adap = $db->get_AnalysisAdaptor();
$analysis = $analysis_adap->fetch_by_logic_name("ensembl");

my $dbe = Bio::EnsEMBL::DBEntry->new(-primary_id => 'test_id',
									 -version    => 1,
									 -dbname     => 'EMBL',
									 -release    => 1,
									 -display_id => 'test_id',
									 -analysis   => $analysis);

$gene->add_DBEntry($dbe);
$gene->display_xref($dbe);

$gene->get_all_Transcripts()->[0]->add_DBEntry($dbe);
$gene->get_all_Transcripts()->[0]->display_xref($dbe);

debug("Storing gene");
$gene_ad->store($gene);

my $dbe_id = $db->dbc->db_handle->selectall_arrayref("SELECT display_xref_id FROM gene")->[0]->[0];

ok($dbe_id && $dbe_id == $dbe->dbID());

$multi->restore();

#
# tests for multiple versions of genes in a database
#

$gene = $ga->fetch_by_stable_id('ENSG00000355555');
debug("fetch_by_stable_id");
ok($gene->dbID == 18275);

$gene->stable_id_version('ENSG00000171455.4');
is($gene->stable_id, 'ENSG00000171455', 'Stable id set with stable_id_version');
is($gene->version, 4, 'Version set with stable_id_version');
is($gene->stable_id_version, 'ENSG00000171455.4', 'Stable id and version from stable_id_version');

$gene->stable_id_version('ENSG00000171456');
is($gene->stable_id, 'ENSG00000171456', 'Stable id set with stable_id_version');
is($gene->version, undef, 'Version undef from stable_id_version');
is($gene->stable_id_version, 'ENSG00000171456', 'Stable id and no version from stable_id_version');

$gene = $ga->fetch_by_stable_id("ENSG00000171456.1");
ok($gene->stable_id eq 'ENSG00000171456', "Fetch by stable_id with version");

$gene = $ga->fetch_by_stable_id_version("ENSG00000171456", 1);
ok($gene->stable_id eq 'ENSG00000171456', "fetch_by_stable_id_version");

$gene = $ga->fetch_by_stable_id("ENSG00000171456.1a");
ok(! defined($gene), "Fetch by stable_id with bad version");

$gene = $ga->fetch_by_stable_id_version("ENSG00000171456", '1a');
ok(! defined($gene), "fetch_by_stable_id_version, with bad version");

@genes = @{$ga->fetch_all_versions_by_stable_id('ENSG00000355555')};
debug("fetch_all_versions_by_stable_id");
ok(scalar(@genes) == 1);

my $sl = $sa->fetch_by_region('chromosome', 'MT_NC_001807');
@genes = @{$sl->get_all_Genes};
ok(scalar(@genes) == 1);

$gene = $ga->fetch_by_transcript_stable_id('ENST00000355555');
debug("fetch_by_transcript_stable_id");
ok($gene->dbID == 18275);

$gene = $ga->fetch_by_translation_stable_id('ENSP00000355555');
debug("fetch_by_translation_stable_id");
ok($gene->dbID == 18275);

$gene = $ga->fetch_by_translation_stable_id('random_ENSP00000355555');
debug("fetch_by_translation_stable_id");
is($gene, undef, "No gene for random translation id");

$gene = $ga->fetch_by_exon_stable_id('ENSE00001109603');
debug("fetch_by_exon_stable_id");
ok($gene->dbID == 18275);

@genes = @{$ga->fetch_all_by_external_name('PAL2_HUMAN')};
debug("fetch_all_by_external_name");
ok(scalar(@genes) == 1 && $genes[0]->dbID == 18264);

$gene = $ga->fetch_by_display_label('PLAGL2');
debug("fetch_by_display_label");
ok($gene->dbID == 18264);

$gene = $ga->fetch_by_dbID(18264);
debug("fetch_by_dbID, current");
ok($gene->is_current == 1);

#$gene = $ga->fetch_by_dbID(18264);
#debug("fetch_by_dbID, non current");
#ok( $gene->is_current == 1 );

# store/update

$gene = $ga->fetch_by_stable_id('ENSG00000355555');
foreach my $t (@{$gene->get_all_Transcripts}) {
  $t->get_all_Exons;
}

$multi->hide("core", "gene", "transcript", "exon", 'xref', 'object_xref', "exon_transcript", "translation", 'meta_coord');

$gene->version(3);
$gene->dbID(undef);
$gene->adaptor(undef);
$ga->store($gene);

$gene->version(4);
$gene->is_current(0);
$gene->dbID(undef);
$gene->adaptor(undef);
$ga->store($gene);

$gene->version(undef);
$gene->is_current(0);
$gene->dbID(undef);
$gene->adaptor(undef);
$ga->store($gene);

$gene = $ga->fetch_by_stable_id('ENSG00000355555');
ok($gene->is_current == 1);

@genes = @{$ga->fetch_all_versions_by_stable_id('ENSG00000355555')};
foreach my $g (@genes) {
  if (defined $g->version && $g->version == 4) {
    ok($g->is_current == 0);
  }
}

$gene->is_current(0);
$ga->update($gene);
my $g1 = $ga->fetch_by_stable_id('ENSG00000355555');
ok(!$g1);

$gene->is_current(1);
$ga->update($gene);
$gene = $ga->fetch_by_stable_id('ENSG00000355555');
ok($gene->is_current == 1);

my $null_versions = 0;
foreach my $g (@genes) {
  if (! defined $g->version) {
    $null_versions++;
  }
}
is ( $null_versions, 1, "Null/undef version stored and retrieved");

$ga->remove_by_Slice($slice);
$geneCount = $ga->count_all_by_Slice($slice);
is($geneCount, 0, "Genes have been removed");

$multi->restore;

SKIP: {
  skip 'No registry support for SQLite yet', 1 if $db->dbc->driver() eq 'SQLite';

  #test the get_species_and_object_type method from the Registry
  my $registry = 'Bio::EnsEMBL::Registry';
  my ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('ENSG00000355555');
  is( $species, 'homo_sapiens', 'Test the get_species_and_object_type method from the Registry, species');
  is( $object_type, 'Gene', 'Test the get_species_and_object_type method from the Registry, object_type');

  ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('ENSG00000355555.1');
  is( $species, 'homo_sapiens', 'Test the get_species_and_object_type method from the Registry with version, species');
  is( $object_type, 'Gene', 'Test the get_species_and_object_type method from the Registry with version, object_type');

  ( $species, $object_type, $db_type ) = $registry->get_species_and_object_type('ENSG00000355555.2a');
  ok( !defined($species), 'Test the get_species_and_object_type method from the Registry with wrong version, species');

}

# Testing compara dba retrieval
{
  dies_ok { $gene->get_all_homologous_Genes(); } 'No Compara DBAdaptor has been configured. No way to retrieve data';
}

# Checking for External DB fetching by Ontology linkage type
{
  my $genes = $ga->fetch_all_by_ontology_linkage_type('GO', 'IDA');
  is(scalar(@{$genes}), 9, 'Expect 9 genes linked to IDAs'); 
  foreach my $gene (@{$genes}) {
    my %linkage_types = 
      map { $_, 1 }
      map { @{$_->get_all_linkage_types()} }
      @{$gene->get_all_DBLinks('GO')};
    ok($linkage_types{'IDA'}, $gene->stable_id().' was linked to an IDA term. Searching through all the links until we find the IDA linkage');
    
    my %undef_linkage_types = 
      map { $_, 1 }
      map { @{$_->get_all_linkage_types()} }
      grep { $_->can('get_all_linkage_types') } 
      @{$gene->get_all_DBLinks(undef)};
    ok($undef_linkage_types{'IDA'}, $gene->stable_id().' was linked to an IDA term found by using an undef external db. Searching through all the links until we find the IDA linkage');
  }
  
}

# Fetching by slice and an external DB
{
  $ga->clear_cache(); # have to clear the cache because otherwise it gets the wrong values back!
  my $local_slice = $sa->fetch_by_region("chromosome", "20", 30_249_935, 31_254_640);
  my $genes = $ga->fetch_all_by_Slice_and_external_dbname_link($local_slice, undef, undef, 'HUGO'); # yes HUGO. Not HGNC. Old data
  is(scalar(@{$genes}), 13, 'Expect 13 genes with HUGO/HGNC links');
  # assume this will be the display xref
  is($_->display_xref()->dbname(), 'HUGO', $_->stable_id().' has a display HUGO') for @{$genes};
}

done_testing();
