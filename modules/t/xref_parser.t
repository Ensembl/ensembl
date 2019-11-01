=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
use strict;
use warnings;

use File::Spec;
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;
use Test::More;
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use XrefParser::Database;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $dba = $multi->get_DBAdaptor("xref");
plan skip_all => "xref database schema is mysql specific - won't work with a different driver"
   unless $dba->dbc->driver eq 'mysql';

my $database = XrefParser::Database->new( $dba->dbc);

$database->populate(
  File::Spec->catdir($multi->curr_dir, "../../misc-scripts/xref_mapping"),
  "with force please",
);


# UniProt species taxon codes happen to match species ids of the core database
my $SPECIES_ID = 6239;
my $SPECIES_NAME = "Caenorhabditis elegans";

my %xref_tables_expected_empty_by_default = (
  checksum_xref=>0,
  coordinate_xref=>0,
  dependent_xref=>0,
  gene_direct_xref=>0,
  go_xref=>0,
  identity_xref=>0,
  object_xref=>0,
  primary_xref=>0,
  transcript_direct_xref=>0,
  translation_direct_xref=>0,
  xref=>0,
);

my $tmp_dir = tempdir(CLEANUP=>1);
sub store_in_temporary_file {
  my ($content, %opts) = @_;
  my $path = join("/", $tmp_dir, $opts{tmp_file_name} || "tmp");
  open(my $fh, ">", $path) or die $path;  
  print $fh $content;
  close($fh);
  return $path;
}

sub test_parser {
  my ($parser, $content, $expected, $test_name, %opts) = @_;
  require_ok($parser);
  $parser->new($database)->run({
   files =>[store_in_temporary_file($content, %opts)],
   source_id => "Source id (unused but sometimes required)",
   species_id => $SPECIES_ID,
   species => $SPECIES_NAME,
  });
  my $expected_table_counts = {%xref_tables_expected_empty_by_default, %$expected};
  subtest "$parser $test_name" => sub {
    plan tests => scalar(keys %$expected_table_counts);
    for my $table (keys %$expected_table_counts){
      my $actual_count = count_rows($dba, $table);
      $dba->dbc->prepare("delete from $table;")->execute() if ($actual_count and not $opts{skip_clean});
      my $expected_count = $expected_table_counts->{$table};
      is($actual_count, $expected_count, "$table has $expected_count rows") or diag "$table has $actual_count rows";
    }
  }
}
test_parser("XrefParser::WormbaseDirectParser", "",  {}, "null case");
my $wormbase_celegans_xrefs_head= <<EOF;
//
// WormBase Caenorhabditis elegans XREFs for WS265
//
// Columns (tab separated) are:
//    1. WormBase Gene sequence name
//    2. WormBase Gene accession
//    3. WormBase Gene CGC name
//    4. WormBase Transcript sequence name
//    5. WormPep protein accession
//    6. INSDC parent sequence accession
//    7. INSDC locus_tag id
//    8. INSDC protein_id
//    9. UniProt accession
//
// Missing or not applicable data (e.g. protein identifiers for non-coding RNAs) is denoted by a "."
//
2L52.1	WBGene00007063	.	2L52.1b	CE50569	BX284602	CELE_2L52.1	CTQ86426	A0A0K3AWR5
2L52.1	WBGene00007063	.	2L52.1a.1	CE32090	BX284602	CELE_2L52.1	CCD61130	A4F336
2L52.1	WBGene00007063	.	2L52.1a.2	CE32090	BX284602	CELE_2L52.1	CCD61130	A4F336
2L52.2	WBGene00200402	.	2L52.2	.	BX284602	CELE_2L52.2	.	.
EOF
my $wormbase_celegans_xrefs_expected_count = {
xref=>14,
gene_direct_xref => 4,
transcript_direct_xref => 7,
translation_direct_xref => 6,
};
test_parser("XrefParser::WormbaseDirectParser", $wormbase_celegans_xrefs_head,  
   $wormbase_celegans_xrefs_expected_count, "Direct xrefs: genes: columns 1,2,3, transcripts: column 4 once as transcript and also once as CDS for coding genes, translations: column 5 and 7. xrefs: sum of these minus one reused CDS and two reused proteins"
);

my $uniprot_elegans_record = <<EOF;
ID   A0A0K3AWR5_CAEEL        Unreviewed;       220 AA.
AC   A0A0K3AWR5;
DT   11-NOV-2015, integrated into UniProtKB/TrEMBL.
DT   11-NOV-2015, sequence version 1.
DT   18-JUL-2018, entry version 14.
DE   SubName: Full=Uncharacterized protein {ECO:0000313|EMBL:CTQ86426.1};
GN   ORFNames=2L52.1 {ECO:0000313|EMBL:CTQ86426.1,
GN   ECO:0000313|WormBase:2L52.1b},
GN   CELE_2L52.1 {ECO:0000313|EMBL:CTQ86426.1};
OS   $SPECIES_NAME.
OC   Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida;
OC   Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis.
OX   NCBI_TaxID=$SPECIES_ID {ECO:0000313|EMBL:CTQ86426.1, ECO:0000313|Proteomes:UP000001940};
RN   [1] {ECO:0000313|EMBL:CTQ86426.1, ECO:0000313|Proteomes:UP000001940}
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RC   STRAIN=Bristol N2 {ECO:0000313|EMBL:CTQ86426.1,
RC   ECO:0000313|Proteomes:UP000001940};
RX   PubMed=9851916; DOI=https://doi.org/10.1126/science.282.5396.2012;
RG   The C. elegans sequencing consortium;
RA   Sulson J.E., Waterston R.;
RT   "Genome sequence of the nematode C. elegans: a platform for
RT   investigating biology.";
RL   Science 282:2012-2018(1998).
CC   -----------------------------------------------------------------------
CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
CC   Distributed under the Creative Commons Attribution (CC BY 4.0) License
CC   -----------------------------------------------------------------------
DR   EMBL; BX284602; CTQ86426.1; -; Genomic_DNA.
DR   RefSeq; NP_001300487.1; NM_001313558.1.
DR   UniGene; Cel.25279; -.
DR   EnsemblMetazoa; 2L52.1b; 2L52.1b; WBGene00007063.
DR   GeneID; 181792; -.
DR   CTD; 181792; -.
DR   WormBase; 2L52.1b; CE50569; WBGene00007063; -.
DR   Proteomes; UP000001940; Chromosome II.
DR   ExpressionAtlas; A0A0K3AWR5; baseline and differential.
PE   4: Predicted;
KW   Complete proteome {ECO:0000313|Proteomes:UP000001940};
KW   Reference proteome {ECO:0000313|Proteomes:UP000001940}.
SQ   SEQUENCE   220 AA;  26028 MW;  E12D5EA7F6FFF373 CRC64;
     MSDNEEVYVN FRGMNCISTG KSASMVPSKR RNWPKRVKKR LSTQRNNQKT IRPPELNKNN
     IEIKDMNSNN LEERNREECI QPVSVEKNIL HFEKFKSNQI CIVRENNKFR EGTRRRRKNS
     GESEDLKIHE NFTEKRRPIR SCKQNISFYE MDGDIEEFEV FFDTPTKSKK VLLDIYSAKK
     MPKIEVEDSL VNKFHSKRPS RACRVLGSME EVPFDVEIGY
//
EOF
test_parser("XrefParser::UniProtParser", $uniprot_elegans_record,  {
  xref => 3,
  primary_xref => 1,
  dependent_xref => 2,
},"Example UniProt record"); 
(my $uniprot_elegans_record_embl = $uniprot_elegans_record) =~ s/DR   EMBL;.*?\n//;
test_parser("XrefParser::UniProtParser", $uniprot_elegans_record_embl,  {
  xref => 1,
  primary_xref => 1,
},"EMBL entries are the dependent xrefs");
my @recognised_sources = (
 "PDB; 3HRI; X-ray; 2.85 A; A/B/C/D/E/F=44-477.",
 "MEROPS; C26.956; -.",
);
for my $l (@recognised_sources) {
  (my $uniprot_elegans_record_extra_line = $uniprot_elegans_record) =~ s/DR(.*?)\n/DR$1\nDR   $l/;
  test_parser("XrefParser::UniProtParser", $uniprot_elegans_record_extra_line,  {
    xref => 4,
    primary_xref => 1,
    dependent_xref => 3,
  }, "Pick up as extra xref + dependent xref: $l" );
}
sub test_elegans_uniprot {
  my ($expected_count, $extra_line) = @_;
  my $uniprot_elegans_record_here = $uniprot_elegans_record;
  $uniprot_elegans_record_here =~ s/DR(.*?)\n/DR$1\nDR  $extra_line/ if $extra_line;
  test_parser("XrefParser::WormbaseCElegansUniProtParser", 
      $uniprot_elegans_record_here,
     {}, 
     "No UniProt entries without corresponding INSDC entries $extra_line"
  );
  test_parser("XrefParser::WormbaseDirectParser",
    $wormbase_celegans_xrefs_head,
    $wormbase_celegans_xrefs_expected_count,
    "Test again to set up the next test",
    skip_clean => 1,
  );
  test_parser("XrefParser::WormbaseCElegansUniProtParser",
    $uniprot_elegans_record_here,
    $expected_count,
    "Correct xrefs and dependent xrefs $extra_line",
  );
}
my $wormbase_and_uniprot_expected_count = {
  %$wormbase_celegans_xrefs_expected_count,
  xref => $wormbase_celegans_xrefs_expected_count->{xref}+1,
  dependent_xref => 1 #protein id still there, no parent sequence ID
};
test_elegans_uniprot($wormbase_and_uniprot_expected_count, "");
for my $l (@recognised_sources) {
  test_elegans_uniprot({
    %$wormbase_and_uniprot_expected_count,
    xref => $wormbase_and_uniprot_expected_count->{xref}+1,
    dependent_xref => $wormbase_and_uniprot_expected_count->{dependent_xref}+1,
  },$l);
}
test_elegans_uniprot({
  %$wormbase_and_uniprot_expected_count,
  dependent_xref => $wormbase_and_uniprot_expected_count->{dependent_xref}+1,
  }, "EMBL; BX284602; CCD61130.1; -; Genomic_DNA",
);

my $refseq_protein_elegans_record = <<EOF;
LOCUS       NP_493629                427 aa            linear   INV 19-AUG-2018
DEFINITION  Uncharacterized protein CELE_2L52.1 [Caenorhabditis elegans].
ACCESSION   NP_493629
VERSION     NP_493629.2
DBLINK      BioProject: PRJNA158
            BioSample: SAMEA3138177
DBSOURCE    REFSEQ: accession NM_061228.2
KEYWORDS    RefSeq.
SOURCE      Caenorhabditis elegans
  ORGANISM  Caenorhabditis elegans
            Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida;
            Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis.
REFERENCE   
  <snipped>
COMMENT     REVIEWED REFSEQ: This record has been curated by WormBase. The
            reference sequence is identical to CCD61130.
FEATURES             Location/Qualifiers
     source          1..427
                     /organism="Caenorhabditis elegans"
                     /strain="Bristol N2"
                     /db_xref="taxon:$SPECIES_ID"
                     /chromosome="II"
     Protein         1..427
                     /product="hypothetical protein"
                     /calculated_mol_wt=49887
     CDS             1..427
                     /gene="2L52.1"
                     /locus_tag="CELE_2L52.1"
                     /standard_name="2L52.1a"
                     /coded_by="NM_061228.2:1..1284"
                     /note="Confirmed by transcript evidence"
                     /db_xref="EnsemblGenomes-Gn:WBGene00007063"
                     /db_xref="EnsemblGenomes-Tr:2L52.1a"
                     /db_xref="GeneID:181792"
                     /db_xref="GOA:A4F336"
                     /db_xref="InterPro:IPR013087"
                     /db_xref="UniProtKB/TrEMBL:A4F336"
                     /db_xref="WormBase:WBGene00007063"
ORIGIN      
        1 msmvrnvsnq sekleilsck wvgclkstev fktveklldh vtadhipevi vnddgseevv
       61 cqwdccemga srgnlqkkke wmenhfktrh vrkakifkcl iedcpvvkss sqeiethlri
      121 shpinpkker lkefksstdh ieptqanrvw tivngevqwk tpprvkkktv iyyddgpryv
      181 fptgcarcny dsdeselesd efwsatemsd neevyvnfrg mncistgksa smvpskrrnw
      241 pkrvkkrlst qrnnqktirp pelnknniei kdmnsnnlee rnreeciqpv sveknilhfe
      301 kfksnqiciv rennkfregt rrrrknsges edlkihenft ekrrpirsck qnisfyemdg
      361 dieefevffd tptkskkvll diysakkmpk ievedslvnk fhskrpsrac rvlgsmeevp
      421 fdveigy
//
EOF

my $refseq_mrna_elegans_record = <<EOF;
LOCUS       NM_001313558             663 bp    mRNA    linear   INV 19-AUG-2018
DEFINITION  Caenorhabditis elegans Uncharacterized protein (2L52.1), partial
            mRNA.
ACCESSION   NM_001313558
VERSION     NM_001313558.1
DBLINK      BioProject: PRJNA158
            BioSample: SAMEA3138177
KEYWORDS    RefSeq.
SOURCE      Caenorhabditis elegans
  ORGANISM  Caenorhabditis elegans
            Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida;
            Rhabditoidea; Rhabditidae; Peloderinae; Caenorhabditis.
REFERENCE   1  (bases 1 to 663)
  <snipped>
COMMENT     REVIEWED REFSEQ: This record has been curated by WormBase. This
            record is derived from an annotated genomic sequence (NC_003280).
            COMPLETENESS: incomplete on both ends.
FEATURES             Location/Qualifiers
     source          1..663
                     /organism="Caenorhabditis elegans"
                     /mol_type="mRNA"
                     /strain="Bristol N2"
                     /db_xref="taxon:$SPECIES_ID"
                     /chromosome="II"
     gene            <1..>663
                     /gene="2L52.1"
                     /locus_tag="CELE_2L52.1"
                     /db_xref="GeneID:181792"
                     /db_xref="WormBase:WBGene00007063"
     CDS             1..663
                     /gene="2L52.1"
                     /locus_tag="CELE_2L52.1"
                     /standard_name="2L52.1a"
                     /note="Confirmed by transcript evidence"
                     /codon_start=1
                     /product="hypothetical protein"
                     /protein_id="NP_001300487.1"
                     /db_xref="EnsemblGenomes-Gn:WBGene00007063"
                     /db_xref="EnsemblGenomes-Tr:2L52.1b"
                     /db_xref="GeneID:181792"
                     /db_xref="UniProtKB/TrEMBL:A0A0K3AWR5"
                     /db_xref="WormBase:WBGene00007063"
                     /translation="MSDNEEVYVNFRGMNCISTGKSASMVPSKRRNWPKRVKKRLSTQ
                     RNNQKTIRPPELNKNNIEIKDMNSNNLEERNREECIQPVSVEKNILHFEKFKSNQICI
                     VRENNKFREGTRRRRKNSGESEDLKIHENFTEKRRPIRSCKQNISFYEMDGDIEEFEV
                     FFDTPTKSKKVLLDIYSAKKMPKIEVEDSLVNKFHSKRPSRACRVLGSMEEVPFDVEI
                     GY"
ORIGIN      
        1 atgtcagata atgaagaagt atatgtgaac ttccgtggaa tgaactgtat ctcaacagga
       61 aagtcggcca gtatggtccc gagcaaacga agaaattggc caaaaagagt gaagaaaagg
      121 ctatcgacac aaagaaacaa tcagaaaact attcgaccac cagagctgaa taaaaataat
      181 atagagataa aagatatgaa ctcaaataac cttgaagaac gcaacagaga agaatgcatt
      241 cagcctgttt ctgttgaaaa gaacatcctg cattttgaaa aattcaaatc aaatcaaatt
      301 tgcattgttc gggaaaacaa taaatttaga gaaggaacga gaagacgcag aaagaattct
      361 ggtgaatcgg aagacttgaa aattcatgaa aactttactg aaaaacgaag acccattcga
      421 tcatgcaaac aaaatataag tttctatgaa atggacgggg atatagaaga atttgaagtg
      481 tttttcgata ctcccacaaa aagcaaaaaa gtacttctgg atatctacag tgcgaagaaa
      541 atgccaaaaa ttgaggttga agattcatta gttaataagt ttcattcaaa acgtccatca
      601 agagcatgtc gagttcttgg aagtatggaa gaagtaccat ttgatgtgga aataggatat
      661 tga
//
EOF
sub test_refseq {
  my ($record, $type) = @_;
  my $file = "RefSeq_${type}.gpff";
  test_parser("XrefParser::RefSeqGPFFParser",$record, {
    xref =>1,
    primary_xref => 1,
  }, "$type: Example record" , tmp_file_name => $file);
  test_parser("XrefParser::WormbaseCElegansRefSeqGPFFParser",$record, {
  }, "$type no entries without WormBase records" , tmp_file_name => $file);
  test_parser("XrefParser::WormbaseDirectParser", $wormbase_celegans_xrefs_head,
      $wormbase_celegans_xrefs_expected_count, "$type test again to set up the next test",
  skip_clean => 1);
  test_parser("XrefParser::WormbaseCElegansRefSeqGPFFParser",$record, {
    %$wormbase_celegans_xrefs_expected_count,
    xref => $wormbase_celegans_xrefs_expected_count->{xref}+1,
    dependent_xref => 1,
  }, "$type RefSeq entries hang off INSDC entries", tmp_file_name => $file);
}

test_refseq($refseq_protein_elegans_record, "protein");
test_refseq($refseq_mrna_elegans_record, "mrna");

done_testing();


