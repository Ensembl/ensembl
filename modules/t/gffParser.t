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
use Test::Warnings;
use IO::String;
use Bio::EnsEMBL::Utils::IO::GFFParser;

{
	my $io = IO::String->new(<<'GFF');
##gff-version 3
##sequence-region   ctg123 1 1497228
##taken-from http://www.sequenceontology.org/gff3.shtml
ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN

ctg123	.	TF_binding_site	1000	1012	.	+	.	ID=tfbs00001;Parent=gene00001

ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123	.	mRNA	1300	9000	.	+	.	ID=mRNA00003;Parent=gene00001;Name=EDEN.3

ctg123	.	exon	1300	1500	.	+	.	ID=exon00001;Parent=mRNA00003
ctg123	.	exon	1050	1500	.	+	.	ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123	.	exon	3000	3902	.	+	.	ID=exon00003;Parent=mRNA00001,mRNA00003
ctg123	.	exon	5000	5500	.	+	.	ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123	.	exon	7000	9000	.	+	.	ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003

ctg123	.	CDS	1201	1500	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	3000	3902	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	5000	5500	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	7000	7600	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
GFF

	my $gff = Bio::EnsEMBL::Utils::IO::GFFParser->new($io);
	my $header = $gff->parse_header();
	is_deeply(
		$header,
		[ '##gff-version 3', '##sequence-region   ctg123 1 1497228', '##taken-from http://www.sequenceontology.org/gff3.shtml'],
		'Checking headers all parse'
	);

	my $actual_gene = $gff->parse_next_feature();
	my $expected_gene = {
		seqid => 'ctg123', start => 1000, end => 9000, strand => 1,
		source => '.', type => 'gene', score => '.', phase => '.',
		attribute => { ID => 'gene00001', Name => 'EDEN' }
	};
	is_deeply($actual_gene, $expected_gene, 'Checking gene record parses');
	
	my $actual_tf = $gff->parse_next_feature();
  my $expected_tf = {
    seqid => 'ctg123', start => 1000, end => 1012, strand => 1,
    source => '.', type => 'TF_binding_site', score => '.', phase => '.',
    attribute => { ID => 'tfbs00001', Parent => 'gene00001' }
  };
  is_deeply($actual_tf, $expected_tf, 'Checking TF record parses');
  
  #SKIP TO EXONS
  $gff->parse_next_feature(); #mrna
  $gff->parse_next_feature(); #mrna
  $gff->parse_next_feature(); #mrna
  
  #EXONS
  {
    my $actual = $gff->parse_next_feature();
    my $expected = {
    seqid => 'ctg123', start => 1300, end => 1500, strand => 1,
    source => '.', type => 'exon', score => '.', phase => '.',
    attribute => { ID => 'exon00001', Parent => 'mRNA00003' }
    };
    is_deeply($actual, $expected, 'Checking Exon 1 record parses');
  }
  {
    my $actual = $gff->parse_next_feature();
    my $expected = {
    seqid => 'ctg123', start => 1050, end => 1500, strand => 1,
    source => '.', type => 'exon', score => '.', phase => '.',
    attribute => { ID => 'exon00002', Parent => ['mRNA00001', 'mRNA00002'] }
    };
    is_deeply($actual, $expected, 'Checking Exon 2 record parses');
  }
}

{
  my $string = <<GFF;
##gff-version 3
##sequence-region   ctg123 1 1497228
##taken-from example parsing file with FASTA directive
ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN
###
##FASTA
>ctg123
AACCTTTGGGCCGGGCCTTAAAA
AACC
GFF
  my $io = IO::String->new(\$string);
  
  my $gff = Bio::EnsEMBL::Utils::IO::GFFParser->new($io);
  my $header = $gff->parse_header();
  is_deeply(
    $header,
    [ '##gff-version 3', '##sequence-region   ctg123 1 1497228', '##taken-from example parsing file with FASTA directive'],
    'Checking headers all parse'
  );
  my $actual_gene = $gff->parse_next_feature();
  my $expected_gene = {
    seqid => 'ctg123', start => 1000, end => 9000, strand => 1,
    source => '.', type => 'gene', score => '.', phase => '.',
    attribute => { ID => 'gene00001', Name => 'EDEN' }
  };
  is_deeply($actual_gene, $expected_gene, 'Checking TF record parses');
  my $id = $gff->parse_next_feature();
  ok(! defined $id, 'No more features');
  my $seq = $gff->parse_next_sequence();
  is_deeply($seq, {header => '>ctg123', sequence => "AACCTTTGGGCCGGGCCTTAAAAAACC"}, "Checking Sequence parses correctly")
}

done_testing();
