use ensembl;

select * from analysis into outfile '/tmp/analysis.table';
select * from clone into outfile '/tmp/clone.table';
select * from contig into outfile '/tmp/contig.table';
select * from contig_equiv into outfile '/tmp/contig_equiv.table';
select * from db_update into outfile '/tmp/db_update.table';
select * from dna into outfile '/tmp/dna.table';
select * from exon into outfile '/tmp/exon.table';
select * from exon_transcript into outfile '/tmp/exon_transcript.table';
select * from feature into outfile '/tmp/feature.table';
select * from fset into outfile '/tmp/fset.table';
select * from fset_feature into outfile '/tmp/fset_feature.table';
select * from gene into outfile '/tmp/gene.table';
select * from geneclone_neighbourhood into outfile '/tmp/geneclone_neighbourhood.table';
select * from ghost into outfile '/tmp/ghost.table';
select * from mapbin into outfile '/tmp/mapbin.table';
select * from meta into outfile '/tmp/meta.table';
select * from repeat_feature into outfile '/tmp/repeat_feature.table';
select * from supporting_feature into outfile '/tmp/supporting_feature.table';
select * from transcript into outfile '/tmp/transcript.table';
select * from translation into outfile '/tmp/translation.table';


