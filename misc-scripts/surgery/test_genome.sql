# Copyright EMBL-EBI 2001
# Author: Alistair Rust
# Creation: 04.10.2001
# Last modified:
#
# File name: test_genome_populate.sql
#
# This file should be used in conjunction with 
# test_genome_table.sql which generates the empty
# database into which the following data is stored.
#
# The current implementation is to grab 2Mbs from:
# - 233M to 224M from chr 2
# - 1 to 1M from chr 20
#
# The sql also modifies the Rule* tables to run a subset
# of jobs
#

use arne_2MB_test;
               
# create a temp table to store clone ids for those clones
# that we're interested in

create  temporary table tmp1( 
        clone_id int(10) NOT NULL );


# find those clones present in some central, exciting
# 1Mbases on chromosome 2

insert  into tmp1
select  distinct( c.clone_id )
from    arne_ens130.contig c, arne_ens130.assembly s,
        arne_ens130.chromosome chr
where   s.contig_id = c.contig_id
and     s.chromosome_id = chr.chromosome_id
and     chr.name = '2'
and     chr_start < 224000000
and     chr_end > 223000001;


# some nice sticky exons

insert  into tmp1
select  distinct(c.clone_id)
from    arne_ens130.contig c, arne_ens130.assembly s,
        arne_ens130.chromosome chr
where   s.contig_id = c.contig_id
and     s.chromosome_id = chr.chromosome_id
and     chr.name = '11'
and     chr_end > 62800000
and     chr_start < 63800000;


# slice out the relevant clones from the arne_ens130
# clone table to create the test genome clone table

insert  into clone
select  h.* 
from    arne_ens130.clone h, tmp1 t
where   h.clone_id = t.clone_id;


# retrieve all contigs on the clones present in
# the first 1Mbases be they on the Golden Path
# or not

insert  into contig
select  h.contig_id, h.name, h.clone_id, h.length, h.offset, h.corder,
        h.dna_id, h.international_name
from    arne_ens130.contig h, tmp1 t
where   h.clone_id = t.clone_id;


# create the relevant dna table for the contigs
# in the test genome

insert  into dna
select  d.*
from    arne_ens130.dna d, contig
where   d.dna_id = contig.dna_id;



# copy the analysis table from the arne_ens130 db

insert  into analysis
select  *
from    arne_ens130.analysis;


#
# copy the static_golden_path table from arne_ens130 db
#

insert  into assembly
select  a.*
from    arne_ens130.assembly a, arne_ens130.chromosome chr
where a.chromosome_id = chr.chromosome_id
and chr.name = '2';

insert  into assembly
select  a.*
from    arne_ens130.assembly a, arne_ens130.chromosome chr
where a.chromosome_id = chr.chromosome_id
and chr.name = '20';

# copy overlapping features, repeats
# genes, exons, etc

# first find gene start end thing

create table gene_global_start_end 
SELECT STRAIGHT_JOIN tr.gene_id,
MIN(IF( a.contig_ori=1,
        (e.contig_start+a.chr_start-a.contig_start),
        (a.chr_start+a.contig_end-e.contig_end))) 
    as start,
MAX(IF( a.contig_ori=1,
        (e.contig_end+a.chr_start-a.contig_start),
        (a.chr_start+a.contig_end-e.contig_start))) 
    as end, 
IF (a.contig_ori=1,e.contig_strand,(-e.contig_strand)) as strand,
chr.name as chromosome

from   arne_ens130.transcript tr, arne_ens130.exon_transcript et,
       arne_ens130.exon e, arne_ens130.assembly a, arne_ens130.chromosome chr
where  tr.transcript_id = et.transcript_id
  and  et.exon_id = e.exon_id
  and  e.contig_id = a.contig_id
  and  a.chromosome_id = chr.chromosome_id
group by tr.gene_id;

# then make gene list in area
# genes have to be completely in

select "Gene start end done.";

create table gene_list
select gene_id from gene_global_start_end g
where g.chromosome = '11'
  and g.end < 63800000
  and g.start > 62800000;

insert into gene_list
select gene_id from gene_global_start_end g
where g.chromosome = '2'
  and g.start >= 223000000
  and g.end < 224000000;

alter table gene_list add index gene_id_idx( gene_id );

# now copy gene and gene_stable_id

insert into gene
select g.* 
from arne_ens130.gene g, gene_list gl
where g.gene_id = gl.gene_id;

insert into gene_stable_id
select gsi.*
from arne_ens130.gene_stable_id gsi, gene_list gl
where gsi.gene_id = gl.gene_id;

# now transcript, transcript_stable_id

insert into transcript
select tr.* 
from arne_ens130.transcript tr, gene_list gl
where tr.gene_id = gl.gene_id;

insert into transcript_stable_id
select tsi.*
from arne_ens130.transcript_stable_id tsi, transcript tr
where tsi.transcript_id = tr.transcript_id;

# translations, translation_stable_id

insert into translation
select tl.*
from arne_ens130.translation tl, transcript tr
where tr.translation_id = tl.translation_id;

insert into translation_stable_id
select tsi.*
from arne_ens130.translation_stable_id tsi, translation tl
where tsi.translation_id = tl.translation_id;

# exon_transcript

insert into exon_transcript
select et.*
from arne_ens130.exon_transcript et, transcript tr
where tr.transcript_id = et.transcript_id;

# exons are tricky, first make a unique list

create table exon_unique
select distinct exon_id
from exon_transcript;

# then uses this list to copy 

insert into exon
select e.*
from arne_ens130.exon e, exon_unique eu
where e.exon_id = eu.exon_id;

insert into exon_stable_id
select esi.*
from arne_ens130.exon_stable_id esi, exon_unique eu
where esi.exon_id = eu.exon_id;

drop table exon_unique;
drop table gene_list;
drop table gene_global_start_end;

# meta table

insert into meta
select * 
from arne_ens130.meta;


# object_xref, identity_xref, xref
# gene_description
# external_db

insert into gene_description
select gd.* 
from arne_ens130.gene_description gd, gene g
where gd.gene_id = g.gene_id;

insert into object_xref
select ox.* 
from arne_ens130.object_xref ox, translation tr
where ox.translation_id = tr.translation_id;

insert into xref
select x.*
from arne_ens130.xref x, object_xref ox
where x.xref_id = ox.xref_id;

insert into identity_xref 
select ix.* 
from arne_ens130.identity_xref ix, object_xref ox
where ix.object_xref_id = ox.object_xref_id;

insert into external_db
select ed.*
from arne_ens130.external_db;

# the interpro we need are those which are in xref







# finally, drop the temporary table

drop    table tmp1;