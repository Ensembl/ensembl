# Select all genes
select   * 
from     gene;


# Select a specific transcript
select   *
from     transcript
where    id =  'dJ402G11.05009.GENSCAN.0.CDS.1';

# Select all + strand exons (id and contig only)

select   id,contig
from     exon
where    strand = '+';

# Select all exons in gene   dJ40E16.02584.GENSCAN.0
select   exon.id
from     exon,transcript,exon_bridge
where    exon_bridge.exon = exon.id                     AND
         exon_bridge.transcript = transcript.id         AND
         transcript.geneid = 'dJ40E16.02584.GENSCAN.0';
         
# how many genscan predictions are there?
select   homolset,count(*) 
from     homol
where    homolset != '' 
group by homolset;

# how many times does a blastp match hit a contig
select    id,count(*) 
from      homol 
where     analysis = 'blastp_swir' 
group by  id; 

# which database hits hit different contigs - this prints ALL links
# this should have a group by contig here somewhere
select    p1.contig,p2.contig,p2.id 
from      homol as p1, homol as p2 
where     p1.contig != p2.contig AND 
          p1.id = p2.id;

# this might be better
# just prints out different contigs
select    count(*),homol.id,homol.contig 
from      homol 
where     id = 'WP:F15G9.4A' 
group by  homol.contig;

# groups by contigs for blastp hits for all homols
select    count(*),homol.id,homol.contig 
from      homol 
where     analysis = 'blastp_swir'
group by  homol.contig;

# similar to above but counting the number of contig-contig links
select    p1.id,count(*) 
from      homol as p1, homol as p2 
where     p1.contig != p2.contig AND 
          p1.id = p2.id 
group by  p1.id;

# show all hits between different contigs for id XL28SR
# n.b. brings up hits twice - don't know how to get around this.
select    p1.* 
from      homol as p1, homol as p2 
where     p1.contig != p2.contig AND 
          p1.id = p2.id          AND 
          p1.id = 'XL28SR';

# count number of exons for gene  dJ437I16.00294.GENSCAN.5
select    transcript.geneid,count(*) 
from      exon_bridge,transcript 
where     transcript.geneid = 'dJ437I16.00294.GENSCAN.5' AND  
          exon_bridge.transcript = transcript.id 
group by  transcript.geneid;

# count number of exons for all genes
select    transcript.geneid,count(*)
from      exon_bridge,transcript 
where     exon_bridge.transcript = transcript.id
group by  transcript.geneid;


# find homol overlaps with exon dJ437I16.00294.GENSCAN.0.1
select    homol.* 
from      homol,exon 
where     homol.analysis != ''                   AND
          exon.id = 'dJ437I16.00294.GENSCAN.0.1' AND 
          homol.strand = exon.strand             AND 
          homol.start >= exon.start              AND 
          homol.end   <= exon.end                AND 
          homol.contig = exon.contig;

# or to just report the different homol ids and number of hits
select    homol.id,count(*)
from      homol,exon 
where     homol.analysis != ''                   AND
          exon.id = 'dJ437I16.00294.GENSCAN.0.1' AND 
          homol.strand = exon.strand             AND 
          homol.start >= exon.start              AND 
          homol.end   <= exon.end                AND 
          homol.contig = exon.contig
group by  homol.id;


# find all homol hits with all exons
# note - can't find how to order these by exon id
select    exon.id,homol.id,count(*)
from      homol,exon 
where     homol.analysis != ''                   AND
          homol.strand = exon.strand             AND 
          homol.start >= exon.start              AND 
          homol.end   <= exon.end                AND 
          homol.contig = exon.contig
group by  homol.id;


# find all homol hits that overlap with genscan predictions
select    p1.id,p2.id,count(*) 
from      homol as p1, homol as p2 
where     p1.contig = p2.contig      AND
          p1.analysis != ''          AND
	  p2.analysis  = ''          AND
          p2.strand = p1.strand      AND 
          p1.start >= p2.start       AND 
          p1.end <= p2.end                     
group by  p1.id;
