# $Id$
# create a new description table that merges in the descriptions
# from the family database

# this script assumes that the family and ensembl-core database live in the
# same server (so you can do cross joins)

# The local table is also called gene_description

select count(*)
from homo_sapiens_core_110.gene_description gd
where gd.description is not null
  and gd.description not in ('', 'unknown', 'UNKNOWN');

# creating new local gene_description table;
CREATE TABLE gene_description AS 
SELECT *
FROM  homo_sapiens_core_110.gene_description gd 
WHERE gd.gene_id IS NULL;

ALTER TABLE gene_description ADD PRIMARY KEY(gene_id);
# inserting existing desc's:
INSERT INTO gene_description 
  SELECT *
  FROM homo_sapiens_core_110.gene_description;

# delete anything that looks remotely unknown:
delete from gene_description where description is null;
delete from gene_description where description = '';
delete from gene_description where description = 'unknown';
delete from gene_description where description = 'UNKNOWN';

# selecting all genes from gene, this time properly as 'unknown':
INSERT INTO gene_description
  SELECT id, 'unknown'
  FROM homo_sapiens_core_110.gene;
# this will fail on the knowns, as they should; only the 'unknowns' do get
# in.

# give stats:
select count(*)
from gene_description
where description = 'unknown';
  
create temporary table new_descriptions as
  select gd.gene_id, f.description
  from family f, family_members fm, gene_description gd
  where gd.description ='unknown'
   and fm.db_name ='ENSEMBLGENE'
   and fm.db_id = gd.gene_id
   and fm.family = f.internal_id
   and f.description <> 'UNKNOWN';

delete from gene_description where description = 'unknown';

insert into gene_description
  select * from new_descriptions;

# give new stats:
select count(*)
from gene_description;

