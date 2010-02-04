# patch_54_55_d.sql
#
# title: Add table to store the dependent xrefs
#
# description:
# For a given object_xref store its master if it is dependent on another xref (master)

CREATE TABLE dependent_xref(
     object_xref_id         INT NOT NULL,
     master_xref_id         INT NOT NULL,
     dependent_xref_id      INT NOT NULL,

     PRIMARY KEY( object_xref_id ),
     KEY dependent ( dependent_xref_id ),
     KEY master_idx (master_xref_id)

) COLLATE=latin1_swedish_ci TYPE=MyISAM;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_d.sql|add_dependent_xref_table');

