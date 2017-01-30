DROP INDEX stable_id_db_type ON stable_id_lookup;
DROP INDEX stable_id_object_type ON stable_id_lookup;
CREATE INDEX stable_id_db_type USING BTREE ON stable_id_lookup(stable_id,db_type);
CREATE INDEX stable_id_object_type USING BTREE ON stable_id_lookup(stable_id,object_type);
