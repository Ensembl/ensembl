#CREATE UNIQUE INDEX stable_id_lookup_idx USING BTREE ON stable_id_lookup(stable_id,db_type,object_type);
CREATE INDEX stable_id_db_type USING BTREE ON stable_id_lookup(stable_id,db_type);
CREATE INDEX stable_id_object_type USING BTREE ON stable_id_lookup(stable_id,object_type);
