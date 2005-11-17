# There were no schema changes on the core databases between release 34 and 35
# 
# Just update the schema_version in the meta table
UPDATE meta set meta_value="35" where meta_key="schema_version";
