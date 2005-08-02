# SQL to patch a release 32 Ensembl database schema to release 33

# Add db_display_name column to external_db
ALTER TABLE external_db ADD COLUMN db_display_name VARCHAR(255);
