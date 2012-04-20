-- Patching in support for comments on web data

ALTER TABLE web_data
ADD COLUMN comment TEXT 
AFTER data;
