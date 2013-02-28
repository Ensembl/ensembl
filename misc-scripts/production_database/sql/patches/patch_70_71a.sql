-- Patching in support for comments on web data

ALTER TABLE analysis_web_data
MODIFY COLUMN db_type enum('cdna','core','funcgen','otherfeatures','rnaseq','vega','presite','sangervega') NOT NULL DEFAULT 'core';
