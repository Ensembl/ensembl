#!/bin/sh

mysql -u ensro -h ecs3 -P 3307 -e 'SELECT week_beginning,hits,size,pi FROM summary WHERE YEAR(week_beginning) >= 2003 AND vhost_id = 3 ORDER BY week_beginning' ensembl_summary_stats
