#!/usr/bin/sh

mysql -u ensro -h ecs3 -P 3307 -e 'select week_beginning,hits,size,pi from summary where (YEAR(week_beginning) = 2007 OR YEAR(week_beginning) = 2006 OR YEAR(week_beginning) = 2003 OR YEAR(week_beginning) = 2004 OR YEAR(week_beginning) = 2005)   and vhost_id = 3 order by week_beginning' ensembl_summary_stats
