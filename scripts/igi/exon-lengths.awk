#! /usr/bin/nawk -f
# $Id$
#  Little script to find statistics on exon lengths from a gtf file.
BEGIN { min=100000; max=-1}
$3=="exon" {
  n++
  len=($5-$4)
  tot+=len; 
  if (len<min) min=len
  if (len>max) max=len
}
END {   
  print "avg: ",  tot/n, "min: ", min, "max: ", max
}
