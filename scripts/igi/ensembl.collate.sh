#!/bin/sh -x
# -*- mode: sh; -*-

nawk -F\t '$3=="exon" || $3 ~ "_codon" '
# | sort -k1,1 -k7,7 -k4,4n 
