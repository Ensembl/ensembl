#!/bin/sh -x
# -*- mode: sh; -*-

nawk -F\t '$3=="exon" || $3 ~ "_codon" ' 
