#!/bin/sh

java -jar /scratch/saxon/saxon7.jar tables.xml xml2html.xsl > tables.html 

java -jar /scratch/saxon/saxon7.jar tables.xml xml2wiki.xsl > tables.txt