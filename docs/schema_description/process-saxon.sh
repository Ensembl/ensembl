#!/bin/sh

java -jar "/nfs/acari/gp1/saxon/saxon8.jar" tables.xml xml2html.xsl > schema_description.html 

java -jar "/nfs/acari/gp1/saxon/saxon8.jar" tables.xml xml2wiki.xsl > schema_description.txt

cp schema_description.html ../../../private-web/htdocs/Docs/

cd ../../../private-web/htdocs/Docs/

cvs commit schema_description.html