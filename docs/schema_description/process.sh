#!/bin/sh

java org.apache.xalan.xslt.Process -in tables.xml -xsl schema-description.xsl -out tables.html

