<?xml version="1.0"?>

<!-- XSLT for creating Wiki-format documentation from tables.xml -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:output method="text"/>

<xsl:template match="/">

=EnsEMBL Core Schema Documentation=

<xsl:call-template name="intro"/>
<xsl:call-template name="index"/>
<xsl:apply-templates/>
<xsl:call-template name="concepts"/>

</xsl:template>


<!-- Introduction -->
<xsl:template name="intro" match="schemadescription/introduction">
==Introduction==

<xsl:value-of select="schemadescription/introduction/text"/>

<xsl:value-of select="schemadescription/introduction/process/@intro"/>

||Step||Process||
<xsl:for-each select="schemadescription/introduction/process/step">||<xsl:value-of select="@number"/>||<xsl:value-of select="."/>||
</xsl:for-each>
This document refers to version &apos;&apos;<xsl:value-of select="schemadescription/@schema-version"/>&apos;&apos; of the EnsEMBL core schema. You are looking at revision &apos;&apos;<xsl:value-of select='replace(schemadescription/@document-version, "\$Revisi?on:\s*(\d+\.\d+)\s*\$", "$1")'/>&apos;&apos;of this document.
<xsl:if test="count(schemadescription/diagram) &gt; 0">

&apos;&apos;&apos;Diagrams:&apos;&apos;&apos;
<xsl:for-each select="schemadescription/diagram">
 * <a href="{@url}"><xsl:value-of select="@url"/> </a> - <xsl:value-of select="@description"/>
</xsl:for-each> 
</xsl:if>
----
</xsl:template>
  
<!-- Index -->
<xsl:template name="index">
&apos;&apos;Quick links to tables:&apos;&apos;
<xsl:for-each select="schemadescription/tablegroup">

&apos;&apos;&apos;<xsl:value-of select="@name"/>&apos;&apos;&apos;
<xsl:if test="string-length(@description) &gt; 0"> 
(<xsl:value-of select="@description"/>) 
</xsl:if>
  
<xsl:for-each select="table">
<xsl:sort select="name"/>
 * <xsl:value-of select="name"/>
</xsl:for-each>
  </xsl:for-each>
</xsl:template>
<xsl:template match="schemadescription">
  <xsl:apply-templates select="tablegroup"/>
</xsl:template>
<!-- Table groups -->
<xsl:template match="schemadescription/tablegroup">
  ----
== <xsl:value-of select="@name"/> ==
<xsl:value-of select="@description"/> 
<xsl:apply-templates/>
</xsl:template>
<!-- Table info -->
<xsl:template match="tablegroup/table">
===<xsl:value-of select="name"/>===
<xsl:value-of select="description"/>
  <!-- see also -->
<xsl:if test="count(see/*) &gt; 0">
&apos;&apos;&apos;See also:&apos;&apos;&apos;
<xsl:if test="count(see/tableref) &gt; 0">

Tables:
<xsl:for-each select="see/tableref">
 * <a href="#{@name}"><xsl:value-of select="@name"/> </a> - <xsl:value-of select="@reason"/> 
</xsl:for-each> 
</xsl:if>
<xsl:if test="count(see/conceptref) &gt; 0">

Concepts:
<xsl:for-each select="see/conceptref">
 * <xsl:value-of select="@name"/>  - <xsl:value-of select="@reason"/> 
</xsl:for-each> 
</xsl:if>
<xsl:if test="count(see/urlref) &gt; 0">
      
Links:
<xsl:for-each select="see/urlref">
 * <xsl:value-of select="@name"/> - <xsl:value-of select="@reason"/> 
</xsl:for-each> 
</xsl:if>
</xsl:if>
<!-- TODO: used -->
</xsl:template>
<!-- Concepts -->
<xsl:template name="concepts" match="schemadescription/concepts">
 ----
==Concepts==
<xsl:for-each select="schemadescription/concepts/concept">
&apos;&apos;<xsl:value-of select="@name"/> &apos;&apos;
: <xsl:value-of select="@description"/> 
</xsl:for-each> 
  ----
</xsl:template>
</xsl:stylesheet>

