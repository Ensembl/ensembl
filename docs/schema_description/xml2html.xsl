<?xml version="1.0"?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

  <xsl:output method="html"/>

  <xsl:template match="/">

    <html>      
    <head>
      <title>EnsEMBL Core Schema Documentation</title>    
      <link rel="stylesheet" href="http://www.ensembl.org/EnsEMBL.css"/>
    </head>
    <body>
 
    <h1>EnsEMBL Core Schema Documentation</h1>

    <xsl:call-template name="intro"/>
    <xsl:call-template name="index"/>
    <xsl:apply-templates/>
    <xsl:call-template name="concepts"/>
    
    </body>
    </html>

  </xsl:template>


  <!-- Introduction -->
  <xsl:template name="intro" match="schemadescription/introduction">
    <h2>Introduction</h2>
    <p><xsl:value-of select="schemadescription/introduction/text"/></p>

     <p><xsl:value-of select="schemadescription/introduction/process/@intro"/></p>
     <p>
     <table border="1" cellpadding="10">
     <tr><th>Step</th><th>Process</th></tr>
     <xsl:for-each select="schemadescription/introduction/process/step">
       <tr><td><xsl:value-of select="@number"/></td><td><xsl:value-of select="."/></td></tr>
     </xsl:for-each>
     </table>
     </p>
     <!-- The replace() function requires an XPath2.0 compatible processor e.g. Saxon -->
    <!-- The ? in the middle of Revision is to prevent CVS keyword substitution happening on /this/ file! -->

    <p>This document refers to version <strong><xsl:value-of select="schemadescription/@schema-version"/></strong> of the EnsEMBL core schema. You are looking at revision <strong><xsl:value-of select='replace(schemadescription/@document-version, "\$Revisi?on:\s*(\d+\.\d+)\s*\$", "$1")'/></strong> of this document.</p>

    <xsl:if test="count(schemadescription/diagram) &gt; 0">
      <p><b>Diagrams:</b></p>
      <ul>
      <xsl:for-each select="schemadescription/diagram">
        <li> <a href="{@url}"><xsl:value-of select="@url"/> </a> - <xsl:value-of select="@description"/> </li>
      </xsl:for-each> 
      </ul>
    </xsl:if>
    <hr/>
  </xsl:template>

  <!-- Index -->
  <xsl:template name="index">
  <p>Quick links to tables:</p>
  <xsl:for-each select="schemadescription/tablegroup">
    <p><b><xsl:value-of select="@name"/></b> 
    <xsl:if test="string-length(@description) &gt; 0"> 
      (<xsl:value-of select="@description"/>) 
    </xsl:if>
    </p>
    <ul>
    <xsl:for-each select="table">
      <xsl:sort select="name"/>
      <li><a href="#{name}"><xsl:value-of select="name"/></a></li>
    </xsl:for-each>
    </ul>
  </xsl:for-each>
  </xsl:template>


  <xsl:template match="schemadescription">

    <xsl:apply-templates select="tablegroup"/>

  </xsl:template>

  <!-- Table groups -->
  <xsl:template match="schemadescription/tablegroup">

    <hr/>
    <h2> <xsl:value-of select="@name"/> </h2>
    <p> <xsl:value-of select="@description"/> </p>
    <xsl:apply-templates/>

  </xsl:template>

  <!-- Table info -->
  <xsl:template match="tablegroup/table">
   
    <h3><a name="{name}"><xsl:value-of select="name"/></a></h3>

    <p><xsl:value-of select="description"/></p>

    <!-- see also -->
    <xsl:if test="count(see/*) &gt; 0">
      <p><b>See also:</b></p>
      <xsl:if test="count(see/tableref) &gt; 0">
        <p>Tables:</p>
        <ul>
        <xsl:for-each select="see/tableref">
          <li> <a href="#{@name}"><xsl:value-of select="@name"/> </a> - <xsl:value-of select="@reason"/> </li>
        </xsl:for-each> 
        </ul>
      </xsl:if>
      <xsl:if test="count(see/conceptref) &gt; 0">
        <p>Concepts:</p>
        <ul>
        <xsl:for-each select="see/conceptref">
          <li> <a href="#{@name}"><xsl:value-of select="@name"/> </a> - <xsl:value-of select="@reason"/> </li>
        </xsl:for-each> 
        </ul>
      </xsl:if>
      <xsl:if test="count(see/urlref) &gt; 0">
        <p>Links:</p>
        <ul>
        <xsl:for-each select="see/urlref">
          <li> <a href="{@url}"><xsl:value-of select="@name"/> </a> - <xsl:value-of select="@reason"/> </li>
        </xsl:for-each> 
        </ul>
      </xsl:if>

    </xsl:if>

    <!-- TODO: used -->
  
  </xsl:template>

  <!-- Concepts -->
  <xsl:template name="concepts" match="schemadescription/concepts">

    <hr/>
    <h2> Concepts </h2>
      <dl>
        <xsl:for-each select="schemadescription/concepts/concept">
          <dt> <p> <strong> <a name="{@name}"><xsl:value-of select="@name"/> </a> </strong> </p> </dt>
          <dd> <p> <xsl:value-of select="@description"/> </p> </dd>
        </xsl:for-each> 
      </dl>
    <hr/>
  </xsl:template>

</xsl:stylesheet>

