/*
 *    AltSplicingToolkit
 *    Author: Gautier Koscielny <koscieln@ebi.ac.uk>
 *
 *    Copyright (c) 1999-2012 The European Bioinformatics Institute and 
 *    Genome Research Limited, and others.  All rights reserved. 
 *
 *    Redistribution and use in source and binary forms, with or without 
 *    modification, are permitted provided that the following conditions 
 *    are met: 
 *
 *    1. Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer. 
 *
 *    2. Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer 
 *       in the documentation and/or other materials provided with the 
 *       distribution. 
 *
 *    3. The name "Ensembl" must not be used to endorse or promote products 
 *       derived from this software without prior written permission.  For 
 *       written permission, please contact helpdesk@ensembl.org 
 *
 *    4. Products derived from this software may not be called "Ensembl" 
 *       nor may "Ensembl" appear in their names without prior written 
 *       permission of the Ensembl developers. 
 *
 *    5. Redistributions in any form whatsoever must retain the following 
 *       acknowledgement: 
 *
 *         "This product includes software developed by Ensembl 
 *         (http://www.ensembl.org/)" 
 *
 *    THIS SOFTWARE IS PROVIDED BY THE ENSEMBL GROUP ``AS IS'' AND ANY 
 *    EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
 *    PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE ENSEMBL GROUP OR ITS 
 *    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
 *    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
 *    OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 *    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 *    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 *    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 */

#ifndef GFFPARSER_H_
#define GFFPARSER_H_

#include <string>
#include <fstream>
#include <iostream>


using namespace std;

namespace gff {

  const int GFF_CHR=0;
  const int GFF_DATASOURCE=1;
  const int GFF_TYPE=2;
  const int GFF_START=3;
  const int GFF_END=4;
  const int GFF_SCORE=5;
  const int GFF_STRAND=6;
  const int GFF_PHASE=7;
  const int GFF_COMMENTS=8;

  const int GFF_GENE_ID=0;
  const int GFF_TRANSCRIPT_ID=1;
  const int GFF_EXON_ID=2;

  class GffHandler {

  public:
    virtual void start() = 0;
    virtual bool newline(string & str) = 0;
    virtual void end() = 0;

    // copy constructor
    const GffHandler &operator=( const GffHandler& inRhs ) { return inRhs; };
    // assignment operator
    //void GffHandler( const GffHandler &inParam );

  };

  class GffParser {

  public:
    GffParser(istream &iStream, GffHandler *handler);
    virtual ~GffParser();

  public:
    void parse();

  protected:
    istream &iStream;
    GffHandler *pHandler; //OK, pfile should point to a derived object

  };



}

#endif /* GFFPARSER_H_ */
